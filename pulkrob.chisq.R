pulkrob.chisq <-
function(model, catvars) {
  if (class(model) == "polr") {
    yhat <- as.data.frame(fitted(model))
  } else if (class(model) == "clm") {
    predprob <- model$model[, 2:ncol(model$model)]
    yhat <- predict(model, newdata = predprob, type = "prob")$fit
  } else warning("Model is not of class polr or clm. Test may fail.")
  formula <- formula(model$terms)
  DNAME <- paste("formula: ", deparse(formula))
  METHOD <- "Pulkstenis-Robinson chi-squared test"
  covars <- model$model[-1]
  covars <- covars[names(covars)==catvars]
  covpat <- epi.cp(covars)
  yhat$score <- apply(sapply(1:ncol(yhat), function(i) { yhat[, i] * i }), 1, sum)
  yhat <- cbind(id=1:nrow(yhat), yhat, covpat=covpat$id)
  medians <- cbind(covpat=covpat$cov.pattern$id,
                   med=sapply(covpat$cov.pattern$id, function(x) median(yhat[yhat$covpat==x, ]$score)))
  yhat <- merge(x = yhat, y = medians, by = "covpat", all.x = TRUE)
  yhat$covpatsplit <- sapply(1:nrow(yhat), function(i) ifelse(yhat[i, ]$score <= yhat[i, ]$med,
                                                              paste0(yhat[i, ]$covpat, "a"),
                                                              paste0(yhat[i, ]$covpat, "b")))
  dfobs <- cbind(id=1:nrow(model$model), model$model[1])
  dfobs <- merge(x = dfobs, y = yhat[, c("id", "covpatsplit")], by = "id", all.x = TRUE)
  dfobsmelt <- melt(dfobs[, -1], id.vars = "covpatsplit")
  observed <- cast(dfobsmelt, covpatsplit ~ value, length)
  observed <- observed[order(c(1, names(observed[, 2:ncol(observed)])))]
  dfexp <- yhat[, !colnames(yhat) %in% c("id", "covpat", "med", "score")]
  dfexpmelt <- melt(dfexp, id.vars = ncol(dfexp))
  expected <- cast(dfexpmelt, covpatsplit ~ variable, sum)
  expected <- expected[order(c(1, names(expected[, 2:ncol(expected)])))]
  stddiffs <- abs(observed[, 2:ncol(observed)] - expected[, 2:ncol(expected)]) / sqrt(expected[, 2:ncol(expected)])
  if (any(expected[, 2:ncol(expected)] < 1))
    warning("At least one cell in the expected frequencies table is < 1. Chi-square approximation may be incorrect.")
  chisq <- sum((observed[, 2:ncol(observed)] - expected[, 2:ncol(expected)])^2 / expected[, 2:ncol(expected)])
  I2 <- nrow(observed)
  J <- length(levels(as.factor(model$model[, 1])))
  k <- length(catvars)
  PARAMETER <- ((I2 - 1) * (J - 1)) - k - 1
  PVAL <- 1 - pchisq(chisq, PARAMETER)
  names(chisq) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = chisq, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = observed, 
                 expected = expected, stddiffs = stddiffs), class = "htest")
}
