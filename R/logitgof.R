logitgof <-
function (obs, exp, g = 10, ord = FALSE) {
  DNAME <- paste(deparse(substitute(obs)), deparse(substitute(exp)), sep = ", ")
  yhat <- exp
  if (is.null(ncol(yhat))) {
    mult <- FALSE
  } else {
    if (ncol(yhat) == 1) {
      mult <- FALSE
    } else mult <- TRUE
  }
  n <- ncol(yhat)
  if (mult) {
    if (!ord) {
      METHOD <- "Hosmer and Lemeshow test (multinomial model)"
    } else {
      METHOD <- "Hosmer and Lemeshow test (ordinal model)"
    }
    qq <- unique(quantile(1 - yhat[, 1], probs = seq(0, 1, 1/g)))
    cutyhats <- cut(1 - yhat[, 1], breaks = qq, include.lowest = TRUE)
    dfobs <- data.frame(obs, cutyhats)
    dfobsmelt <- melt(dfobs, id.vars = 2)
    observed <- cast(dfobsmelt, cutyhats ~ value, length)
    observed <- observed[order(c(1, names(observed[, 2:ncol(observed)])))]
    dfexp <- data.frame(yhat, cutyhats)
    dfexpmelt <- melt(dfexp, id.vars = ncol(dfexp))
    expected <- cast(dfexpmelt, cutyhats ~ variable, sum)
    expected <- expected[order(c(1, names(expected[, 2:ncol(expected)])))]
    stddiffs <- abs(observed[, 2:ncol(observed)] - expected[, 2:ncol(expected)]) / sqrt(expected[, 2:ncol(expected)])
    if (ncol(expected) != ncol(observed)) stop("Observed and expected tables have different number of columns. Check you entered the correct data.")
    if (any(expected[, 2:ncol(expected)] < 1))
      warning("At least one cell in the expected frequencies table is < 1. Chi-square approximation may be incorrect.")
    if (any(observed[, 2:ncol(observed)] < 0)  || any(expected[, 2:ncol(expected)] < 0))
      stop("Observed or expected values < 0. Check that observed and fitted values entered correctly.")
    chisq <- sum((observed[, 2:ncol(observed)] - expected[, 2:ncol(expected)])^2 / expected[, 2:ncol(expected)])
    if (!ord) {
        PARAMETER <- (nrow(expected) - 2) * (ncol(yhat) - 1) 
      } else {
        PARAMETER <- (nrow(expected) - 2) * (ncol(yhat) - 1) + ncol(yhat) - 2
      }
    ### need to add a check that number of groups >= number of covariate patterns
  } else {
    METHOD <- "Hosmer and Lemeshow test (binary model)"
    if (is.factor(obs)) {
      y <- as.numeric(obs) - 1
    } else {
      y <- obs
    }
    qq <- unique(quantile(yhat, probs = seq(0, 1, 1/g)))
    cutyhat <- cut(yhat, breaks = qq, include.lowest = TRUE)
    observed <- xtabs(cbind(y0 = 1 - y, y1 = y) ~ cutyhat)
    expected <- xtabs(cbind(yhat0 = 1 - yhat, yhat1 = yhat) ~ cutyhat)
    stddiffs <- abs(observed - expected) / sqrt(expected)
    if (any(expected < 1))
      warning("At least one cell in the expected frequencies table is < 1. Chi-square approximation may be incorrect.")
    if (any(observed < 0) || any(expected < 0))
      stop("Observed or expected values < 0. Check that observed and fitted values entered correctly.")
    chisq <- sum((observed - expected)^2/expected)
    PARAMETER <- nrow(expected) - 2
  }
  if (g != nrow(expected))
    warning(paste("Not possible to compute", g, "rows. There might be too few observations."))
  PVAL <- 1 - pchisq(chisq, PARAMETER)
  names(chisq) <- "X-squared"
  names(PARAMETER) <- "df"
  structure(list(statistic = chisq, parameter = PARAMETER, 
                 p.value = PVAL, method = METHOD, data.name = DNAME, observed = observed, 
                 expected = expected, stddiffs = stddiffs), class = "htest")
}
