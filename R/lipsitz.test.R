lipsitz.test <-
  function (model, g = 10)  {
    oldmodel <- model
    if (class(oldmodel) == "polr") {
      yhat <- as.data.frame(fitted(oldmodel))
    } else if (class(oldmodel) == "clm") {
      predprob <- oldmodel$model[, 2:ncol(oldmodel$model), drop = F]
      yhat <- as.data.frame(predict(oldmodel, newdata = predprob, type = "prob")$fit)
    } else warning("Model is not of class polr or clm. Test may fail.")
    formula <- formula(oldmodel$terms)
    DNAME <- paste("formula: ", deparse(formula))
    METHOD <- "Lipsitz goodness of fit test for ordinal response models"
    obs <- oldmodel$model[1]
    if (g < 6) warning("g < 6. Running this test when g < 6 is not recommended.")
    if (g >= nrow(obs) / (5 * ncol(yhat))) warning("g >= n/5c. Running this test when g >= n/5c is not recommended.")
    yhat$score <- apply(sapply(1:ncol(yhat), function(i) {
      yhat[, i] * i
    }), 1, sum)
    yhat$tmp <- 1:nrow(yhat)
    yhat <- yhat[order(yhat$score), ]
    cutyhats <- cut(1:nrow(yhat), breaks = g, include.lowest = T)
    cutyhats <- cutyhats[order(yhat$tmp)]
    yhat <- yhat[order(yhat$tmp), ]
    yhat$score <- NULL
    yhat$tmp <- NULL
    dfobs <- data.frame(obs, cutyhats)
    dfobsmelt <- melt(dfobs, id.vars = 2)
    observed <- cast(dfobsmelt, cutyhats ~ value, length)
    if (g != nrow(observed)) {
      warning(paste("Not possible to compute", g, "rows. There might be too few observations."))
    }
    oldmodel$model <- cbind(oldmodel$model, cutyhats = dfobs$cutyhats)
    oldmodel$model$grp <- as.factor(vapply(oldmodel$model$cutyhats, 
                                           function(x) which(observed[, 1] == x), 1))
    newmodel <- update(oldmodel, . ~ . + grp, data = oldmodel$model)
    if (class(oldmodel) == "polr") {
      LRstat <- oldmodel$deviance - newmodel$deviance
    } else if (class(oldmodel) == "clm") {
      LRstat <- abs(-2 * (newmodel$logLik - oldmodel$logLik))
    }
    PARAMETER <- g - 1
    PVAL <- 1 - pchisq(LRstat, PARAMETER)
    names(LRstat) <- "LR statistic"
    names(PARAMETER) <- "df"
    structure(list(statistic = LRstat, parameter = PARAMETER, 
                   p.value = PVAL, method = METHOD, data.name = DNAME, newmoddata = oldmodel$model, 
                   predictedprobs = yhat), class = "htest")
  }
