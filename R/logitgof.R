logitgof <-
function (obs, exp, g = 10) 
{
    DNAME <- paste(deparse(substitute(obs)), deparse(substitute(exp)), 
        sep = ", ")
    yhat <- exp
    if (is.null(ncol(yhat))) {
        mult <- FALSE
    }
    else {
        if (ncol(yhat) == 1) {
            mult <- FALSE
        }
        else mult <- TRUE
    }
    n <- ncol(yhat)
    if (mult) {
        METHOD <- "Hosmer and Lemeshow test (multinomial generalisation)"
        qq <- unique(quantile(1 - yhat[, 1], probs = seq(0, 1, 
            1/g)))
        cutyhats <- cut(1 - yhat[, 1], breaks = qq, include.lowest = TRUE)
        dfobs <- data.frame(obs, cutyhats)
        dfobsmelt <- melt(dfobs, id.vars = 2)
        observed <- cast(dfobsmelt, cutyhats ~ value, length)
        dfexp <- data.frame(yhat, cutyhats)
        dfexpmelt <- melt(dfexp, id.vars = ncol(dfexp))
        expected <- cast(dfexpmelt, cutyhats ~ variable, sum)
        chisq <- sum((observed[, 2:ncol(observed)] - expected[, 
            2:ncol(expected)])^2/expected[, 2:ncol(expected)])
        PARAMETER <- (g - 2) * (ncol(yhat) - 1)
    }
    else {
        METHOD <- "Hosmer and Lemeshow test (binary model)"
    	if (is.factor(obs)) {
      	    y <- as.numeric(obs) - 1
    	} else {
      	    y <- obs
    	}
        qq <- unique(quantile(yhat, probs = seq(0, 1, 1/g)))
        cutyhat <- cut(yhat, breaks = qq, include.lowest = TRUE)
        observed <- xtabs(cbind(y0 = 1 - y, y1 = y) ~ cutyhat)
        expected <- xtabs(cbind(yhat0 = 1 - yhat, yhat1 = yhat) ~ 
            cutyhat)
        chisq <- sum((observed - expected)^2/expected)
        PARAMETER <- g - 2
    }
    PVAL <- 1 - pchisq(chisq, PARAMETER)
    names(chisq) <- "X-squared"
    names(PARAMETER) <- "df"
    structure(list(statistic = chisq, parameter = PARAMETER, 
        p.value = PVAL, method = METHOD, data.name = DNAME, observed = observed, 
        expected = expected), class = "htest")
}
