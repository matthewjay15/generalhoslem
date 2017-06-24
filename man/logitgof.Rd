\name{logitgof}
\alias{logitgof}
\title{Hosmer-Lemeshow Tests for Logistic Regression Models}
\description{Performs the Hosmer and Lemeshow goodness of fit tests for binary, multinomial and ordinal logistic regression models.}
\usage{logitgof(obs, exp, g = 10, ord = FALSE)}
\arguments{
  \item{obs}{ a vector of observed values. See details.}
  \item{exp}{ expected values fitted by the model. See details.}
  \item{g}{ number of quantiles of risk, 10 by default.}
  \item{ord}{ logical indicating whether to run the ordinal version, FALSE by default.}
}
\details{
\bold{The Hosmer-Lemeshow tests}
The Hosmer-Lemeshow tests are goodness of fit tests for binary, multinomial and ordinal logistic regression models. \code{logitgof} is capable of performing all three. Essentially, they compare observed with expected frequencies and compute a test statistic which is roughly distributed according to the chi-squared distribution. The degrees of freedom depend upon the number of quantiles used and the number of outcome categories. A non-significant p value indicates that there is no evidence that the observed and expected frequencies differ.

\bold{Binary version}
If \code{obs} is a vector of 1s and 0s or a factor vector with 2 levels, then the binary version of the test is run. \code{exp} must be the fitted values obtained from the model, which can be accessed using the \code{fitted()} function.

\bold{Multinomial version}
If \code{obs} is a factor with three or more levels and \code{ord = FALSE}, the multinomial version of the test is run. If  using the \code{mlogit} package to run a model, ensure \code{outcome = FALSE} in the \code{fitted()} function. See examples.

\bold{Ordinal version}
If \code{obs} is a factor with three or more levels and \code{ord = TRUE}, the ordinal version of the test is run. See examples for how to extract fitted values from models constructed using \code{polr} in package \code{MASS} or \code{clm} in \code{ordinal}.

Note that Fagerland and Hosmer recommend running the Hosmer-Lemeshow test for ordinal models alongisde the Lipsitz test (\code{lipsitz.test}) and Pulkstenis-Robinson tests (\code{pulkrob.chisq} and \code{pulkrob.deviance}), as each detects different types of lack of fit.
}
\value{
A list of class \code{htest} containing:
\item{statistic}{ the value of the relevant test statistic.}
\item{parameter}{ the number of degrees of freedom used.}
\item{p.value}{ the p-value.}
\item{method}{ a character string indicating whether the binary or multinomial version of the test was performed.}
\item{data.name}{ a character string containing the names of the data passed to \code{obs} and \code{exp}.}
\item{observed}{ a table of observed frequencies with \code{g} rows. Either an \code{xtabs} generated table (used in the binary version) or a \code{cast} generated data frame (multinomial version).}
\item{expected}{ a table of expected frequencies with \code{g} rows. Either an \code{xtabs} generated table or a \code{cast} generated data frame.}
\item{stddiffs}{ a table of the standardised differences. See Hosmer, Lemeshow and Sturdivant (2013), p 162.}
}
\references{
\itemize{
\item Fagerland MW, Hosmer DW, Bofin AM. Multinomial goodness-of-fit tests for logistic regression models. Statistics in Medicine 2008;27(21):4238-53.
\item Fagerland MW, Hosmer DW. A goodness-of-fit test for the proportional odds regression model. Statistics in Medicine 2013;32:2235-2249.
\item Fagerland MW, Hosmer DW. Tests for goodness of fit in ordinal logistic regression models. Journal of Statistical Computation and Simulation 2016. DOI: 10.1080/00949655.2016.1156682.
\item Hosmer DW, Lemeshow S, Sturdivant RX. Applied Logistic Regression, 3rd Edition. 2013. New York, USA: John Wiley and Sons.
}
}
\author{Matthew Alexander Jay, with code adapted from the \code{hoslem.test} function written by Peter Solymos for the \code{ResourceSelection} package.}
\seealso{
\code{\link{lipsitz.test}}, \code{\link{pulkrob.chisq}}.
}
\examples{
### The mtcars dataset is a terrible example to use due to its small
### size - some of the models will return warnings as a result.
## Binary model
# 1/0 coding
data(mtcars)
mod1 <- glm(vs ~ cyl + mpg, data = mtcars, family = binomial)
logitgof(mtcars$vs, fitted(mod1))

# factor name coding
mtcars$engine <- factor(ifelse(mtcars$vs==0, "V", "S"), levels = c("V", "S"))
mod2 <- glm(engine ~ cyl + mpg, data = mtcars, family = binomial)
logitgof(mtcars$engine, fitted(mod2))

## Multinomial model
# with nnet
library(nnet)
mod3 <- multinom(gear ~ mpg + cyl, data = mtcars)
logitgof(mtcars$gear, fitted(mod3))

# with mlogit
library(mlogit)
data("Fishing", package = "mlogit")
Fish <- mlogit.data(Fishing, varying = c(2:9), shape = "wide", choice = "mode")
mod4 <- mlogit(mode ~ 0 | income, data = Fish)
logitgof(Fishing$mode, fitted(mod4, outcome = FALSE))

## Ordinal model
# polr in package MASS
mod5 <- polr(as.factor(gear) ~ mpg + cyl, data = mtcars)
logitgof(mtcars$gear, fitted(mod5), g = 5, ord = TRUE)

# clm in package ordinal
library(ordinal)
mtcars$gear <- as.factor(mtcars$gear)
mod6 <- clm(gear ~ mpg + cyl, data = mtcars)
predprob <- data.frame(mpg = mtcars$mpg, cyl = mtcars$cyl)
fv <- predict(mod6, newdata = predprob, type = "prob")$fit
logitgof(mtcars$gear, fv, g = 5, ord = TRUE)
}
\keyword{regression}
