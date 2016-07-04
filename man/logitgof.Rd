\name{logitgof}
\alias{logitgof}
\title{Hosmer-Lemeshow Test for Binary and Multinomial Logistic Models}
\description{Performs the Hosmer and Lemeshow goodness-of-fit test for binary logistic regression models as well as its multinomial generalisation.}
\usage{logitgof(obs, exp, g = 10)}
\arguments{
  \item{obs}{ a vector of observed values. See details.}
  \item{exp}{ expected values fitted by the model.}
  \item{g}{ number of quantiles of risk, 10 by default.}
}
\details{
The Hosmer-Lemeshow test is a goodness-of-fit test for binary logistic regression models which was generalised by Fagerland, Hosmer and Bofin (2008) to accommodate multinomial logistic models. \code{logitgof} is capable of performing both tests, which are decribed in Hosmer, Lemeshow and Sturdivant (2013). Essentially, they compare observed with expected frequencies and compute a test statistic which is roughly distributed according to the chi-squared distribution. The degrees of freedom depend upon the number of quantiles used and the number of outcome categories. A non-significant p value indicates that there is no evidence that the observed and expected frequencies differ.

The vector of observed values passed to \code{logitgof} can be either a factor with a number of levels or a vector coded 1/0. If \code{obs} is a factor with three or more levels, the multinomial version of the test is run. If \code{obs} is a vector of 1s and 0s or a factor vector with 2 levels, then the binary version of the test is run.
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
}
\references{
\itemize{
\item Fagerland MW, Hosmer DW, Bofin AM. Multinomial goodness-of-fit tests for logistic regression models. Statistics in medicine 2008;27(21):4238-53.
\item Hosmer DW, Lemeshow S, Sturdivant RX. Applied Logistic Regression, 3rd Edition. 2013. New York, USA: John Wiley and Sons.
}
}
\author{Matthew Alexander Jay, with code adapted from the \code{hoslem.test} function written by Peter Solymos for the \code{ResourceSelection} package.}
\examples{
library(nnet)
library(foreign)
ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
ml$prog2 <- relevel(ml$prog, ref = "academic")

# Multinomial model
m <- multinom(prog2 ~ ses + write, data = ml)
summary(m)
logitgof(ml$prog2, fitted(m))
logitgof(ml$prog2, fitted(m))$exp
logitgof(ml$prog2, fitted(m), g = 5)$exp
logitgof(ml$prog2, fitted(m), g = 5)

# Binary logistic model
# Outcome with factor names
m2 <- glm(female ~ ses + write, data=ml, family="binomial")
logitgof(ml$female, fitted(m2))

# Outcome with 1/0 coding
ml$f2 <- 1
ml[ml$female=="male", ]$f2 <- 0
m3 <- glm(f2 ~ ses + write, data=ml, family="binomial")
logitgof(ml$f2, fitted(m3))
}
\keyword{regression}
