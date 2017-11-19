\name{pulkrob}
\alias{pulkrob.chisq}
\alias{pulkrob.deviance}
\title{Pulkstenis-Robinson goodness of fit tests for ordinal response models.}
\description{Performs the Pulkstenis-Robinson goodness of fit chi-squared and deviance tests for ordinal logistic regression models.}
\usage{pulkrob.chisq(model, catvars)
pulkrob.deviance(model, catvars)
}
\arguments{
  \item{model}{ an ordinal response model. Must be an object of class \code{polr} or \code{clm}.}
  \item{catvars}{ a character vector containing the names of the categorical covariates.}
}
\details{
The Pulkstenis-Robinson tests are goodness of fit tests for ordinal logistic regression models. They are capable of accommodating models with continuous as well as categorical predictors. The data are partitioned according to observed covariate patterns using the categorical covariates only. Any unobserved covariate patterns are discarded. Only categorical predictors are used to avoid partitioning among an unacceptably high number of covariate patterns. Each subject is assigned an ordinal response score by summing the predicted probabilities of each subject for each outcome level multiplied by equally spaced integer weights. The covariate patterns are then split into two at the median score within each.

Based on this partitioning, observed and expected frequencies are calculated and the test statistic computed. This statistic is distributed by the chi-squared distribution with \eqn{(2I - 1)(J - 1) - k - 1} degress of freedom, where \code{I} is the number of covariate patterns, \code{J} is the number of of response categories and k is the number of categorical variables in the model.

It is recommended (Fagerland and Hosmer, 2016) that the Pulkstenis-Robinson tests be run alongside the Hosmer-Lemeshow test (\code{logitgof}) and the Lipsitz test (\code{lipsitz.test}).
}
\value{
A list of class \code{htest} containing:
\item{statistic}{ the chi-squared or deviance statistic.}
\item{parameter}{ degrees of freedom used.}
\item{p.value}{ the p-value.}
\item{method}{ a character string indicating the name of the test.}
\item{data.name}{ a character string indicating the model formula used.}
\item{observed}{ a \code{cast} generated data frame of observed frequencies.}
\item{expected}{ a \code{cast} generated data frame of expected frequencies.}
\item{stddiffs}{ a table of the standardised differences. See Hosmer, Lemeshow and Sturdivant (2013), p 162.}
}
\references{
\itemize{
\item Fagerland MW, Hosmer DW. Tests for goodness of fit in ordinal logistic regression models. Journal of Statistical Computation and Simulation 2016. DOI: 10.1080/00949655.2016.1156682.
\item Hosmer DW, Lemeshow S, Sturdivant RX. Applied Logistic Regression, 3rd Edition. 2013. New York, USA: John Wiley and Sons.
\item Pulkstenis E, Robinson TJ. Goodness-of-fit tests for ordinal response regression models. Statistics in Medicine 2004;23:999-1014.
}
}
\author{Matthew Alexander Jay, with code from \code{epiR::epi.cp()} by Mark Stevenson.}
\seealso{
\code{\link{logitgof}}, \code{\link{lipsitz.test}}.
}
\examples{
data(mtcars)

# using polr
mod1 <- polr(as.factor(gear) ~ mpg + cyl + vs, data = mtcars)
pulkrob.chisq(mod1, c("vs"))
pulkrob.deviance(mod1, c("vs"))

# using clm - ensure outcome variable is a factor before fitting the model
library(ordinal)
mtcars$gear <- as.factor(mtcars$gear)
mod2 <- clm(gear ~ mpg + cyl + vs, data = mtcars)
pulkrob.chisq(mod2, c("vs"))
pulkrob.deviance(mod2, c("vs"))
}
\keyword{regression}
