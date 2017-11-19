\name{lipsitz.test}
\alias{lipsitz.test}
\title{Lipsitz goodness of fit test for ordinal logistic models.}
\description{Performs the Lipsitz goodness of fit test for ordinal logistic regression models.}
\usage{lipsitz.test(model, g = 10)}
\arguments{
  \item{model}{ an ordinal logistic regression model. Must be an object of class \code{polr} or \code{clm}.}
  \item{g}{ number of quantiles of risk, defaults to 10.}
}
\details{
The Lipsitz test is a goodness of fit test for ordinal response logistic regression models. It involves firstly binning the observed data into percentiles in the same manner as for the multinomial Hosmer-Lemeshow test. Strictly speaking, the Lipsitz test first requires an ordinal score be assigned to each subject. This score is the sum of the predicted probabilities for each outcome level within each subject multiplied by equally spaced integer weights. This, however, has been shown to be idenitical to the more straightforward binning method (Fagerland and Hosmer, 2013). The user can specify the number of groups by assigning an integer value to \code{g}, which is 10 by default.

Given this partitioning of the data, dummy variables, \code{I}, are derived such that, for each group, \code{I = 1} if the subject is in region \code{g} and \code{I = 0} if not. The model is then re-fit with these dummy variables. If the model has good fit, then the coefficients for all these dummy variables simultaneously = 0. Lipsitz et al (1996) suggest that likelihood ratio, Wald or score tests can be used; \code{lipsitz.test} just uses the likelihood ratio test with \code{g-1} degrees of freedom.

Note that the outcome variable MUST be converted to a factor before running the model. Using \code{as.factor()} within the model function will cause an error because of the way in which \code{lipsitz.test} uses the \code{update()} function to re-fit the model.

It is recommended (Fagerland and Hosmer, 2016) that the Lipsitz test be run alongside the ordinal Hosmer-Lemeshow test (\code{logitgof}) and the Pulkstenis-Robinson tests (\code{pulkrob.chisq} and \code{pulkrob.deviance}).
}
\value{
A list of class \code{htest} containing:
\item{statistic}{ the value of the likelihood ratio statistic.}
\item{parameter}{ degrees of freedom used.}
\item{p.value}{ the p-value.}
\item{method}{ a character string indicating the name of the test.}
\item{data.name}{ a character string indicating the model formula used.}
\item{newmoddata}{ a data frame containing the data used in fitting the updated model (the original data plus the dummy variables).}
\item{predictedprobs}{ a data frame of predicted probabilities from the original model.}
}
\references{
\itemize{
\item Fagerland MW, Hosmer DW. A goodness-of-fit test for the proportional odds regression model. Statistics in Medicine 2013;32:2235-2249.
\item Fagerland MW, Hosmer DW. Tests for goodness of fit in ordinal logistic regression models. Journal of Statistical Computation and Simulation 2016. DOI: 10.1080/00949655.2016.1156682.
\item Lipsitz SR, Fitzmaurice GM, Molenberghs G. Goodness-of-Fit Tests for Ordinal Response Regression Models. Journal of the Royal Statistical Society 1996;45(2):175-190.
}
}
\author{Matthew Alexander Jay}
\seealso{
\code{\link{logitgof}}, \code{\link{pulkrob.chisq}}.
}
\examples{
data(mtcars)
mtcars$gear <- as.factor(mtcars$gear) # Outcome must be converted to factor 
                                      # before fitting model otherwise 
                                      # lipsitz.test() will fail.
mod1 <- polr(gear ~ mpg + cyl, data = mtcars)
lipsitz.test(mod1)
  
library(ordinal)
mod2 <- clm(gear ~ mpg, data = mtcars)
lipsitz.test(mod2)
}
\keyword{regression}
