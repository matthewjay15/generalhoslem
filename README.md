# generalhoslem
Version 1.3.0 (28 Oct 2017). Directly incorporates the epi.cp() function from the epiR package. This removes the dependency on epiR. Some typographical/grammar errors in documentation fixed.

Version 1.3.1 (19 Nov 2017). Bug fix - could not estimate lipsitz.test() where clm model had only one predictor.

Version 1.3.2 (2 Dec 2017). Bug fix - same problem with clm models with one predictor in pulkrob.chisq() and pulkrob.deviance().

Version 1.3.3 (20 Sep 2018). Bug fix - problem with pulkrob.chisq() and pulkrob.deviance() returning an error where outcome variable has "1" as category.


Known issues to be fixed (1 June 2019):

There is an issue where the Pulkstenis-Robinson tests give the wrong results. The source of the error is in the epi.cp() code where it sorts the patterns. Removing the sorting fixes the issue.

It has also been spotted that the ordinal Hosmer-Lemeshow test and Lipsitz test give slightly different results to the implementation in Stata (which is not mine). This has something to do with the way in which the contingency tables are constructed but it is not at all clear what is different as the implementation is (should be) identical. The difference is very minor.
