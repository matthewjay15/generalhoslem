# generalhoslem
Version 1.3.0 (28 Oct 2017). Directly incorporates the epi.cp() function from the epiR package. This removes the dependency on epiR. Some typographical/grammar errors in documentation fixed.

Version 1.3.1 (19 Nov 2017). Bug fix - could not estimate lipsitz.test() where clm model had only one predictor.

Version 1.3.2 (2 Dec 2017). Bug fix - same problem with clm models with one predictor in pulkrob.chisq() and pulkrob.deviance().

Version 1.3.3 (20 Sep 2018). Bug fix - problem with pulkrob.chisq() and pulkrob.deviance() returning an error where outcome variable has "1" as category.

Version 1.3.4 (3 June 2019). Bug fix - problem with Pulkstenis-Robinson tests giving the wrong chisq2 and deviance values.

Known issues (1 June 2019):

It has been spotted that the binary and ordinal Hosmer-Lemeshow test and Lipsitz test give slightly different results to their implementations in Stata (which are not mine). I think this has something to do with the way in which the contingency tables are constructed but it is not at all clear what is different as the implementation is (should be) identical. I have altered my code to exactly match what the Stata versions do in terms of sorting the groups but this has not fixed the issue.
