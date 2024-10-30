# PEPBVS 2.0

* Added function `peptest` for hypothesis testing

* Added function `estimation.pep` which can be used for estimation under Bayesian model averaging

* Added function `posteriorpredictive.pep` which simulates values from the posterior predictive distribution under Bayesian model averaging 

* Added function `comparepriors.lm` which returns selected models under different choices of prior on the model parameters and the model space

* The argument `xnew` of S3 method `predict` has to be a data frame instead of a matrix.

* The older functions `full_enumeration_pep` and `mc3_pep` were unified and replaced by the new function `pep.lm`.
The main differences are the following: a) `pep.lm` takes as arguments a formula and a data frame instead of a
matrix (input data matrix) and a vector (response vector) and b) with the new argument
`algorithmic.choice="automatic"`, the choice of the model selection algorithm (full enumeration or MC3)
can also be done automatically and does not necessarily have to be chosen by the user.  
Finally, the new function returns few additional list elements: `fullmodel` representing the full model that was used,
`mapp` containing the mapping between Xi's and the original explanatory variable names
and `allvisitedmodsM` containing all `visited' models after the burnin period (relevant for MC3). The older functions `full_enumeration_pep` and `mc3_pep` return a warning message
that are deprecated.