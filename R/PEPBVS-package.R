#' _PACKAGE
#' @name PEPBVS-package
#' @aliases PEPBVS-package
#' @title Bayesian variable selection using power--expected--posterior prior
#'
#' @description Performs Bayesian variable selection under normal linear
#' models for the data with the model parameters following as prior distributions 
#' either the PEP or the intrinsic (a special case of the former). 
#' The prior distribution on model space is the 
#' uniform over all models or the uniform on model dimension (a special case 
#' of the beta--binomial prior). Posterior model probabilities and marginal
#' likelihoods can be derived in closed--form expressions under this setup. 
#' The selection is performed by either implementing a full enumeration 
#' and evaluation of all possible models (for model spaces of 
#' small--to--moderate dimension) or using the MC3 algorithm 
#' (for model spaces of large dimension). Complementary functions for hypothesis testing, 
#' estimation and predictions under Bayesian model averaging, 
#' as well as plotting and printing the results are also available. Selected
#' models can be compared to those arising from other well--known priors.
#'
#' @references Bayarri, M., Berger, J., Forte, A. and Garcia--Donato, G. (2012) 
#' Criteria for Bayesian Model Choice with Application to Variable Selection. 
#' The Annals of Statistics, 40(3): 1550–1577. \doi{10.1214/12-AOS1013}
#'
#' Fouskakis, D. and Ntzoufras, I. (2022) Power--Expected--Posterior 
#' Priors as Mixtures of g--Priors in Normal Linear Models. 
#' Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
#'
#' Fouskakis, D. and Ntzoufras, I. (2020) Bayesian Model Averaging Using 
#' Power--Expected--Posterior Priors. 
#' Econometrics, 8(2): 17. \doi{10.3390/econometrics8020017}
#'
#' Garcia--Donato, G. and Forte, A. (2018) Bayesian Testing, 
#' Variable Selection and Model Averaging in Linear Models using R with 
#' BayesVarSel. The R Journal, 10(1): 155–174. 
#' \doi{10.32614/RJ-2018-021}
#'
#' Kass, R. and Raftery, A. (1995) Bayes Factors. 
#' Journal of the American Statistical Association, 90(430): 773–795. 
#' \doi{10.1080/01621459.1995.10476572}
#'
#' Ley, E. and Steel, M. (2012) Mixtures of g--Priors for Bayesian Model 
#' Averaging with Economic Applications. 
#' Journal of Econometrics, 171(2): 251–266. 
#' \doi{10.1016/j.jeconom.2012.06.009}
#'
#' Liang, F., Paulo, R., Molina, G., Clyde, M. and Berger, J. (2008) 
#' Mixtures of g Priors for Bayesian Variable Selection. 
#' Journal of the American Statistical Association, 103(481): 410–423.
#' \doi{10.1198/016214507000001337}
#'
#' Raftery, A., Madigan, D. and Hoeting, J. (1997) Bayesian Model Averaging 
#' for Linear Regression Models. 
#' Journal of the American Statistical Association, 92(437): 179–191.
#' \doi{10.1080/01621459.1997.10473615}
#'
#' Zellner, A. (1976) Bayesian and Non--Bayesian Analysis of the Regression 
#' Model with Multivariate Student--t Error Terms. 
#' Journal of the American Statistical Association, 71(354): 400–405. 
#' \doi{10.1080/01621459.1976.10480357}
#' 
#' Zellner, A. and Siow, A. (1980) Posterior Odds Ratios for Selected 
#' Regression Hypotheses.
#' Trabajos de Estadistica Y de Investigacion Operativa, 31: 585-603. 
#' \doi{10.1007/BF02888369}
#'
#' @importFrom Matrix rankMatrix
#' @importFrom mvtnorm rmvt
#' @importFrom Rcpp evalCpp
#' @importFrom stats var
#' @importFrom stats lm.fit
#' @importFrom stats model.frame
#' @importFrom stats model.response
#' @importFrom stats terms
#' @useDynLib PEPBVS , .registration=TRUE 
NULL







