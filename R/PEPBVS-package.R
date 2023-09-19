#' @docType package
#' @name PEPBVS-package
#' @aliases PEPBVS-package
#' @title Bayesian variable selection using power-expected-posterior prior
#'
#' @description Performs Bayesian variable selection under normal linear
#' models for the data with the model parameters following as prior either
#' the PEP or the intrinsic (a special case of the former). 
#' The prior distribution on model space is the 
#' uniform on model space or the uniform on model dimension (a special case 
#' of the beta-binomial prior). Posterior model probabilities and marginal
#' likelihoods can be derived in closed-form expressions under this setup. 
#' The selection can be done either with full enumeration of all possible 
#' models (for small–to–moderate model spaces) or using the MC3 algorithm 
#' (for large model spaces). Complementary functions for making predictions, 
#' as well as plotting and printing the results are also available.
#'
#' @references Fouskakis, D. and Ntzoufras, I. (2022) Power-Expected-Posterior 
#' Priors as Mixtures of g-Priors in Normal Linear Models. 
#' Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
#'
#' Fouskakis, D. and Ntzoufras, I. (2020) Bayesian Model Averaging Using 
#' Power-Expected-Posterior Priors. 
#' Econometrics, 8(2): 17. \doi{10.3390/econometrics8020017}
#'
#' @importFrom Matrix rankMatrix
#' @importFrom Rcpp evalCpp
#' @importFrom stats var
#' @useDynLib PEPBVS , .registration=TRUE 
NULL







