#' Bayesian variable selection for Gaussian linear models using PEP through 
#' exhaustive search or with the MC3 algorithm
#'
#' Given a formula and a data frame, performs Bayesian 
#' variable selection using either full enumeration and evaluation of 
#' all models in the model space (for model spaces of 
#' small--to--moderate dimension) or the MC3 algorithm (for 
#' model spaces of large dimension). Normal linear models are assumed for the data with the 
#' prior distribution on the model parameters (beta coefficients and 
#' error variance) being the PEP or the intrinsic. 
#' The prior distribution on the model space can be the uniform on models 
#' or the uniform on the model dimension (special case of the beta--binomial prior).
#' The model space consists of all possible models including an intercept term.
#'
#' @param formula A formula, defining the full model.
#' @param data A data frame (of numeric values), containing the data.
#' @param algorithmic.choice A character, the type of algorithm to be used
#' for selection: full enumeration and evaluation of all models or the MC3 algorithm. 
#' One of ``automatic'' (the choice is done automatically based on the number
#' of explanatory variables in the full model), ``full enumeration'' 
#' or ``MC3''. Default value=\code{"automatic"}.
#' @param intrinsic Logical, indicating whether the PEP 
#' (\code{FALSE}) or the intrinsic --- which   
#' is a special case of it --- (\code{TRUE}) should be used as prior on the  
#' regression parameters. Default value=\code{FALSE}.
#' @param reference.prior Logical, indicating whether the reference prior
#' (\code{TRUE}) or the dependence Jeffreys prior (\code{FALSE}) is used as 
#' baseline. Default value=\code{TRUE}.
#' @param beta.binom Logical, indicating whether the beta--binomial 
#' distribution (\code{TRUE}) or the uniform distribution (\code{FALSE})    
#' should be used as prior on the model space. Default value=\code{TRUE}.
#' @param ml_constant.term Logical, indicating whether the constant
#' (marginal likelihood of the null/intercept--only model) should be
#' included in computing the marginal likelihood of a model (\code{TRUE})  
#' or not (\code{FALSE}). Default value=\code{FALSE}.
#' @param burnin Non--negative integer, the burnin period for the MC3 algorithm.
#' Default value=1000.
#' @param itermc3 Positive integer (larger than \code{burnin}),
#' the (total) number of iterations for the MC3 algorithm. Default value=11000.
#'
#' @return \code{pep.lm} returns an object of class \code{pep}, 
#' i.e., a list with the following elements:
#' \item{models}{A matrix containing information about the models examined. 
#' In particular, in row \eqn{i} after representing model \eqn{i} with variable inclusion 
#' indicators, its marginal likelihood (in log scale), the R2, its dimension  
#' (including the intercept), the corresponding Bayes factor, 
#' posterior odds and its posterior probability are contained. The models
#' are sorted in decreasing order of the posterior probability. For the 
#' Bayes factor and the posterior odds, the comparison is made with the model 
#' with the highest posterior probability. The number of rows of this first list element
#' is \eqn{2^{p}} with full enumeration of all possible models, or equal to 
#' the number of unique models `visited' by the algorithm, if MC3 was run.
#' Further, for MC3, the posterior probability of a model corresponds to
#' the estimated posterior probability as this is computed by the relative
#' Monte Carlo frequency of the `visited' models by the MC3 algorithm.}
#' \item{inc.probs}{A named vector with the posterior inclusion probabilities of the 
#' explanatory variables.}
#' \item{x}{The input data matrix (of dimension \eqn{n\times p}), i.e., matrix containing 
#' the values of the \eqn{p} explanatory variables (without the intercept).}
#' \item{y}{The response vector (of length \eqn{n}).}
#' \item{fullmodel}{Formula, representing the full model.}
#' \item{mapp}{For \eqn{p\geq 2}, a matrix (of dimension \eqn{p\times 2}) containing the mapping between 
#' the explanatory variables and the Xi's, where the \eqn{i}--th explanatory variable
#' is denoted by Xi. If \eqn{p<2}, \code{NULL}.}
#' \item{intrinsic}{Whether the prior on the model parameters was PEP or intrinsic.}
#' \item{reference.prior}{Whether the baseline prior was the reference prior
#' or the dependence Jeffreys prior.}
#' \item{beta.binom}{Whether the prior on the model space was beta--binomial or
#'  uniform.}
#' When MC3 is run, there is the additional list element \code{allvisitedmodsM}, a matrix of 
#' dimension (\code{itermcmc}-\code{burnin}) \eqn{\times \,(p+2)} containing all `visited' models 
#' (as variable inclusion indicators together with their corresponding 
#' marginal likelihood and R2) by the MC3 algorithm after the burnin period.
#'
#' @details
#' The function works when \eqn{p\leq n-2}, where \eqn{p} is the number of explanatory variables
#' of the full model and \eqn{n} is the sample size.
#'
#' The reference model is the null model (i.e., intercept--only model).
#'
#' The case of missing data (i.e., presence of \code{NA}'s either in the  
#' response or the explanatory variables) is not currently supported. Further,
#' the data needs to be quantitative.
#'
#' All models considered (i.e., model space) include an intercept term.
#'
#' If \eqn{p>1}, the explanatory variables cannot have an exact linear relationship 
#' (perfect multicollinearity).
#'
#' The reference prior as baseline corresponds to hyperparameter values
#' \eqn{d0=0} and \eqn{d1=0}, while the dependence Jeffreys prior corresponds to 
#' model--dependent--based values for the hyperparameters \eqn{d0} and \eqn{d1},
#' see Fouskakis and Ntzoufras (2022) for more details.
#'
#' For computing the marginal likelihood of a model, Equation 16 of 
#' Fouskakis and Ntzoufras (2022) is used.
#'
#' When \code{ml_constant.term=FALSE} then the log marginal likelihood of a
#' model in the output is shifted by -logC1
#' (logC1: log marginal likelihood of the null model).
#'
#' When the prior on the model space is beta--binomial 
#' (i.e., \code{beta.binom=TRUE}), the following special case is used: uniform 
#' prior on model dimension. 
#'
#' If \code{algorithmic.choice} equals ``automatic'' then the choice of 
#' the selection algorithm is as follows: if \eqn{p < 20}, full enumeration 
#' and evaluation of all models in the model space is performed, 
#' otherwise the MC3 algorithm is used.
#' To avoid potential memory or time constraints, if \code{algorithmic.choice} 
#' equals ``full enumeration'' but \eqn{p \geq 20} then the MC3 algorithm is 
#' used instead (once issuing a warning message).
#' 
#' The MC3 algorithm was first introduced by Madigan and York (1995)
#' while its current implementation is described in the Appendix 
#' of Fouskakis and Ntzoufras (2022).
#'
#' @references Fouskakis, D. and Ntzoufras, I. (2022) Power--Expected--Posterior 
#' Priors as Mixtures of g--Priors in Normal Linear Models. 
#' Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
#' 
#' Madigan, D. and York, J. (1995) Bayesian Graphical Models for Discrete Data.
#' International Statistical Review, 63(2): 215–232. \doi{10.2307/1403615}
#'
#' @examples
#' data(UScrime_data)
#' res <- pep.lm(y~.,data=UScrime_data)
#' resu <- pep.lm(y~.,data=UScrime_data,beta.binom=FALSE)
#' resi <- pep.lm(y~.,data=UScrime_data,intrinsic=TRUE)
#' set.seed(123)
#' res2 <- pep.lm(y~.,data=UScrime_data,algorithmic.choice="MC3",itermc3=2000)
#' resj2 <- pep.lm(y~.,data=UScrime_data,reference.prior=FALSE,
#'                algorithmic.choice="MC3",burnin=20,itermc3=1800) 
#'
#' @export
pep.lm <- function(formula, data, 
                   algorithmic.choice="automatic", intrinsic=FALSE, 
                   reference.prior=TRUE, beta.binom=TRUE, 
                   ml_constant.term=FALSE, 
                   burnin=1000, itermc3=11000){
# check arguments' type
if(!inherits(formula,"formula"))    # check if formula is of type formula
  stop("formula should be of type formula.")
if(!is.data.frame(data)){
  message("Error: data should be a data frame.\n") 
  stop("Please respecify and call the function again.")
}
if (!(is.character(algorithmic.choice)&&
             algorithmic.choice%in%c("automatic","full enumeration","MC3"))){
  message("Error: algorithmic.choice should be one of automatic,
                  full enumeration or MC3.\n")
  stop("Please respecify and call the function again.")
}
mf <- model.frame(formula,data,na.action=NULL)
y <- model.response(mf, "numeric")
predctrs <- attr(terms(formula,data=data),"term.labels")
p <- length(predctrs)
if(p>=1){
 x <- mf[,predctrs]
 if(p==1)
   x <- matrix(x,ncol=1)
 if (algorithmic.choice=="automatic"){
  if (p<20)
    result <- full_enumeration_pepn(x, y, intrinsic, 
                                   reference.prior,
                                   beta.binom, ml_constant.term) else
    result <- mc3_pepn(x, y, 
                      intrinsic, reference.prior,
                      beta.binom, ml_constant.term,
                      burnin,itermc3)
 }else if (algorithmic.choice=="full enumeration"){
   if (p<20)
    result <- full_enumeration_pepn(x, y, intrinsic, 
                                   reference.prior,
                                   beta.binom, ml_constant.term)  else{
    message("The number of covariates does not allow full enumeration, thus 
the MC3 algorithm (with default parameters) will be applied instead.\n")
    result <- mc3_pepn(x, y, 
                      intrinsic, reference.prior,
                      beta.binom, ml_constant.term,
                      burnin,itermc3)
   }
 }
 else if (algorithmic.choice=="MC3")
    result <- mc3_pepn(x, y, 
                      intrinsic, reference.prior,
                      beta.binom, ml_constant.term,
                      burnin,itermc3)
}else{
  result <- full_enumeration_pepn(y=y,intrinsic=intrinsic,
                                 reference.prior=reference.prior,
                                 beta.binom=beta.binom,
                                 ml_constant.term=ml_constant.term)
}
if (length(result)==8)
 result <- pep(models=result$models, inc.probs=result$inc.probs, 
             x = result$x, y = result$y, fullmodel = formula(terms(formula,data=data)),
             mapp=result$mapp,
             intrinsic = result$intrinsic, reference.prior=result$reference.prior, 
             beta.binom = result$beta.binom) else
 result <- pep(models=result$models, inc.probs=result$inc.probs, 
             x = result$x, y = result$y, fullmodel = formula(terms(formula,data=data)),
             mapp=result$mapp,
             intrinsic = result$intrinsic, reference.prior = result$reference.prior, 
             beta.binom = result$beta.binom, allvisitedmodsM = result$allvisitedmodsM)
return(result)
}                                 # end function pep.lm