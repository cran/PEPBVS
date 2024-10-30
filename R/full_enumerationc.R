# Bayesian variable selection through exhaustive search
#
# Given a response vector and an input data matrix, performs Bayesian variable 
# selection using full enumeration of the model space. Normal linear 
# models are assumed for the data with the prior distribution on the model 
# parameters (beta coefficients and error variance) being the PEP or the intrinsic. 
# The prior distribution on the model space can be the uniform on the model 
# space or the uniform on the model dimension (special case of the beta-binomial prior).
# The model space consists of all possible models including an intercept term.
#
# @param x A matrix of numeric (of size nxp), input data matrix. 
# This matrix contains the values of the p explanatory variables 
# without an intercept column of 1's.
# @param y A vector of numeric (of length n), response vector.
# @param intrinsic Boolean, indicating whether the PEP 
# (\code{FALSE}) or the intrinsic - which   
# is a special case of it - (\code{TRUE}) should be used as prior on the  
# regression parameters. Default value=\code{FALSE}.
# @param reference.prior Boolean, indicating whether the reference prior
# (\code{TRUE}) or the dependence Jeffreys prior (\code{FALSE}) is used as 
# baseline. Default value=\code{TRUE}.
# @param beta.binom Boolean, indicating whether the beta-binomial 
# distribution (\code{TRUE}) or the uniform distribution (\code{FALSE})    
# should be used as prior on the model space. Default value=\code{TRUE}.
# @param ml_constant.term Boolean, indicating whether the constant
# (marginal likelihood of the null/intercept-only model) should be
# included in computing the marginal likelihood of a model (\code{TRUE})  
# or not (\code{FALSE}). Default value=\code{FALSE}.
#
# @return \code{full_enumeration_pep} returns an object of class \code{pep}, 
# i.e. a list with the following elements:
# \item{models}{A matrix containing information about the models examined. 
# In particular, in row i after representing the model i with variable inclusion 
# indicators, its marginal likelihood (in log scale), the R2, its dimension  
# (including the intercept), the corresponding Bayes factor, 
# posterior odds and its posterior probability are contained. The models
# are sorted in decreasing order of the posterior probability. For the 
# Bayes factor and the posterior odds, the comparison is done to the model 
# with the largest posterior probability.}
# \item{inc.probs}{A named vector with the posterior inclusion probabilities of the 
# explanatory variables.}
# \item{x}{The input data matrix (of size nxp).}
# \item{y}{The response vector (of length n).}
# \item{mapp}{For p>=2, a matrix (of size px2) containing the mapping between 
# the explanatory variables and the Xi's, since the i-th explanatory variable
# is denoted by Xi. If p<2, \code{NULL}.}
# \item{intrinsic}{Whether the prior on the model parameters was PEP or intrinsic.}
# \item{reference.prior}{Whether the baseline prior was the reference prior
# or the dependence Jeffreys prior.}
# \item{beta.binom}{Whether the prior on the model space was beta-binomial or
#  uniform.}
# 
# @details
# The function works when p<=n-2 where p is the number of explanatory variables
# and n is the sample size.
#
# It is suggested to use this function (i.e. enumeration of the model space) 
# when p is up to 20.
#
# The reference model is the null model (i.e. intercept-only model).
#
# The case of missing data (i.e. presence of \code{NA}'s either in the  
# input data matrix or the response vector) is not currently supported.
#
# All models considered (i.e. model space) include an intercept term.
#
# If p>1, the input data matrix needs to be of full rank.
#
# The reference prior as baseline corresponds to hyperparameter values
# d0=0 and d1=0, while the dependence Jeffreys prior corresponds to 
# model-dependent-based values for the hyperparameters d0 and d1,
# see Fouskakis and Ntzoufras (2022) for more details.
#
# For computing the marginal likelihood of a model, Equation 16 of 
# Fouskakis and Ntzoufras (2022) is used.
#
# When \code{ml_constant.term=FALSE} then the log marginal likelihood of a
# model in the output is shifted by -logC1
# (logC1: log marginal likelihood of the null model).
#
# When the prior on the model space is beta-binomial 
# (i.e. \code{beta.binom=TRUE}), the following special case is used: uniform 
# prior on model dimension.
#
# @references Fouskakis, D. and Ntzoufras, I. (2022) Power-Expected-Posterior 
# Priors as Mixtures of g-Priors in Normal Linear Models. 
# Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
#
# @examples
# data(UScrime_data)
# y <- UScrime_data[,"y"]
# X <- UScrime_data[,-15]
# res <- full_enumeration_pep(X,y)
# resu <- full_enumeration_pep(X,y,beta.binom=FALSE)
# resi <- full_enumeration_pep(X,y,intrinsic=TRUE)
# resiu <- full_enumeration_pep(X,y,intrinsic=TRUE,beta.binom=FALSE)
# resj <- full_enumeration_pep(X,y,reference.prior=FALSE)
#
# @seealso \code{\link{pep.lm}}, \code{\link{mc3_pep}}
# 
full_enumeration_pepn <- function(x,y, intrinsic=FALSE, reference.prior=TRUE,
                                 beta.binom=TRUE, ml_constant.term=FALSE){
# The function takes a response vector y and an input matrix x. Normal linear
# models are assumed for the data. Performs bayesian variable  
# comparison/selection through exhaustive search/full enumeration of all 
# possible models. The prior distribution of the beta coefficients can be one 
# of PEP (intrinsic=FALSE, default) or intrinsic (intrinsic=TRUE) while the prior 
# distribution on the different models can be one of uniform (beta.binom=FALSE) 
# or beta binomial (beta.binom=TRUE, default). The additional parameters d0, d1 
# (default to 0) are involved in the prior distributions, while k0 (default to 1) 
# is the dimension of the null model. Returns a list of ...
#
if(!missing(x)){                       # in the presence of at least one input variable
 if(is.data.frame(x))                  # if x is a data frame
  x <- as.matrix(x)                    # turn it to a matrix
 if(is.vector(x))                      # if x is a vector
  x <- matrix(x,ncol=1)                # turn it to a 1-column matrix
 namesinitial <- colnames(x)
 if (is.matrix(y)&&(ncol(y)==1))       # if the response vector is given as an
                                       # 1-column matrix
   y <- as.vector(y)                   # turn it to a vector
 # check that the function arguments are of the correct type
 checkinputvartype(x,y,intrinsic,reference.prior,
                   beta.binom,ml_constant.term)
 # check dimension compatibility between input matrix and response vector
 checkdimenscomp(x,y)
 # check that n>=p+2 for feasibility
 compareparamvssamplesize(x)
 checkmissingness(x,y)                 # check that there are no NAs
 p <- ncol(x)                          # number of explanatory variables
 if (p>1)                              # if more than one
  checkmatrixrank(x)                   # check if input matrix is of full rank
 n <- nrow(x)                          # sample size
 namesg <- colnames(x) <- paste("X",1:p,sep="")  # add column names (X1,X2,...)
                                       # names of explanatory variables
 
 # models in full enumeration
 nmodels <- 2^p                        # number of possible models
 gamma <- integer.base.b(0:(nmodels-1))# matrix with all possible models 
                                       # gamma[i,j]=1 if variable Xj is included 
                                       # in model i and 0 differently) 
 # Check if the centering is necessary now that the prediction was removed
 # and if not also remove the command below... 
 xi <- x                               # store it before scaling so that you
                                       # can save it afterwards
 x <- scale(x, center=TRUE, scale = FALSE)
                                       # compute for each model its log marginal 
                                       # likelihood and R^2 (code in C++)
 res1 <- full_enumeration_pepc(x,gamma,y,intrinsic=intrinsic,reference.prior)
 marginalv <- res1[[1]]                # log marginal likelihoods for the models
 Rsqv <- res1[[2]]                     # & corresponding R^2's
 if (ml_constant.term){                # add constant term logC1 (to each 
                                       # marginal likelihood) if set to TRUE 
  if (reference.prior)
   d0 <- 0 else
   d0 <- 1
  marginalv <- marginalv + constant.marginal.likelihood(y, d0)
 }
                                       # matrix with model info & marginal likelihoods
                                       # R^2's and models' dimension (including
                                       # the intercept)
 result <- cbind(gamma,marginalv,Rsqv,rowSums(gamma)+1)
 # --------------------------------------------- 
 # compute log posterior prob 
 # --------------------------------------------- 
 logpost <- marginalv                  # log posterior prob will be equal to the 
                                       # marginal likelihood for uniform model prior
 if (beta.binom){                      
  betabinomv <- -lchoose(p,0:p)
  logpriorprobs <- betabinomv[rowSums(gamma)+1] # log prior for the models
  logpost <- logpost+logpriorprobs     # log posterior
 }
 # --------------------------------------------- 
 # sort results according to posterior prob 
 # --------------------------------------------- 
                                       # sort log posterior in decreasing order
 index <- order(logpost, decreasing=TRUE) # of value
 result <- result[index,]              # and the whole result/matrix accordingly
 logpost <- logpost[index] 
 # --------------------------------------------- 
 # Calculation of BFs and Posterior probs  
 # --------------------------------------------- 
 logBF <- result[,p+1] -result[1,p+1]  # log BF w.r.t the MAP model
 BF <- exp(logBF)                      # BF
 if (beta.binom)                       # for beta-binomial model prior
   PO <- exp(logpost-logpost[1]) else  # POs  
   PO <- BF                            # PO = BF for uniform model prior
 sumPOs <- sum(PO)                     # POs sum
 post.prob <- PO/sumPOs                # posterior probs
 # compute negative & inverse already above??
                                       # add to the result the following info:
                                       # 1/BF, 1/PO and posterior prob
 models <- cbind(result, 1/BF, 1/PO, post.prob)
                                       # add column names
 colnames(models)<- c(namesg,"PEP log-marginal-like", "R2", "dim", 
                      "BF", "PO", "Post.Prob") 
 # --------------------------------------------- 
 # Calculation of Inclusion Probabilities  
 # --------------------------------------------- 
 g <- result[,1:p]
 if (p==1)
  inc.probs <- sum(g*post.prob)else
 inc.probs <- colSums(g*post.prob)
 names(inc.probs) <- namesg
 nmsmapp <- NULL
 if(p>=2){
  nmsmapp <- cbind(namesinitial,namesg)
  colnames(nmsmapp) <- c("Initial names", "Names")
 }
}else{                                 # in case the input matrix is completely 
                                       # missing i.e. (p=0)--treat it differently
 if (is.matrix(y)&&(ncol(y)==1))
   y <- as.vector(y)
                                       # check the remaining arguments
 checkinputvartype(y=y,intrinsic=intrinsic,reference.prior=reference.prior,
                   beta.binom=beta.binom,constant.term=ml_constant.term)
 if (sum(is.na(y))>0){
  message("Error: the output vector should not contain missing values.\n")
  stop("Please respecify and call the function again.")
 }
 if (length(y)<2){    # if the sample size does not exceed by 2 at least
                      # the number of explanatory variables
  message("Error: the number of parameters (2) exceeds the sample size.
       This setup is not supported by the package.\n")
  stop("Please respecify and call the function again.") 
 }
                                       # save the result (i.e. marglik = R^2=0 
                                       # and the rest, that is 
                                       # dim,BF,PO,postprob = 1)
 models <- cbind(0,0,1,1,1,1)
 if (ml_constant.term){                # add constant term if set to TRUE
   if (reference.prior)
    d0 <- 0 else
    d0 <- 1
   models[1,1] <- models[1,1] + constant.marginal.likelihood(y, d0)
 }
 colnames(models)<- c("PEP log-marginal-like", "R2", "dim", 
                      "BF", "PO", "Post.Prob")
 inc.probs <- NULL                     # inclusion probabilities not involving  
                                       # the intercept do not exist
 xi <- NULL
 nmsmapp <- NULL
}                                      # generate result - object of type pep
result <- list(models=models, inc.probs=inc.probs, x = xi, y = y, 
             mapp=nmsmapp,
             intrinsic = intrinsic, reference.prior=reference.prior, beta.binom = beta.binom)
return(result)
}                                      # end function full_enumeration_pep




#' @param ... Deprecated.
#'
#' @rdname PEPBVS-deprecated
#' @export

full_enumeration_pep <- function(...){
 .Deprecated("pep.lm",
             msg=paste("The function full_enumeration_pep is deprecated,", 
                       "please use the function pep.lm instead."))
}
