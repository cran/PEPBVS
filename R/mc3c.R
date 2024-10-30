# Bayesian variable selection with MC3 algorithm
#
# Given a response vector and an input data matrix, performs Bayesian variable 
# selection using the MC3 algorithm. Normal linear models are assumed for 
# the data with the prior distribution on the model parameters 
# (beta coefficients and error variance) being the PEP or the intrinsic. The prior 
# distribution on the model space can be the uniform on the model 
# space or the uniform on the model dimension (special case of the beta-binomial prior).
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
# @param burnin Non-negative integer, the burnin period for the MC3 algorithm.
# Default value=1000.
# @param itermc3 Positive integer (larger than \code{burnin}),
# the (total) number of iterations for the MC3 algorithm. Default value=11000.
#
# @return \code{mc3_pep} returns an object of class \code{pep}, 
# as this is described in detail in \code{\link{full_enumeration_pep}}. The
# difference is that here the number of rows of the first list element
# is not \eqn{2^{p}} but the number of unique models `visited' by the 
# MC3 algorithm. Further, the posterior probability of a model corresponds to
# the estimated posterior probability as this is computed by the relative
# Monte Carlo frequency of the `visited' models by the MC3 algorithm. Finally,
# there is the additional list element \code{allvisitedmodsM}, a matrix of 
# size (\code{itermcmc}-\code{burnin}) x (p+2) containing all `visited' models 
# (as variable inclusion indicators together with their corresponding 
# marginal likelihood and R2) by the MC3 algorithm after the burnin period.  
#
# @details
# The function works when p<=n-2 where p is the number of explanatory variables
# and n is the sample size.
#
# It is suggested to use this function (i.e. MC3 algorithm) when p is 
# larger than 20.
#
# The reference model is the null model (i.e. intercept-only model).
#
# The case of missing data (i.e. presence of \code{NA}'s either in the  
# input matrix or the response vector) is not currently supported.
#
# The intercept term is included in all models.
#
# If p>1, the input matrix needs to be of full rank.
#
# The reference prior as baseline corresponds to hyperparameter values
# d0=0 and d1=0, while the dependence Jeffreys prior corresponds to 
# model-dependent-based values for the hyperparameters d0 and d1,
# see Fouskakis and Ntzoufras (2022) for more details.
#
# The MC3 algorithm was first introduced by Madigan and York (1995)
# while its current implementation is described in the Appendix 
# of Fouskakis and Ntzoufras (2022).
#
# When \code{ml_constant.term=FALSE} then the log marginal likelihood of a
# model in the output is shifted by -logC1
# (logC1: marginal likelihood of the null model).
#
# When the prior on the model space is beta-binomial 
# (i.e. \code{beta.binom=TRUE}), the following special case is used: uniform 
# prior on model size.
#
# @references Fouskakis, D. and Ntzoufras, I. (2022) Power-Expected-Posterior 
# Priors as Mixtures of g-Priors in Normal Linear Models. 
# Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
# 
# Madigan, D. and York, J. (1995) Bayesian Graphical Models for Discrete Data.
# International Statistical Review, 63(2): 215–232. \doi{10.2307/1403615}
#
# @examples
# data(UScrime_data)
# y <- UScrime_data[,"y"]
# X <- UScrime_data[,-15]
# set.seed(123)
# res <- mc3_pep(X,y,itermc3=3000)
# resu <- mc3_pep(X,y,beta.binom=FALSE,itermc3=3000)
# resj <- mc3_pep(X,y,reference.prior=FALSE,burnin=500,itermc3=2200)
# 
# @seealso \code{\link{pep.lm}}, \code{\link{full_enumeration_pep}}
# 
mc3_pepn <- function(x,y, intrinsic=FALSE, reference.prior=TRUE, beta.binom=TRUE, 
                    ml_constant.term=FALSE,burnin=1000,itermc3=11000){
# The function takes a response vector y and an input matrix x. Normal linear
# models are assumed for the data. Performs bayesian variable  
# comparison/selection through application of the MC3 algorithm. The number of
# is itermc3 (default to 11000).
# The prior distribution of the beta coefficients can be one 
# of PEP (intrinsic=FALSE, default) or intrinsic (intrinsic=TRUE) while the prior 
# distribution on the different models can be one of uniform (beta.binom=FALSE) 
# or beta binomial (beta.binom=TRUE, default). The additional parameters d0, d1 
# (default to 0) are involved in the prior distributions.
# Returns ...
#
if(!missing(x)){                        # in the presence of at least one input variable
 if(is.data.frame(x))                   # if x is a data frame
  x <- as.matrix(x)                     # turn it to a matrix
 if(is.vector(x))                       # if x is a vector
  x <- matrix(x,ncol=1)                 # turn it to a 1-column matrix
 namesinitial <- colnames(x)
 if (is.matrix(y)&&(ncol(y)==1))        # if the response is given as an
                                        # 1-column matrix
   y <- as.vector(y)                    # turn it to a vector
 # check that the function arguments are of the correct type
 checkinputvartype(x,y,intrinsic,reference.prior,
                   beta.binom,ml_constant.term,
                   burnin,itermc3)
 # check dimension compatibility between input matrix and response vector
 checkdimenscomp(x,y)
 # check that n>=p+2 for feasibility
 compareparamvssamplesize(x)
 checkmissingness(x,y)                  # check that there are no NAs
 p <- ncol(x)                           # number of explanatory variables
 if (p>1)                               # if more than one
  checkmatrixrank(x)                    # check if input matrix is of full rank
 namesg <- colnames(x) <- paste("X",1:p,sep="")   # add column names (X1,X2,...)
                                        # names of explanatory variables
 # Null model
 gamma.null <- rep(0,p)                 # vector of indicator variables for null model         
                                        # marginal likelihood for null model
                                        # (or posterior in case of uniform prior on models)
 logmarg.null <- 0 
 Rsq.null <- 0                          # R^2 for null: 0   
 
 # start from null model
 current_gamma <- gamma.null            # current model
                                        # current log marginal likelihood
 current_logpost <- current_marginal <- logmarg.null        
 current_Rsq <- Rsq.null
                                        # store the result
 resultnull <- c(current_gamma,current_marginal,current_Rsq)
 xi <- x                                # store it before scaling so that you
                                        # can save it afterwards
 x <- scale(x, center=TRUE, scale = FALSE)
                                        # apply MC3 algorithm and return all
                                        # 'visited' models together with their log
                                        # marginal likelihood and R^2
 resultc <- mc3_pepc(x,y,beta.binom,itermc3-1,current_Rsq, 
                    current_marginal,current_logpost,current_gamma,
                    intrinsic,reference.prior)
                                        # combine the results of the 'visited' models
 result <- rbind(resultnull,resultc)    # with the results of the starting/null model
 if (burnin>0)                          # if there is burnin period
   result <- result[-(1:burnin),]       # ignore the first burnin iterations
 diff <- itermc3-burnin                 # nbr of remaining iterations 
 if(diff==1)                            # if only 1 remaining iteration
   result <- matrix(result,nrow=1)      # turn the vector to a 1-row matrix
 if (ml_constant.term){                 # add logC1 if constant term set to TRUE
   if (reference.prior)
    d0 <- 0 else
    d0 <- 1
   result[,p+1] <- result[,p+1] + constant.marginal.likelihood(y, d0)
 }
 resultallvisitedmods <- result
 # --------------------------------------------- 
 # Estimation of Inclusion Probabilities  
 # --------------------------------------------- 
 if (p==1){                             # if only 1 input variable
   inc.probs  <- mean(result[,1])       # mean of the 1st column 
 }else if (diff==1)                     # if only 1 remaining iteration 
   inc.probs  <- result[,1:p] else      # take the 1st row   
 inc.probs  <- apply(result[,1:p],2,mean)# mean (proportion of times MC3 'visited'
                                        # that input variable) of the p columns
 names(inc.probs) <- namesg             # add names
 # write each model represented with inclusion indicator variables
 # as a vector of character/string (e.g. "0100...1")
 if (p==1){                             # if only 1 input variable
   rownames(result) <- apply(matrix(result[,1],ncol=1),1,
                             paste,collapse='')
 }else if (diff==1)                     # if only 1 remaining iteration
   rownames(result) <- paste(result[,1:p],collapse='') else
 rownames(result) <- apply(result[,1:p],1,paste,collapse='') # general case
                                        # compute frequencies that MC3 visited
                                        # the different models
 ftbs <- sort(table(rownames(result))/(itermc3-burnin),decreasing=T)
 if (length(ftbs)==1)                   # if only 1 visited model
                                        # turn the vector to a 1-row matrix
   result <- matrix(result[names(ftbs),],nrow=1) else
 result <- result[names(ftbs),]         # reduce matrix to the unique visited models 
 # models' dimension (including intercept)
 if(p==1){                              # if only 1 input variable
  mdim <- matrix(result[,1],ncol=1)+1
 }else if (dim(result)[1]==1)           # if only 1 model
  mdim <- sum(result[,1:p])+1 else 
 mdim <- rowSums(result[,1:p])+1        # general case
 # --------------------------------------------- 
 # Calculation of BFs and Posterior odds  
 # --------------------------------------------- 
 logBF <- result[,p+1] -result[1,p+1]   # log BF w.r.t the MAP model
 BF <- exp(logBF)                       # BF
 PO <- ftbs/ftbs[1]                     # estimated POs (ratio based on estimated
                                        # posterior probs)
                                        # add to the result the following info:
                                        # model dimension, 1/BF, 1/PO and posterior prob
 result <- cbind(result,mdim,1/BF, 1/PO, ftbs)
                                        # add column names
 colnames(result) <- c(namesg,"PEP log-marginal-like", "R2","dim",
                       "BF", "PO", "Post.Prob")
 rownames(result) <- NULL 
 colnames(resultallvisitedmods) <- c(namesg,"PEP log-marginal-like", "R2")
 nmsmapp <- NULL
 if(p>=2){
  nmsmapp <- cbind(namesinitial,namesg)
  colnames(nmsmapp) <- c("Initial names", "Names")
 }
}else{                                  # in case the input matrix is completely 
                                        # missing i.e. (p=0)--treat it differently

 if (is.matrix(y)&&(ncol(y)==1))
   y <- as.vector(y)
                                        # check the remaining arguments
 checkinputvartype(y=y,intrinsic=intrinsic,reference.prior,
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
 result <- cbind(0,0,1,1,1,1)
 if (ml_constant.term){                 # add constant term if set to TRUE 
   if (reference.prior)
    d0 <- 0 else
    d0 <- 1
   result[1,1] <- result[1,1] + constant.marginal.likelihood(y, d0)
 }
 colnames(result)<- c("PEP log-marginal-like", "R2", "dim", 
                      "BF", "PO", "Post.Prob")
 inc.probs <- NULL                      # inclusion probabilities not involving  
                                        # the intercept do not exist
 xi <- NULL
 nmsmapp <- NULL
 resultallvisitedmods <- NULL
}                                       # generate result - object of type pep
resultlst <- list(models=result, inc.probs=inc.probs, x = xi, y = y, 
             mapp=nmsmapp,
             intrinsic = intrinsic, reference.prior=reference.prior, 
             beta.binom = beta.binom, allvisitedmodsM = resultallvisitedmods)
return(resultlst)
}                                       # end function mc3_pep


#' @param ... Deprecated.
#'
#' @rdname PEPBVS-deprecated
#' @export

mc3_pep <- function(...){
 .Deprecated("pep.lm",
             msg=paste("The function mc3_pep is deprecated,", 
                       "please use the function pep.lm instead."))
}
