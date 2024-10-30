#' Plots for object of class pep
#'
#' Generates four plots related to an object of class pep. In particular,
#' the first one is a plot of the residuals against fitted values under 
#' Bayesian model averaging. The second plots the cumulative posterior 
#' probability of the top models (those with cumulative posterior probability 
#' larger than 0.99). The third plot depicts the marginal likelihood 
#' (in log scale) of a model against its dimension while 
#' the fourth plot shows the posterior inclusion probabilities
#' of the explanatory variables (with those exceeding 0.5 marked in red).
#'
#' @param x An object of class pep (e.g., output of \code{pep.lm}).
#' @param ... Additional graphical parameters to be passed to plotting functions.
#'
#' @return No return value, used for figure generation.
#'
#' @details
#' Let \eqn{k} be the number of models with cumulative posterior probability up 
#' to 0.99. Then, the second plot depicts the cumulative posterior probability 
#' of the top \eqn{(k+1)} models.
#'
#' In the special case of no explanatory variables, the fourth plot with the
#' posterior inclusion probabilities is not generated.
#'
#' @examples
#' data(UScrime_data)
#' res <- pep.lm(y~.,data=UScrime_data)
#' plot(res)
#'
#' @seealso \code{\link{image.pep}}
#' @importFrom graphics axis par
#' @method plot pep
#' @export
plot.pep <- function(x,...){
# Takes the result of full enumeration or mc3 algorithm res and generates
# 4 plots:
# a) residuals against fitted values under BMA model
# b) cumulative prob summing up to 0.99 across the ordered models
# c) log marginal likelihood vs model dimension
# d) inclusion probs for the different explanatory variables 
#
 if(!inherits(x,"pep"))               # check if x is an object of type pep
   stop("The input should be an object of type pep.")
 # 1st plot: residuals against fitted values
 y <- x$y
                                      # fitted values under (full) BMA
 BMA.ypred <- predict.pep(x,cumul.prob=1)

 oldpar <- par(mfrow=c(2,2),cex.main=1) 
 on.exit(par(oldpar)) 
 resid <- BMA.ypred-y                 # compute residuals
                                      # generate the plot
 plot(BMA.ypred,resid,xlab="Fitted values under BMA",
      ylab="Residuals",main="Residuals vs Fitted",...)
 # 2nd plot: cumulative prob up to 0.99 for the different models
 modelsmr <- x$models                 # matrix with models' info
 post.prob <- modelsmr[,"Post.Prob"]  # vector with posterior probs
 cumprob <- cumsum(post.prob)         # cumulative prob
 ind <- which(cumprob>0.99)[1]        # first model for which cumulative prob > 0.99
 cumprobf <- cumprob[1:ind]
 if (ind>1)                           # if more than one models sum above 0.99
   plot(cumprobf,                     # show cumulative prob for them
     xlab="Model order",ylab="Cumulative probability",
     main="Model Probabilities",cex=0.8,...) else
   plot(1,cumprob[1],                 # else show prob for the top model
     xlab="Model order",ylab="Cumulative probability",
     main="Model Probabilities",xlim=c(0,1),ylim=c(0,1),...)
 # 3rd plot: log marginal likelihood vs model dimension
 plot(modelsmr[,"dim"],modelsmr[,"PEP log-marginal-like"],
     xlab="Model dimension",ylab="Log(Marginal)",
     main="Model Complexity",cex=0.8,pch=20,...)
 # 4th plot: inclusion probs
 inc.probs <- x$inc.probs             # inclusion probs
 p <- length(inc.probs)               # number of explanatory variables
 if (p>0){                            # if there are input variables
  vcol <- rep("black",p)              # color vector
  vcol[inc.probs>0.5] <- "red"        # red to be used for those variables with
                                      # inclus prob > 0.5
                                      # generate the plot
  plot(inc.probs,xaxt='n',type="h",col=vcol,
     xlab="Variable",ylab="Marginal Inclusion Probability",
     main="Inclusion probabilities",...)  
                                      # add x-axis manually
  axis(1,at=1:p,labels=names(inc.probs),
      lwd.ticks=0,                    # no ticks
      cex.axis=0.6,las=2)             # x-axis names cex, perpendicular to axis
 }
}                                     # end function plot.pep


#' Printing object of class pep
#'
#' For each of the top models (shown in columns), the following information is
#' printed: the model representation using variable inclusion indicators, 
#' its marginal likelihood (in log scale), the R2, the model dimension, the Bayes 
#' factor, posterior odds (comparison made with the highest posterior 
#' probability model) and posterior probability. An additional 
#' column with the posterior inclusion probabilities of the explanatory variables is 
#' also printed. 
#'
#' @param x An object of class pep (e.g., output of \code{pep.lm}).
#' @param n.models Positive integer, the number of top models for which information 
#' is provided. Default value=5.
#' @param actual.PO Logical, relevant for the MC3 algorithm. If \code{TRUE}
#' then apart from the estimated posterior odds, the actual posterior
#' odds of the MAP model versus the top models (i.e., ratios based on the marginal likelihood 
#' times prior probability) are also printed --- which could be used as a 
#' convergence indicator of the algorithm. Default value=\code{FALSE}.
#' @param digits Positive integer, the number of digits for printing numbers. 
#' Default value=\code{max(3L, getOption("digits") - 3L)}.
#' @param ... Additional parameters to be passed to 
#' \code{print.default}.
#'
#' @return No return value, used for printing the results on the R console.
#'
#' @details
#' The number of models for which information is provided, is computed as the minimum
#' between the number asked by the user and the number of models present in
#' the object \code{x}.
#'
#' @examples
#' data(UScrime_data)
#' res <- pep.lm(y~.,data=UScrime_data)
#' print(res)
#'
#' @method print pep
#' @export
print.pep <- function(x, n.models=5, actual.PO=FALSE, 
                      digits = max(3L, getOption("digits") - 3L),...){
# Takes the result of full enumeration or mc3 algorithm res
# and the number of top models n.models (5 by default) we are interested in 
# and prints relevant info about them. For each model, its log(marginal), R2,
# dimension, post.prob, BF, PO as well as which explanatory variables 
# are included is given. The inclusion probs are also printed as a separate 
# column. The number of significant digits used for printing is determined by
# the argument digits (default value max(3L, getOption("digits") - 3L)). In case 
# the actual POs are further asked by the user (actual.PO=T)-(this would be relevant 
# in the case of MC3 algorithm where estimated POs are available)- then these 
# are computed (ratio of products of actual marginal by actual prior) and printed as well. 
#
 if(!inherits(x,"pep"))               # check if x is an object of type pep
   stop("The object to be printed should be an object of type pep.")
                                      # check the remaining arguments
 checkinputvartypeprint(n.models, actual.PO, digits)
 modelsmr <- x$models                 # matrix with models' info
 inc.probs <- x$inc.probs             # inclusion probs
 ncolres <- dim(modelsmr)[2]          # number of columns
 p <- length(inc.probs)               # number of explanatory variables
 diff <- ncolres-p                    # difference
                                      # add that many NAs (to fill in the difference)
 inc.probsa <- c(inc.probs, rep(NA, diff))
                                      # number of top models (note: take the min
                                      # because the object res might contain less
                                      # rows/models than the top asked)
 n.models <- min(n.models,dim(modelsmr)[1])
                                      # reduced matrix with models' info
                                      # focusing on the top n.models models
 resred <- modelsmr[1:n.models,]
 if(n.models==1){                     # special case where resred is a vector
   resred <- matrix(resred,nrow=1)    # turn it to a one-row matrix
 }
 resa <- cbind(inc.probsa,t(resred))  # add inclusion probs info as first column
 rownames(resa) <- colnames(modelsmr) # add row names
                                      # and column names
 colnames(resa) <- c("Incl Probs", paste("model", 1:n.models))
 # The following is relevant in the case of MC3 algorithm
 # where the actual POs are further asked by the user
 if(actual.PO==T){                    # (and not only the estimated ones)
                                      # log BF for top models
  logBF <- resred[,"PEP log-marginal-like"]-resred[1,"PEP log-marginal-like"]
  beta.binom <- x$beta.binom
  if(beta.binom==T){                  # if model prior is beta-binomial
   betabinomv <- -lchoose(p,0:p)      # log prior probs (shifted)
   dims <- resred[,"dim"]             # models' dimension
                                      # log prior odds
   logpriorodds <- betabinomv[dims]-betabinomv[dims[1]]
   logpostodds <- logBF+logpriorodds  # log POs
   PO2 <- exp(logpostodds)            # take the exponential
  }else{    
   PO2 <- exp(logBF)                  # POs if model prior is uniform
  }
  resa <- rbind(resa,c(NA,1/PO2))     # add actual POs as additional row
  rownames(resa)[ncolres+1] <- "PO2"  # and row name
 }                                    # print the result
                                      # number of significant digits
 print.default(format(resa, digits = digits), # format the result
    print.gap = 2L,                   # space between adjacent columns
    quote = FALSE,...)                # do not add surrounding quotes
}                                     # end function print.pep

#' (Point) Prediction under PEP approach
#'
#' Computes predicted or fitted values under the PEP approach. Predictions
#' can be based on Bayesian model averaging, maximum a posteriori model or
#' median probability model. For the Bayesian model averaging, a subset of the 
#' top models (either based on explicit number or on their cumulative probability) 
#' can be used for prediction.
#'
#' @param object An object of class pep (e.g., output of \code{pep.lm}).
#' @param xnew An optional data frame of numeric, the new data 
#' on the explanatory variables to be used 
#' for prediction. The data frame needs to contain information about all 
#' explanatory variables available in the full model; if not an error message
#' is output. If omitted, fitted values are 
#' computed.
#' @param estimator A character, the type of prediction. One of 
#' ``BMA'' (Bayesian model averaging, default), 
#' ``MAP'' (maximum a posteriori model) or ``MPM'' (median probability model).
#' Default value=\code{"BMA"}.
#' @param n.models Positive integer, the number of (top) models that
#' prediction is based on or \code{NULL}. Relevant for \code{estimator="BMA"}.
#' Default value=\code{NULL}.
#' @param cumul.prob Numeric between zero and one, cumulative probability of
#' top models to be used for prediction. Relevant for \code{estimator="BMA"}. 
#' Default value=0.99.
#' @param ... Additional parameters to be passed, currently none.
#'
#' @return \code{predict} returns a vector with the predicted (or fitted)
#' values for the different observations.
#'
#' @details
#' When \code{xnew} is missing or \code{xnew} is equal to the initial
#' data frame used for fitting, then fitted 
#' values are computed (and returned).
#'
#' For prediction, Equation 9 of Fouskakis and Ntzoufras (2020) is used.
#'
#' The case of missing data (i.e., presence of NA’s) 
#' and non--quantitative data in the new data frame
#' \code{xnew} is not currently supported.
#'
#' Let \eqn{k} be the number of models with cumulative posterior probability up 
#' to the given value of \code{cumul.prob}. Then, for Bayesian model averaging 
#' the prediction is based on the top \eqn{(k+1)} models if they exist, otherwise
#' on the top \eqn{k} models.
#'
#' When both \code{n.models} and \code{cumul.prob} are provided --- once 
#' specifying the number of models for the given cumulative probability as 
#' described above --- the minimum between the two numbers is used for prediction.
#'
#' @references Fouskakis, D. and Ntzoufras, I. (2022) Power--Expected--Posterior 
#' Priors as Mixtures of g--Priors in Normal Linear Models. 
#' Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
#'
#' Fouskakis, D. and Ntzoufras, I. (2020) Bayesian Model Averaging Using 
#' Power--Expected--Posterior Priors. 
#' Econometrics, 8(2): 17. \doi{10.3390/econometrics8020017}
#'
#' @examples
#' data(UScrime_data)
#' X <- UScrime_data[,-15]
#' set.seed(123)
#' res <- pep.lm(y~.,data=UScrime_data[1:45,],intrinsic=TRUE,
#'                algorithmic.choice="MC3",itermc3=4000)
#' resf <- predict(res)
#' resf2 <- predict(res,estimator="MPM")
#' resp <- predict(res,xnew=X[46:47,])
#' @method predict pep
#' @export
predict.pep <- function(object,xnew,estimator="BMA",
                        n.models=NULL,cumul.prob=0.99,...){
# Takes an object of class pep and (possibly) a new data matrix xnew and
# computes fitted or predicted values. The type of prediction (estimator) can 
# be one of "BMA" (default value), "MAP" or "MPM". For "BMA", a subset of the 
# (top) models can also be used- either giving directly their number through the 
# argument n.models (default value - NULL) or indirectly through their 
# cumulative probability cumul.prob (0.99 by default).
#
 if(!inherits(object,"pep"))          # check if x is an object of type pep
   stop("object should be an object of type pep.")
 x <- object$x                        # input matrix
 y <- object$y                        # response vector
 if (!is.null(x)){
  if(missing(xnew))                   # if new data matrix is missing
   xnew <- x else{                    # use the original input matrix
   if(!is.data.frame(xnew)){          # if new data is not a data frame
    message("Error: xnew should be a data frame.\n") 
    stop("Please respecify and call the function again.")
   }  
   tt <- terms(object$fullmodel)
   Terms <- delete.response(tt)
   xnew <- as.matrix(model.frame(Terms, xnew,na.action=NULL))
  }
                                      # check additional arguments
  checkinputvartypepredict(x,xnew,estimator,n.models,cumul.prob) 
  p <- ncol(x)                        # number of explanatory variables
  intrinsic <- object$intrinsic       # whether intrinsic or PEP prior
  reference.prior <- object$reference.prior   # hyperparameters d0 & d1
  modelsmr <- object$models           # model info
 }
 if (is.null(x)){                     # if no input variable
                                      # the predicted value will be the mean response
  res2M <- matrix(rep(mean(y),length(y)),ncol=1)
  res2v <- as.vector(res2M)           # and turn it to a vector
 }else{
  gamma <- modelsmr[,1:p]             # models as inclusion indicator variables 
  if (p==1)                           # if only one input variable                     
   gamma <- matrix(gamma,ncol=1)      # turn the vector to a 1-column matrix
  if (dim(modelsmr)[1]==1)            # if only one model returned (e.g. from MC3)
                                      # turn the vector to a 1-row matrix
  gamma <- matrix(gamma,nrow=1)
                                      # center input matrix
  x <- scale(x, center=TRUE, scale = FALSE) 
                                      # use the same centering for the new data matrix
  xnew <- scale(xnew, center=attr(x,"scaled:center"), scale = FALSE)
  if (estimator=="BMA"){              # prediction under BMA
   post.prob <- modelsmr[,"Post.Prob"]# posterior probability of models
   cumprob <- cumsum(post.prob)       # cumulative posterior probability
                      # how many models have cumulative prob <= to the given one
                      # including the following model
   indcumprob <- which(cumprob>cumul.prob)[1]
   if (is.na(indcumprob))             # if following model does not exist (e.g.
                                      # when cumul.prob=1)
      indcumprob <- length(cumprob)   # take the available models
   if (!is.null(n.models))            # if the number of models is given as well
     if(n.models<indcumprob)          # then find the min between that and the
                                      # nbr of models with the given cumul.prob
       indcumprob <- n.models
   if (indcumprob>1){ # if more than one models will be used for prediction
    gamma <- gamma[1:indcumprob,]     # get those models
                                      # and normalize their posterior prob
    post.prob <- post.prob[1:indcumprob]
                                      # so that it sums to 1
    post.probn <- post.prob/sum(post.prob) 
    if (p==1)                         # if only 1 input variable     
     gamma <- matrix(gamma,ncol=1)    # turn the vector to a 1-column matrix
   }else{
     gamma <- matrix(gamma[1,],nrow=1)# turn the vector (top model inclusion indicators) 
                                      # to a 1-row matrix
     post.probn <- 1                  # set its normalized posterior prob to 1
   }                                  # get list of fitted/predicted values 
                                      # under each model (code in C++)
   res2 <- predict_pepc(x, gamma, y, xnew, intrinsic, reference.prior)
   res2M <- do.call(cbind,res2)       # turn the list to a matrix
   res2M <- res2M%*%post.probn        # fitted/predicted values under BMA
   res2v <- as.vector(res2M)          # and turn it to a vector
  }else if (estimator=="MAP"){        # prediction using MAP model
   gamma <- matrix(gamma[1,],nrow=1)  # use only top model for prediction
                                      # make the prediction (code in C++)
   res2 <- predict_pepc(x, gamma, y, xnew, intrinsic, reference.prior)
   res2M <- do.call(cbind,res2)       # turn the list to a matrix
   res2v <- as.vector(res2M)          # and then turn it to a vector
  }else if (estimator=="MPM"){        # prediction using median probability model
   inc.probs <- object$inc.probs      # inclusion probs
   indv <- which(inc.probs>0.5)       # input variables with inclusion prob>0.5
   # if (length(indv)==0){              # if there is none
   #                                   # print a message
   #  stop("No input variable with incl prob > 0.5")
   # }else{
    gammav <- rep(0,p)                # vector with inclusion indicators initialization
    gammav[indv] <- 1                 # set 1 to the included variables
    gamma <- matrix(gammav,nrow=1)    # turn the vector to a 1-row matrix 
                                      # make the prediction (code in C++)
    res2 <- predict_pepc(x, gamma, y, xnew, intrinsic, reference.prior)
    res2M <- do.call(cbind,res2)      # turn the list to a matrix
    res2v <- as.vector(res2M)         # and then to a vector
   # }                  
  }                                   # end else (MPM)     
 }                                    # end else (for if (is.null(x)))
 return(res2v)
}                                     # end function predict.pep

#' Model averaged estimates
#'
#' Simulates values from the (joint) posterior distribution of the 
#' beta coefficients under Bayesian model averaging.
#' 
#' @param object An object of class pep (e.g., output of \code{pep.lm}).
#' @param ssize Positive integer, the number of values to be simulated from
#' the (joint) posterior distribution of the beta coefficients.
#' Default value=10000.
#' @param estimator A character, the type of estimation. One of 
#' ``BMA'' (Bayesian model averaging, default), 
#' ``MAP'' (maximum a posteriori model) or ``MPM'' (median probability model).
#' Default value=\code{"BMA"}.
#' @param n.models Positive integer, the number of (top) models where
#' the average is based on or \code{NULL}. Relevant for \code{estimator="BMA"}.
#' Default value=\code{NULL}.
#' @param cumul.prob Numeric between zero and one, cumulative probability of
#' top models to be used for computing the average. Relevant for \code{estimator="BMA"}. 
#' Default value=0.99.
#'
#' @return \code{estimation.pep} returns a matrix (of dimension 
#' \code{ssize} \eqn{\times \, (p+1)}) --- 
#' where the rows correspond
#' to the simulations and the columns to the beta coefficients
#' (including the intercept) --- containing the 
#' simulated data.
#'
#' @details
#' For the computations, Equation 10 of Garcia--Donato and Forte (2018) 
#' is used. That (simplified) formula arises when changing the prior on the
#' model parameters to the reference prior. This change of prior is
#' justified in Garcia--Donato and Forte (2018). The resulting formula is a mixture
#' distribution and the simulation is implemented as follows: firstly the 
#' model (component) based on its posterior probability is chosen and 
#' subsequently the values of the beta coefficients included in the chosen model are
#' drawn from the corresponding multivariate Student distribution, while the
#' values of the beta coefficents outside the chosen model are set to zero.
#'
#' Let \eqn{k} be the number of models with cumulative posterior probability up 
#' to the given value of \code{cumul.prob}. Then, for Bayesian model averaging 
#' the summation is based on the top \eqn{(k+1)} models if they exist, otherwise
#' on the top \eqn{k} models.
#'
#' When both \code{n.models} and \code{cumul.prob} are provided --- once 
#' specifying the number of models for the given cumulative probability as 
#' described above --- the minimum between the two numbers is used for estimation.
#'
#' @references Garcia--Donato, G. and Forte, A. (2018) Bayesian Testing, 
#' Variable Selection and Model Averaging in Linear Models using R with 
#' BayesVarSel. The R Journal, 10(1): 155–174. 
#' \doi{10.32614/RJ-2018-021}
#'
#' @examples
#' data(UScrime_data)
#' res <- pep.lm(y~.,data=UScrime_data)
#' set.seed(123)
#' estM1 <- estimation.pep(res,ssize=2000)
#' estM2 <- estimation.pep(res,ssize=2000,estimator="MPM")
#' @export
estimation.pep <- function(object,ssize=10000,estimator="BMA",
                        n.models=NULL,cumul.prob=0.99){
if(!inherits(object,"pep"))          # check if x is an object of type pep
   stop("object should be an object of type pep.")
if (!((ssize%%1==0)&&(ssize>0)))
   message("Error: the argument ssize should be a positive integer.\n")
x <- object$x                        # input matrix
checkinputvartypepredict(x,x,estimator,n.models,cumul.prob) 
p <- ncol(x)                         # number of explanatory variables
n <- nrow(x)                         # sample size
y <- object$y                        # response vector
modelsmr <- object$models            # models and associated results
if (!is.null(x)){                    # if p>=1
 gamma <- modelsmr[,1:p]
 if (p==1)                           # special case 1: if only one explanatory
                                     # variable
  gamma <- matrix(gamma,ncol=1)
 if (dim(modelsmr)[1]==1)            # special case 2: if only one model
                                     # returned (e.g. from MC3)
  gamma <- matrix(gamma,nrow=1)
 if (estimator=="BMA"){               # estimation under BMA
   post.prob <- modelsmr[,"Post.Prob"]# posterior probability of models
   cumprob <- cumsum(post.prob)       # cumulative posterior probability
                      # how many models have cumulative prob <= to the given one
                      # including the following model
   indcumprob <- which(cumprob>cumul.prob)[1]
   if (is.na(indcumprob))             # if following model does not exist (e.g.
                                      # when cumul.prob=1)
      indcumprob <- length(cumprob)   # take the available models
   if (!is.null(n.models))            # if the number of models is given as well
     if(n.models<indcumprob)          # then find the min between that and the
                                      # nbr of models with the given cumul.prob
       indcumprob <- n.models
   if (indcumprob>1){ # if more than one models will be used for prediction
                                      # and normalize their posterior prob
    post.prob <- post.prob[1:indcumprob]
                                      # so that it sums to 1
    post.probn <- post.prob/sum(post.prob) 
   }else{
     post.probn <- 1                  # set its normalized posterior prob to 1
   }
 }else if (estimator=="MAP"){         # MAP is equivalent to having n.models=1,
                                      # i.e. top model
   indcumprob <- 1
   post.probn <- 1
 }else if (estimator=="MPM"){         # MPM is also 1 model but need to find it
   indcumprob <- 1
   post.probn <- 1
   inc.probs <- object$inc.probs      # inclusion probs
   indv <- which(inc.probs>0.5)       # input variables with inclusion prob>0.5
   gammav <- rep(0,p)                 # vector with inclusion indicators initialization
   if (length(indv)>0)
    gammav[indv] <- 1                 # set 1 to the included variables   
   gamma <- matrix(gammav,nrow=1)     # turn the vector to a 1-row matrix
 }                                    # center input data matrix
   x <- scale(x, center=TRUE, scale = FALSE)
   rsample <- sample(1:indcumprob,ssize,replace=T,prob=post.probn)
   freqtable <- table(rsample)
   cumfreqtable <- cumsum(freqtable)
   cumfreqtableaug <- c(0,cumfreqtable)
   lg <- length(freqtable)
   rsunique <- as.numeric(names(freqtable))   
   ressmpls <- matrix(0,nrow=ssize,ncol=(p+1))
   for(i in 1:lg){
     inds <- which(gamma[rsunique[i],]==1)
     xm <- cbind(1,x[,inds])
     lmfit <- lm.fit(xm,y)
     betaest <- lmfit$coefficients
     sse <- sum((lmfit$residuals)^2)
     dfs <- n-length(inds)-1
     meanv <- betaest
     # sigmam <- solve(t(xm)%*%xm)*sse/dfs  # unstable
     # replace by the following commands
     # see https://genomicsclass.github.io/book/pages/qr_and_regression.html
     # and BayesVarSel R code
     Rinv <- qr.solve(qr.R(lmfit$qr))
     iXtX <- Rinv %*% t(Rinv)
     sigmam <- iXtX*sse/dfs
     ressmpls[(cumfreqtableaug[i]+1):cumfreqtableaug[i+1],c(1,inds+1)] <- 
            mvtnorm::rmvt(n=freqtable[i],sigma=sigmam,df=dfs,delta=meanv,type="shifted")
   } 
   colnames(ressmpls) <- c("Inter", colnames(x))
}else{                            # if no explanatory variable
   ressmpls <- matrix(0,nrow=ssize,ncol=1)
   n <- length(y)
   xm <- matrix(1,nrow=n,ncol=1)  # a column of 1's
   lmfit <- lm.fit(xm,y)          # fit intercept-only model
   betaest <- lmfit$coefficients  # b0 hat
   sse <- sum((lmfit$residuals)^2)# SSE
   dfs <- n-1        # degress of freedom
   meanv <- betaest               # mean
   sigmam <- solve(t(xm)%*%xm)*sse/dfs # sigma value
                                  # simulate from corresponding Student
   ressmpls[,1] <- mvtnorm::rmvt(n=ssize,sigma=sigmam,df=dfs,delta=meanv,type="shifted")
   colnames(ressmpls) <- "Inter"  # add column name
}
   return(ressmpls)
}

#' Posterior predictive distribution under Bayesian model averaging
#'
#' Simulates values from the posterior predictive distribution under 
#' Bayesian model averaging.
#'
#' @param object An object of class pep (e.g., output of \code{pep.lm}).
#' @param xnew An optional data frame of numeric, the new data 
#' on the explanatory variables to be used 
#' for prediction. The data frame needs to contain information about all 
#' explanatory variables available in the full model; if not an error message
#' is output.
#' If omitted, the data frame employed for fitting the full model is used.
#' @param ssize Positive integer, the number of values to be simulated from
#' each posterior predictive distribution.
#' Default value=10000.
#' @param estimator A character, the type of prediction. One of 
#' ``BMA'' (Bayesian model averaging, default), 
#' ``MAP'' (maximum a posteriori model) or ``MPM'' (median probability model).
#' Default value=\code{"BMA"}.
#' @param n.models Positive integer, the number of (top) models where
#' the average is based on or \code{NULL}. Relevant for \code{estimator="BMA"}.
#' Default value=\code{NULL}.
#' @param cumul.prob Numeric between zero and one, cumulative probability of
#' top models to be used for computing the average. Relevant for \code{estimator="BMA"}. 
#' Default value=0.99.
#'
#' @return \code{posteriorpredictive.pep} returns a matrix (of dimension 
#' \code{ssize} \eqn{\times} \code{nrow(xnew)}) --- containing the 
#' simulated data. More specifically, column \eqn{i} contains the simulated
#' values from the posterior predictive corresponding to the \eqn{i}--th new 
#' observation (i.e., \eqn{i}--th row of \code{xnew}).
#'
#' @details
#' For the computations, Equation 11 of Garcia--Donato and Forte (2018) 
#' is used. That (simplified) formula arises when changing the prior on the
#' model parameters to the reference prior. This change of prior is
#' justified in Garcia--Donato and Forte (2018). The resulting formula is a mixture
#' distribution and the simulation is implemented as follows: firstly the 
#' model (component) based on its posterior probability is chosen and 
#' subsequently the value for the response is
#' drawn from the corresponding Student distribution.
#'
#' The case of missing data (i.e., presence of NA’s) and non--quantitative data
#' in the new data frame
#' \code{xnew} is not currently supported.
#'
#' Let \eqn{k} be the number of models with cumulative posterior probability up 
#' to the given value of \code{cumul.prob}. Then, for Bayesian model averaging 
#' the prediction is based on the top \eqn{(k+1)} models if they exist, otherwise
#' on the top \eqn{k} models.
#'
#' When both \code{n.models} and \code{cumul.prob} are provided --- once 
#' specifying the number of models for the given cumulative probability as 
#' described above --- the minimum between the two numbers is used for prediction.
#'
#' @references Garcia--Donato, G. and Forte, A. (2018) Bayesian Testing, 
#' Variable Selection and Model Averaging in Linear Models using R with 
#' BayesVarSel. The R Journal, 10(1): 155–174. 
#' \doi{10.32614/RJ-2018-021}
#'
#' @examples
#' data(UScrime_data)
#' X <- UScrime_data[,-15]
#' set.seed(123)
#' res <- pep.lm(y~.,data=UScrime_data[1:45,],intrinsic=TRUE,
#'                algorithmic.choice="MC3",itermc3=4000)
#' resf <- posteriorpredictive.pep(res,ssize=2000,n.models=5)
#' resf2 <- posteriorpredictive.pep(res,ssize=2000,estimator="MPM")
#' resp <- posteriorpredictive.pep(res,xnew=X[46:47,],ssize=2000,n.models=5)
#' @importFrom stats delete.response
#' @export
posteriorpredictive.pep <- function(object,xnew,ssize=10000,estimator="BMA",
                                    n.models=NULL,cumul.prob=0.99){
if(!inherits(object,"pep"))          # check if x is an object of type pep
   stop("object should be an object of type pep.")
if (!((ssize%%1==0)&&(ssize>0)))
   message("Error: the argument ssize should be a positive integer.\n")
x <- object$x                        # input matrix
y <- object$y                        # response vector
if (!is.null(x)){                    # if p>=1
 if(missing(xnew))                   # if new data matrix is missing
  xnew <- x else{                    # use the original input matrix
  if(!is.data.frame(xnew)){          # if new data is not a data frame
   message("Error: xnew should be a data frame.\n") 
   stop("Please respecify and call the function again.")
  }
  tt <- terms(object$fullmodel)
  Terms <- delete.response(tt)
  xnew <- as.matrix(model.frame(Terms, xnew,na.action=NULL))
 }
 checkinputvartypepredict(x,xnew,estimator,n.models,cumul.prob) 
 p <- ncol(x)                         # number of explanatory variables
 n <- nrow(x)                         # sample size
 modelsmr <- object$models
 gamma <- modelsmr[,1:p]
 if (p==1)                           # special case 1: if only one explanatory
                                     # variable
  gamma <- matrix(gamma,ncol=1)
 if (dim(modelsmr)[1]==1)            # special case 2: if only one model
                                     # returned (e.g. from MC3)
  gamma <- matrix(gamma,nrow=1)
 if (estimator=="BMA"){               # estimation under BMA
   post.prob <- modelsmr[,"Post.Prob"]# posterior probability of models
   cumprob <- cumsum(post.prob)       # cumulative posterior probability
                      # how many models have cumulative prob <= to the given one
                      # including the following model
   indcumprob <- which(cumprob>cumul.prob)[1]
   if (is.na(indcumprob))             # if following model does not exist (e.g.
                                      # when cumul.prob=1)
      indcumprob <- length(cumprob)   # take the available models
   if (!is.null(n.models))            # if the number of models is given as well
     if(n.models<indcumprob)          # then find the min between that and the
                                      # nbr of models with the given cumul.prob
       indcumprob <- n.models
   if (indcumprob>1){ # if more than one models will be used for prediction
                                      # and normalize their posterior prob
    post.prob <- post.prob[1:indcumprob]
                                      # so that it sums to 1
    post.probn <- post.prob/sum(post.prob) 
   }else{
     post.probn <- 1                  # set its normalized posterior prob to 1
   }
 }else if (estimator=="MAP"){         # MAP is equivalent to having n.models=1,
                                      # i.e. top model
   indcumprob <- 1
   post.probn <- 1
 }else if (estimator=="MPM"){         # MPM is also 1 model but need to find it
   indcumprob <- 1
   post.probn <- 1
   inc.probs <- object$inc.probs      # inclusion probs
   indv <- which(inc.probs>0.5)       # input variables with inclusion prob>0.5
   gammav <- rep(0,p)                 # vector with inclusion indicators initialization
   if (length(indv)>0)
    gammav[indv] <- 1                 # set 1 to the included variables   
   gamma <- matrix(gammav,nrow=1)     # turn the vector to a 1-row matrix
 }
  x <- scale(x, center=TRUE, scale = FALSE)
                                      # use the same centering for the new data matrix
  xnew <- scale(xnew, center=attr(x,"scaled:center"), scale = FALSE)
  ressmpls <- matrix(0,nrow=ssize,ncol=nrow(xnew))
  for(k in 1:nrow(xnew)){
   xv <- xnew[k,]
   rsample <- sample(1:indcumprob,ssize,replace=T,prob=post.probn)
   freqtable <- table(rsample)
   cumfreqtable <- cumsum(freqtable)
   cumfreqtableaug <- c(0,cumfreqtable)
   lg <- length(freqtable)
   rsunique <- as.numeric(names(freqtable))      
   for(i in 1:lg){
     inds <- which(gamma[rsunique[i],]==1)
     xm <- cbind(1,x[,inds])
     lmfit <- lm.fit(xm,y)
     betaest <- lmfit$coefficients
     sse <- sum((lmfit$residuals)^2)
     dfs <- n-length(inds)-1
     xvm <- c(1,xv[inds])
     meanv <- (xvm%*%betaest)[1]
     # hgamma <- 1-(xvm%*%solve(xvm%*%t(xvm)+t(xm)%*%xm)%*%xvm)[1]
     # same here - replace by the following (more stable) commands
     Xstar <- rbind(xm, xvm)
     qrXstar <- qr(Xstar)
     Rinv <- qr.solve(qr.R(qrXstar))
     iXstartXstar <- Rinv %*% t(Rinv)
     hgamma <- 1-(xvm%*%iXstartXstar%*%xvm)[1]
     sigmam <- matrix(sse/(dfs*hgamma),nrow=1)
     ressmpls[(cumfreqtableaug[i]+1):cumfreqtableaug[i+1],k] <- 
            mvtnorm::rmvt(n=freqtable[i],sigma=sigmam,df=dfs,delta=meanv,type="shifted")
   } 
  }
 }else{                            # if no explanatory variable
     ressmpls <- matrix(0,nrow=ssize,ncol=1)
     n <- length(y)
     xm <- matrix(1,nrow=n,ncol=1)  # a column of 1's
     lmfit <- lm.fit(xm,y)          # fit intercept-only model
     betaest <- lmfit$coefficients  # b0 hat
     sse <- sum((lmfit$residuals)^2)# SSE
     dfs <- n-1                     # degress of freedom
     xvm <- 1
     meanv <- (xvm%*%betaest)[1]
     hgamma <- 1-(xvm%*%solve(xvm%*%t(xvm)+t(xm)%*%xm)%*%xvm)[1]
     sigmam <- matrix(sse/(dfs*hgamma),nrow=1)
     ressmpls[,1] <- 
            mvtnorm::rmvt(n=ssize,sigma=sigmam,df=dfs,delta=meanv,type="shifted")
    
  }
  return(ressmpls)
}

#' Heatmap for top models
#'
#' Generates a heatmap where the rows correspond to the (top) models
#' and the columns to the input/explanatory variables. The value depicted
#' in cell \eqn{(i,j)} corresponds to the posterior inclusion probability 
#' of variable \eqn{i} if this
#' is included in model \eqn{j} and zero otherwise.
#'
#' @param x An object of class pep (e.g., output of \code{pep.lm}).
#' @param n.models Positive integer, number of models to be shown on the 
#' heatmap. Default value=20.
#' @param ... Additional parameters to be passed to \code{heatmap}.
#'
#' @return No return value, used for heatmap generation.
#'
#' @details
#' The number of models to be displayed on the heatmap is computed as the minimum
#' between the number asked by the user and the number of models present in
#' the object \code{x}.
#'
#' The color code is as follows: the darker the blue in the figure, the 
#' higher the posterior inclusion probability is, while white means that the 
#' variable is not included in the model.
#'
#' In the special case of no explanatory variables, no heatmap is produced
#' and a message is printed.
#'
#' @examples
#' data(UScrime_data)
#' set.seed(123)
#' resu <- pep.lm(y~.,data=UScrime_data,beta.binom=FALSE,
#'                algorithmic.choice="MC3",itermc3=5000)
#' image(resu)
#' image(resu,n.models=10)
#'
#' @seealso \code{\link{plot.pep}}
#' @importFrom grDevices colorRampPalette
#' @importFrom stats heatmap
#' @method image pep
#' @export
image.pep <- function(x,n.models=20,...){
# Takes the result of full enumeration or mc3 algorithm res
# and the number of top models n.models (20 by default) we are interested in 
# and generates a heatmap of a matrix where the rows correspond to the top 
# models, the columns to the explanatory variables while the matrix value (i,j) is the
# inclusion probability of variable j in case it is included in model i or
# 0 otherwise.
#
 if(!inherits(x,"pep"))               # check if x is an object of type pep
   stop("The input to make the image plot should be an object of type pep.")
                                      # check if n.models is a positive integer
 if (!((n.models%%1==0)&&(n.models>0))){
   stop("The argument n.models should be a positive integer.")
 }
 inc.probs <- x$inc.probs             # inclusion probs
 p <- length(inc.probs)               # number of explanatory variables
 modelsmr <- x$models
                                      # take the min between the 'asked' number
                                      # of models and the available
 n.models <- min(n.models,dim(modelsmr)[1])
 if (p==0)                            # if no input variables
  stop("An image plot cannot be generated without explanatory variables.")

 Mi <- modelsmr[1:n.models,1:p]       # matrix of variable inclusion 
                                      # indicators for the top n.models
 if(n.models==1)                      # if only 1 model will be shown
  Mi <- matrix(Mi,nrow=1)             # turn the vector into a 1-row matrix
 Misc <- t(t(Mi)*inc.probs)           # replace the 1's in the previous matrix
                                      # by the actual inclusion probs
 rownames(Misc) <- 1:n.models         # add row names
                                      # color code
 colr <- c("white",colorRampPalette(c("lightblue","blue","darkblue"))(256))
 # The next conditions have been introduced for 
 # visualization purposes 
 if((p==1)&&(n.models>1)){            # if only 1 input variable (but several models)
   Misc <- cbind(Misc,0)              # add a second column with 0's
                                      # column names: X1 & NA
   colnames(Misc) <- c(names(inc.probs),NA)
 }  
 if((n.models==1)&&(p>1)){            # if only 1 model (but several input variables)
   if (sum(Misc[1,])==0)              # if this top model is the null model
    colr <- "white"                   # use only white (because everything will be 0)
   Misc <- rbind(Misc,0)              # add a second row with 0's
   rownames(Misc)[2] <- NA            # with rowname NA
   n.models <- 2                      # set the number of models to 2
 }
 if((n.models==1)&&(p==1)){           # if only 1 model and only 1 input variable
   if (Misc[1,1]==0)                  # if this top model is the null model
    colr <- "white"                   # use only white (because everything will be 0)
   Misc <- rbind(Misc,0)              # add a second row with 0
   Misc <- cbind(Misc,0)              # add a second column with 0's
   rownames(Misc) <- c(1,NA)          # row names: 1 & NA 
                                      # column names: X1 & NA    
   colnames(Misc) <- c(names(inc.probs),NA)
   n.models <- 2                      # set the number of models to 2                      
 }
                                      # generate heatmap
 heatmap(Misc[n.models:1,],           # invert row order so that the top model
                                      # appears on the top of the plot
         Rowv=NA,Colv=NA,             # without row or column clustering/
                                      # dendrogram and 
         scale="none",                # without scaling the matrix
         col = colr,...)
}                                     # end function image.pep

#' Bayes factor for model comparison
#'
#' Given two models to be compared (the one nested to the other), 
#' computes the corresponding Bayes factor.
#'
#' @param formula1 One of the two formulas/models to be compared. 
#' @param formula2 The second formula/model. The one model 
#' needs to be nested to the other.
#' @param data A data frame (of numeric values), containing the data.
#' @param intrinsic Logical, indicating whether the PEP 
#' (\code{FALSE}) or the intrinsic --- which   
#' is a special case of it --- (\code{TRUE}) should be used as prior on the  
#' regression parameters. Default value=\code{FALSE}.
#' @param reference.prior Logical, indicating whether the reference prior
#' (\code{TRUE}) or the dependence Jeffreys prior (\code{FALSE}) is used as 
#' baseline. Default value=\code{TRUE}.
#'
#' @return \code{peptest} returns the Bayes factor, i.e., a numeric value.
#' For the ratio, the marginal likelihood of the more complex model (nominator)
#' with respect to that of the simpler one (denominator) is computed. 
#' Both marginal likelihoods are computed with respect
#' to the intercept--only model (reference model).
#'
#' @details 
#' This function can be used to perform hypothesis testing indirectly. More
#' specifically, for the interpretation of the result (Bayes factor), the table in 
#' Kass and Raftery (1995) can be used.
#'
#' The function works when \eqn{p\leq n-2}, where \eqn{p} is the number of explanatory 
#' variables in the more complex model and \eqn{n} is the sample size.
#'
#' The case of missing data (i.e., presence of \code{NA}'s either in the  
#' data matrix corresponding to the explanatory variables of the more complex 
#' model or the response vector) is not currently supported. Further, the
#' explanatory variables of the more complex model need to be quantitative.
#'
#' If \eqn{p>1}, the explanatory variables of the more complex model 
#' cannot have an exact linear relationship (perfect multicollinearity).
#'
#' @references Kass, R. and Raftery, A. (1995) Bayes Factors. 
#' Journal of the American Statistical Association, 90(430): 773–795. 
#' \doi{10.1080/01621459.1995.10476572}
#'
#' @examples
#' data(UScrime_data)
#' resBF1 <- peptest(y~1,y~M+Ed,UScrime_data)
#' resBF1i <- peptest(y~1,y~M+Ed,UScrime_data, intrinsic=TRUE)
#' resBF2j <- peptest(y~M+Ed+Po1+Po2,y~M+Ed,UScrime_data,
#'                    reference.prior=FALSE)
#' resBF2ij <- peptest(y~M+Ed+Po1+Po2,y~M+Ed,UScrime_data,
#'                     intrinsic=TRUE, reference.prior=FALSE)
#'
#' @export
peptest <- function(formula1, formula2, data, intrinsic=FALSE, 
                    reference.prior=TRUE){
formulalist <- list(formula1, formula2)
# check arguments' type
if(!((is.list(formulalist))&&(length(formulalist)==2))){
  message("Error: formulalist should be a list with 2 components.\n") 
  stop("Please respecify and call the function again.")
}
for (j in 1:2)
  if(!inherits(formulalist[[j]],"formula"))    # check if formula is of type formula
      stop("formula1 and formula2 should be of type formula.")

if(!is.data.frame(data)){
  message("Error: data should be a data frame.\n") 
  stop("Please respecify and call the function again.")
}
predctrslist <- lapply(formulalist, 
                       function(x){attr(terms(x,data=data),"term.labels")})
predctrsnbr <- unlist(lapply(predctrslist,length))
simplemodind <- which.min(predctrsnbr)
simplemodpredctrs <- predctrslist[[simplemodind]]
compositemodpredctrs <- predctrslist[[-simplemodind]]
isnestedmod <- all(simplemodpredctrs%in%compositemodpredctrs)
if(!(isnestedmod==TRUE)){
  message("Error: one of the given models should be nested within the other.\n")
  stop("Please respecify and call the function again.")
}

# Reference model
k0 <- 1
R0 <- 0

# Remaining models
diffpredctrs <- compositemodpredctrs
ndiffpredctrs <- length(diffpredctrs)
gamma <- matrix(0,nrow=2,ncol=ndiffpredctrs)
colnames(gamma) <- diffpredctrs
for (i in 1:2)
  gamma[i,predctrslist[[i]]] <- 1
mf <- model.frame(formulalist[[-simplemodind]],data,na.action=NULL)
y <- model.response(mf, "numeric")
predctrs <- attr(terms(formulalist[[-simplemodind]],data=data),"term.labels")
x <- as.matrix(mf[,predctrs])

# check that the function arguments are of the correct type
checkinputvartype(x,y,intrinsic,reference.prior,
                  TRUE,FALSE)
# check that n>=p+2 for feasibility
compareparamvssamplesize(x)
checkmissingness(x,y)                 # check that there are no NAs 
p <- ncol(x)                          # number of explanatory variables
if (p>1)                              # if more than one
 checkmatrixrank(x)

if (length(simplemodpredctrs)==0){
  gamma <- matrix(gamma[-simplemodind,],nrow=1)
  res1 <- test_pepc(x,gamma,y,intrinsic=intrinsic,reference.prior,k0,R0)$marglikel
  finres <- exp(res1)
  names(finres) <- "H1 model to H0 model"
}else{
  res1 <- test_pepc(x,gamma,y,intrinsic=intrinsic,reference.prior,k0,R0)$marglikel
  finres <- exp(res1[-simplemodind]-res1[simplemodind])
  names(finres) <- "H1 model to H0 model"
}
return(finres)
}

# Constructor of the class pep
pep <- function(models, inc.probs, x, y, fullmodel, mapp, 
                intrinsic, reference.prior, 
                beta.binom, allvisitedmodsM){
 if(!missing(allvisitedmodsM))
  obj <- list(models = models, inc.probs = inc.probs, x = x, y = y, 
              fullmodel = fullmodel, mapp = mapp,
              intrinsic = intrinsic, reference.prior=reference.prior, 
              beta.binom = beta.binom, allvisitedmodsM = allvisitedmodsM)else
  obj <- list(models = models, inc.probs = inc.probs, x = x, y = y, 
              fullmodel = fullmodel, mapp = mapp, 
              intrinsic = intrinsic, reference.prior=reference.prior, 
              beta.binom = beta.binom)
 class(obj) <- "pep"
 obj
}
