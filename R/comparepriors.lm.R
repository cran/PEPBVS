#' Selected models under different choices of prior on the model parameters and 
#' the model space
#'
#' Given a formula and a data frame, computes the maximum a posteriori (MAP) model  
#' and median probability model (MPM) 
#' for different choices of prior on the model parameters and the model space.
#' Normal linear models are assumed for the data with the prior distribution on the model 
#' parameters being one or more of the following: PEP, intrinsic, Zellner’s 
#' \eqn{g}--prior, Zellner and Siow, benchmark, robust, hyper--\eqn{g} and related hyper--\eqn{g}--\eqn{n}.
#' The prior distribution on the model space can be either the uniform on models
#' or the uniform on the model dimension (special case of the beta--binomial prior).
#' The model space consists of all possible models including an intercept term.
#' Model selection is performed by using either full enumeration 
#' and evaluation of all models  
#' (for model spaces of small--to--moderate dimension) or a
#' Markov chain Monte Carlo (MCMC) scheme (for 
#' model spaces of large dimension). 
#'
#' @param formula A formula, defining the full model.
#' @param data A data frame (of numeric values), containing the data. 
#' @param algorithmic.choice A character, the type of algorithm to be used
#' for selection: full enumeration and evaluation of all models or an MCMC scheme. 
#' One of ``automatic'' (the choice is done automatically based on the number
#' of explanatory variables in the full model), ``full enumeration'' 
#' or ``MCMC''. Default value=\code{"automatic"}.
#' @param priorbetacoeff A vector of character containing the different priors on the model
#' parameters. The character can be one of ``PEP'', ``intrinsic'', ``Robust'', ``gZellner'', 
#' ``ZellnerSiow'', ``FLS'', ``hyper--g'' and ``hyper--g--n''. 
#' \cr Default value=\code{
#' c("PEP","intrinsic","Robust", "gZellner","ZellnerSiow",}
#' \code{"FLS","hyper-g","hyper-g-n")},
#' i.e., all supported priors are tested.
#' @param reference.prior A vector of logical indicating the baseline prior that is used for 
#' PEP/intrinsic. It can be TRUE (reference prior is used), FALSE (dependence Jeffreys prior 
#' is used) or both. Default value=\code{c(TRUE,FALSE)}, i.e., both baseline priors are tested.
#' @param priormodels A vector of character containing the different priors on the model
#' space. The character can be one of ``beta--binomial'' and ``uniform''. 
#' \cr Default value=\code{c("beta-binomial","uniform")}, i.e., both supported priors are tested.
#' @param burnin Non--negative integer, the burnin period for the MCMC scheme.
#' Default value=1000.
#' @param itermcmc Positive integer (larger than \code{burnin}),
#' the (total) number of iterations for the MCMC scheme. Default value=11000.
#'
#' @return \code{comparepriors.lm} returns a list with two elements:
#' \item{MAPmodels}{A data frame containing the MAP models for all different combinations of 
#' prior on the model parameters and the model space. In particular, in row \eqn{i} the following
#' information is presented: prior on the model parameters, prior on the model space, 
#' hyperparameter value, MAP model (corresponding to the specific combination 
#' of priors on model parameters and model space) 
#' represented with variable inclusion indicators, and the R package used. When an MCMC scheme
#' has been used, there are two additional columns: one depicting the specific algorithm 
#' that has been used and one with the MC standard error (to assess convergence).
#' With an MCMC scheme, the MAP model output corresponds to the most frequently
#' `visited'.}
#' \item{MPMmodels}{Same as the first element containing the MPM models instead.}
#'
#' @details
#' The different priors on the model parameters are implemented using 
#' different packages: for PEP and intrinsic, the current package is used. 
#' For hyper--\eqn{g} and related hyper--\eqn{g}--n priors, the R package \pkg{BAS} is used. Finally, 
#' for the Zellner’s \eqn{g}--prior (``gZellner''), the Zellner and Siow 
#' (``ZellnerSiow''), the robust and the benchmark (``FLS'') prior,
#' the results are obtained using \pkg{BayesVarSel}.
#' 
#' The prior distribution on the model space can be either the uniform on models
#' or the beta--binomial. For the beta--binomial prior, 
#' the following special case is used: uniform prior on model dimension.
#'
#' When an MCMC scheme is used, the R package \pkg{BAS} uses the birth/death random walk 
#' in Raftery et al. (1997) combined with a random swap move, \pkg{BayesVarSel} 
#' uses Gibbs sampling while \pkg{PEPBVS} implements the MC3 algorithm
#' described in the Appendix of Fouskakis and Ntzoufras (2022). 
#'
#' To assess MCMC convergence, Monte Carlo (MC) standard error is 
#' computed using batch means estimator (implemented in the R package \pkg{mcmcse}).
#' For computing a standard error, the number (\code{itermcmc}-\code{burnin}) 
#' needs to be larger than 100.  
#' This quantity cannot be computed for the cases treated by \pkg{BAS} --- 
#' since all `visited' models are not available in the function output --- and thus for those cases 
#' \code{NA} is depicted in the relevant column instead. 
#'
#' Similar to \code{\link{pep.lm}}, if \code{algorithmic.choice} equals ``automatic'' then 
#' model selection is implemented as follows: if \eqn{p < 20} (where \eqn{p} is the 
#' number of explanatory variables in the full model without the intercept), full enumeration 
#' and evaluation of all models is performed, otherwise an MCMC scheme is used.
#' To avoid potential memory or time constraints, if \code{algorithmic.choice} 
#' equals ``full enumeration'' but \eqn{p \geq 20}, once issuing a warning message,
#' an MCMC scheme is used instead.
#'
#' Similar constraints to \code{\link{pep.lm}} hold for the data, i.e., 
#' the case of missing data is not currently supported, the explanatory 
#' variables need to be quantitative and cannot have an exact linear relationship, 
#' and \eqn{p\leq n-2} (\eqn{n} being the sample size). 
#'
#' @references Bayarri, M., Berger, J., Forte, A. and Garcia--Donato, G. (2012) 
#' Criteria for Bayesian Model Choice with Application to Variable Selection. 
#' The Annals of Statistics, 40(3): 1550–1577. \doi{10.1214/12-AOS1013}
#'
#' Fouskakis, D. and Ntzoufras, I. (2022) Power--Expected--Posterior 
#' Priors as Mixtures of g--Priors in Normal Linear Models. 
#' Bayesian Analysis, 17(4): 1073-1099. \doi{10.1214/21-BA1288}
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
#' @examples
#' data(UScrime_data)
#' resc <- comparepriors.lm(y~.,UScrime_data,
#'                          priorbetacoeff = c("PEP","hyper-g-n"),
#'                          reference.prior = TRUE,priormodels = "beta-binomial")
#'
#' @importFrom BAS bas.lm beta.binomial uniform
#' @importFrom BayesVarSel Bvs GibbsBvs
#' @importFrom mcmcse mcse
#' @export
comparepriors.lm <- function(formula, data,
                             algorithmic.choice = "automatic",
                             priorbetacoeff = c("PEP","intrinsic","Robust",
                                                "gZellner","ZellnerSiow","FLS",
                                                "hyper-g","hyper-g-n"),
                             reference.prior = c(TRUE,FALSE),
                             priormodels = c("beta-binomial","uniform"),
                             burnin=1000, itermcmc=11000){
#library(BAS)
#library(BayesVarSel)
#library(mcmcse)
# check arguments' type
if(!inherits(formula,"formula"))    # check if formula is of type formula
  stop("formula should be of type formula.")
if(!is.data.frame(data)){
  message("Error: data should be a data frame.\n") 
  stop("Please respecify and call the function again.")
}
                      # check if the argument priorbetacoeff is one or more of 
                      # PEP, intrinsic, Robust, gZellner, ZellnerSiow, FLS,
                      # hyper-g, or hyper-g-n
priorbetacoeff <- unique(priorbetacoeff)
if (!(is.character(priorbetacoeff)&&
      all(priorbetacoeff%in%c("PEP","intrinsic","Robust",
                          "gZellner","ZellnerSiow","FLS",
                          "hyper-g","hyper-g-n")))){
  message("Error: priorbetacoeff should be one or more of PEP, intrinsic, 
           Robust, gZellner, ZellnerSiow, FLS, hyper-g, or hyper-g-n.\n")
  stop("Please respecify and call the function again.")
}
                      # check if the argument priormodels is one or more of 
                      # beta-binomial or uniform
priormodels <- unique(priormodels)
if (!(is.character(priormodels)&&
      all(priormodels%in%c("beta-binomial","uniform")))){
  message("Error: priormodels should be one or more of beta-binomial or uniform.\n")
  stop("Please respecify and call the function again.")
}
if (!(is.character(algorithmic.choice)&&
             algorithmic.choice%in%c("automatic","full enumeration","MCMC"))){
  message("Error: algorithmic.choice should be one of automatic,
                  full enumeration or MC3.\n")
  stop("Please respecify and call the function again.")
}
predctrs <- attr(terms(formula,data=data),"term.labels")
p <- length(predctrs)
algorithm <- "full enumeration" 
if(p>=20){
   algorithm <- "MC3"
   if (algorithmic.choice=="full enumeration")
    message("The number of covariates does not allow full enumeration, 
             thus an MCMC scheme will be applied instead.\n")
}
if (algorithmic.choice=="MCMC")
   algorithm <- "MC3"
if(algorithm=="MC3"){
  if (!((burnin%%1==0)&&(burnin>=0)&&(burnin<itermcmc))){
   message("Error: the argument burnin should be a non-negative integer 
        smaller than the total number of iterations.\n")
   stop("Please respecify and call the function again.")
  }
  if (!((itermcmc%%1==0)&&(itermcmc>0))){
   message("Error: the argument itermc3 should be a positive integer.\n")
   stop("Please respecify and call the function again.")
  }
}
# R package PEPBVS
v1pepa <- priorbetacoeff[which(priorbetacoeff%in%c("PEP","intrinsic"))]
v1pep <- !(v1pepa %in% "PEP") 
names(v1pep) <- v1pepa
lg1 <- length(v1pep)
if (lg1>0){
 reference.prior <- unique(reference.prior)
 if (!all(is.logical(reference.prior))){
  message("Error: reference.prior should be boolean: TRUE, FALSE or both.\n")
  stop("Please respecify and call the function again.")
 }
 v2pep <- reference.prior
 names(v2pep)[which(v2pep==TRUE)] <- "reference"
 names(v2pep)[which(v2pep==FALSE)] <- "Jeffreys"
 lg2 <- length(v2pep)
}else{
 v2pep <- NA
 lg2 <- 0
 mf <- model.frame(formula,data,na.action=NULL)
 y <- model.response(mf, "numeric")
 if(p>=1){
  x <- as.matrix(mf[,predctrs])
  checkinputvartype(x,y,TRUE,TRUE,
                  TRUE,FALSE)
  compareparamvssamplesize(x)
  checkmissingness(x,y)
  if (p>1)                              # if more than one
   checkmatrixrank(x)
 }else{
   checkinputvartype(y=y,intrinsic=TRUE,reference.prior=TRUE,
                   beta.binom=TRUE,constant.term=FALSE)
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
}
}
v3pepa <- priormodels
v3pep <- !(v3pepa %in% "uniform") 
names(v3pep) <- v3pepa
lg3 <- length(v3pep)
res1pepMAP <- matrix(0,nrow=lg1*lg2*lg3,ncol=p)
res1pepMPM <- matrix(0,nrow=lg1*lg2*lg3,ncol=p)
mcerrorpepMAP <- mcerrorpepMPM <- vector(length=lg1*lg2*lg3)
if(lg1>0){
 l <- 1
 for(i in 1:lg1)
  for(j in 1:lg2)
   for(k in 1:lg3){
     res1 <- pep.lm(formula, data, algorithmic.choice=algorithm,
                    intrinsic=v1pep[i], reference.prior=v2pep[j], 
                    beta.binom=v3pep[k],burnin=burnin, itermc3=itermcmc)
     if (algorithm=="MC3"){
      Mvisitedmodels <- res1$allvisitedmodsM
      if ((p>=1)&&(dim(Mvisitedmodels)[1]>100)){
       Mauxil <- Mvisitedmodels[,1:p]
       if(p==1)
        Mauxil <- matrix(Mauxil,ncol=1)
       visitedmodelsv <- apply(Mauxil,1,paste,collapse='')
       MAPm <-paste(res1[[1]][1,1:p],collapse='')
       MPMm <-paste(as.numeric(res1$inc.probs > 0.5),collapse='')
       vaux <- ifelse(visitedmodelsv==MAPm,1,0)
       vaux2 <- ifelse(visitedmodelsv==MPMm,1,0)
       mcerrorpepMAP[l] <- paste(round(mcmcse::mcse(vaux)$se,5))
       mcerrorpepMPM[l] <- paste(round(mcmcse::mcse(vaux2)$se,5))
      }else{
        message("The number of Monte Carlo simulations is quite small to 
compute the standard error.\n")
       mcerrorpepMAP[l] <- mcerrorpepMPM[l] <- NA
      }
     }
     res1pepMAP[l,] <- res1[[1]][1,1:p]
     res1pepMPM[l,] <- res1$inc.probs > 0.5
     l <- l+1
   }
}
infoM1 <- data.frame(cbind(rep(names(v1pep),each=lg2*lg3),
                           rep(names(v3pep),lg1*lg2),
                           rep(rep(names(v2pep),each=lg3),lg1)))

# R package BayesVarSel
v1bayesvarsel <- priorbetacoeff[which(priorbetacoeff%in%c("Robust","gZellner",
                                                          "ZellnerSiow","FLS"))]
parameterinfov <- c("a=1/2, b=1, rl=1/kl","g=n","a=1","g=max(n,p^2)")
names(parameterinfov) <- c("Robust","gZellner","ZellnerSiow","FLS")
v2bayesvarsel <- priormodels
names(v2bayesvarsel) <- v2bayesvarsel
v2bayesvarsel[v2bayesvarsel=="uniform"] <- "Constant" 
v2bayesvarsel[v2bayesvarsel=="beta-binomial"] <- "ScottBerger" 
lg1 <- length(v1bayesvarsel)
lg2 <- length(v2bayesvarsel)
res2bayesvarselMAP <- matrix(0,nrow=lg1*lg2,ncol=p)
res2bayesvarselMPM <- matrix(0,nrow=lg1*lg2,ncol=p)
mcerrorbayesvarselMAP <- mcerrorbayesvarselMPM <- vector(length=lg1*lg2)
if(lg1>0){
 l <- 1
 for(i in 1:lg1)
  for(j in 1:lg2){
   if (algorithm=="full enumeration")
    res2 <- BayesVarSel::Bvs(formula, data, prior.betas = v1bayesvarsel[i], 
                prior.models = v2bayesvarsel[j], n.keep=1, time.test=FALSE) else{
    res2 <- BayesVarSel::GibbsBvs(formula, data, prior.betas = v1bayesvarsel[i], 
                prior.models = v2bayesvarsel[j], time.test=FALSE,
                n.burnin = burnin, n.iter = itermcmc - burnin)
    Mvisitedmodels <- res2$modelslogBF
    if((itermcmc - burnin)==1)
     nms <- paste(Mvisitedmodels[,1:p],collapse='') else
     nms <- apply(Mvisitedmodels[,1:p],1,paste,collapse='')
    indmax <- which(nms==names(which.max((table(nms)))))[1]
    res2$HPMbin <- Mvisitedmodels[indmax,1:p]
    if (dim(Mvisitedmodels)[1]>100){
     visitedmodelsv <- apply(Mvisitedmodels[,1:p],1,paste,collapse='')
     MAPm <-paste(res2$HPMbin,collapse='')
     MPMm <-paste(as.numeric(res2$inclprob > 0.5),collapse='')
     vaux <- ifelse(visitedmodelsv==MAPm,1,0)
     vaux2 <- ifelse(visitedmodelsv==MPMm,1,0)
     mcerrorbayesvarselMAP[l] <- paste(round(mcmcse::mcse(vaux)$se,5))
     mcerrorbayesvarselMPM[l] <- paste(round(mcmcse::mcse(vaux2)$se,5))
    }else{
     message("The number of Monte Carlo simulations is quite small to 
compute the standard error.\n")
     mcerrorbayesvarselMAP[l] <- mcerrorbayesvarselMPM[l] <- NA
    }
   }
   res2bayesvarselMAP[l,] <- res2$HPMbin
   res2bayesvarselMPM[l,] <- res2$inclprob > 0.5
   l <- l+1
  }
}
interm<- parameterinfov[rep(v1bayesvarsel,each=lg2)]
names(interm) <- NULL
infoM2 <- data.frame(cbind(rep(v1bayesvarsel,each=lg2),
                           rep(names(v2bayesvarsel),lg1),interm))
names(infoM2) <- names(infoM1)
# R package BAS
v1bas <- priorbetacoeff[which(priorbetacoeff%in%c("hyper-g","hyper-g-n"))]
parameterinfov <- c("a=3","a=3")
names(parameterinfov) <- c("hyper-g","hyper-g-n")
v2bas <- priormodels
lg1 <- length(v1bas)
lg2 <- length(v2bas)
res3basMAP <- matrix(0,nrow=lg1*lg2,ncol=p)
res3basMPM <- matrix(0,nrow=lg1*lg2,ncol=p)
if(lg1>0){
 if (algorithm=="full enumeration")
   algorithm <- "BAS" else
   algorithm <- "MCMC"
 l <- 1
 for(i in 1:lg1)
  for(j in 1:lg2){
   if(v2bas[j]=="uniform"){
    if(algorithm=="BAS"){
     res3 <- BAS::bas.lm(formula, data, prior = v1bas[i], 
                modelprior = BAS::uniform(), method=algorithm, n.models=2^p)
    }else{
     res3 <- BAS::bas.lm(formula, data, prior = v1bas[i], 
                modelprior = BAS::uniform(), method=algorithm,
                burnin.iterations = burnin, MCMC.iterations = itermcmc - burnin, 
                n.models=itermcmc+1)
    }
   }
   else{
    if(algorithm=="BAS"){
     res3 <- BAS::bas.lm(formula, data, prior = v1bas[i], 
                modelprior = BAS::beta.binomial(1,1), method=algorithm, n.models=2^p)
    }else{
     res3 <- BAS::bas.lm(formula, data, prior = v1bas[i], 
                modelprior = BAS::beta.binomial(1,1), method=algorithm,
                burnin.iterations = burnin, MCMC.iterations = itermcmc - burnin, 
                n.models=itermcmc+1)
    }
   }
   res3basMAP[l,] <- summary(res3)[2:(p+1),2]
   res3basMPM[l,] <- res3$probne0[-1] > 0.5
   l <- l+1
  }
}
interm <- parameterinfov[rep(v1bas,each=lg2)]
names(interm) <- NULL
infoM3 <- data.frame(cbind(rep(v1bas,each=lg2),
                           rep(v2bas,lg1),interm))
names(infoM3) <- names(infoM1)
infoM <- data.frame(cbind(rbind(infoM1,infoM2,infoM3),
                  c(rep("PEPBVS",nrow(infoM1)),rep("BayesVarSel",nrow(infoM2)),
                            rep("BAS",nrow(infoM3)))))
names(infoM) <- NULL
resMAP <- data.frame(rbind(res1pepMAP,res2bayesvarselMAP,res3basMAP))
resMPM <- data.frame(rbind(res1pepMPM,res2bayesvarselMPM,res3basMPM))
colnames(resMAP) <- colnames(resMPM) <- predctrs
if ((algorithmic.choice=="MCMC")||(p>=20)){
 resMAP[,p+1] <- c(mcerrorpepMAP,mcerrorbayesvarselMAP,rep(NA,nrow(infoM3)))
 resMPM[,p+1] <- c(mcerrorpepMPM,mcerrorbayesvarselMPM,rep(NA,nrow(infoM3)))
 colnames(resMAP)[p+1] <- colnames(resMPM)[p+1] <- "MCerror"
}
resMAPa <- data.frame(cbind(infoM[,1:3],resMAP,infoM[,4]))
resMPMa <- data.frame(cbind(infoM[,1:3],resMPM,infoM[,4]))
nbcols <- ncol(resMAPa)
colnames(resMAPa)[c(1:3,nbcols)] <- colnames(resMPMa)[c(1:3,nbcols)] <- 
              c("prior", "model prior", "baseline/prior param", "R package")
if ((algorithmic.choice=="MCMC")||(p>=20)){
 resMAPa[,nbcols+1] <- resMPMa[,nbcols+1] <-
              c(rep("MC3",nrow(infoM1)),rep("Gibbs",nrow(infoM2)),
                            rep("MCMC",nrow(infoM3)))
 colnames(resMAPa)[nbcols+1] <- colnames(resMPMa)[nbcols+1] <- "algorithm"
}
res <- list(resMAPa,resMPMa)
names(res) <- c("MAPmodels", "MPMmodels")
return(res)
}