checkinputvartype <- function(x,y,intrinsic,reference.prior,
                              beta.binom,constant.term,burnin,itermc3){ 
# Takes the arguments of full_enumeration_pep & mc3_pep and checks
# whether they are of the correct type or not. If one of them is not, the procedure
# is interrupted and a message is output.
#
 if(!missing(x)){     # in case the input matrix is not completely missing (p=0)
  if(!((is.matrix(x))&&(is.numeric(x)))){
                      # check if it is a matrix of numeric
   message("Error: x should be a matrix or a data frame of numeric.\n") 
   stop("Please respecify and call the function again.")
  } 
 }                    # check if the response vector is a vector of numeric
 if (!((is.vector(y))&&(is.numeric(y)))){
   message("Error: y should be a numeric vector.\n")
   stop("Please respecify and call the function again.")
 }                    # check if the argument intrinsic is boolean
 if (!is.logical(intrinsic)){
   message("Error: the argument intrinsic should be boolean.\n")
   stop("Please respecify and call the function again.")
 }                    # check if the argument intrinsic is boolean
 if (!is.logical(reference.prior)){
   message("Error: the argument reference.prior should be boolean.\n")
   stop("Please respecify and call the function again.")
 }                    # check if the argument beta.binom is boolean
 if (!is.logical(beta.binom)){
   message("Error: the argument beta.binom should be boolean.\n")
   stop("Please respecify and call the function again.")
 }                    # check if the argument constant.term is boolean
 if (!is.logical(constant.term)){
   message("Error: the argument constant.term should be boolean.\n")
   stop("Please respecify and call the function again.")
 }
 if(!missing(burnin)) # for the MC3 algorithm
                      # check that the burnin period is a non-negative integer
                      # smaller than itermc3
  if (!((burnin%%1==0)&&(burnin>=0)&&(burnin<itermc3))){
   message("Error: the argument burnin should be a non-negative integer 
        smaller than the total number of iterations.\n")
   stop("Please respecify and call the function again.")
 }
 if(!missing(itermc3))# for the MC3 algorithm
                      # check that the number of total iterations is a 
                      # positive integer
  if (!((itermc3%%1==0)&&(itermc3>0))){
   message("Error: the argument itermc3 should be a positive integer.\n")
   stop("Please respecify and call the function again.")
 }                    
}                     # end function checkinputvartype

checkdimenscomp <- function(x,y){
# Takes an input matrix x of size n1xp and a response vector y
# of length n2 and checks if they are compatible in terms of 
# dimension, i.e. if n1=n2.
#
 n1 <- nrow(x)        # number of rows of x
 n2 <- length(y)      # sample size
 if (n1!=n2){
  message("Error: the number of rows of the input matrix 
       does not match the output vector length.\n")
  stop("Please respecify and call the function again.") 
 }
}                     # end function checkdimenscomp

compareparamvssamplesize <- function(x){
# Takes an input matrix x and checks if the total number of
# parameters (i.e. p+2- p input variables + intercept + error variance) 
# is smaller or equal to the sample size n.
#
 p <- ncol(x)         # number of explanatory variables
 n <- nrow(x)         # sample size
 if (p>n-2){          # if the sample size does not exceed by 2 at least
                      # the number of explanatory variables
  message("Error: the number of parameters (p+2) exceeds the sample size.
       This setup is not supported by the package.\n")
  stop("Please respecify and call the function again.") 
 }
}                     # end function compareparamvssamplesize

checkmissingness <- function(x,y){
# Takes an input matrix x and a response vector y and checks
# if they contain missing values (NAs). Currently, the presence of NAs
# is not supported.
#
 if (sum(is.na(y))>0){
  message("Error: the output vector should not contain missing values.\n")
  stop("Please respecify and call the function again.")
 }
 if (sum(is.na(x))>0){
  message("Error: the input matrix should not contain missing values.\n")
  stop("Please respecify and call the function again.")
 }
}                     # end function checkmissingness

checkmatrixrank <- function(x){
# Takes an input matrix x and checks if it is of full rank.
#
 if (Matrix::rankMatrix(x)[1]!=ncol(x)){
   message("Error: the input matrix is not full rank.\n")
   stop("Please respecify and call the function again.")
 }
}                     # end function checkmatrixrank

checkinputvartypeprint <- function(n.models, actual.PO, digits){
# Takes the arguments of print.pep and checks whether they are of the 
# correct type or not. If one of them is not, the printing call
# is interrupted and a message is output.
#
 if (!((n.models%%1==0)&&(n.models>0))){
                      # check that the number of models to be printed
                      # is a positive integer
   message("Error: the argument n.models should be a positive integer.")
   stop("Please respecify and call the function again.")
 }                    # check if the argument actual.PO is boolean
 if (!is.logical(actual.PO)){
   message("Error: the argument actual.PO should be boolean.\n")
   stop("Please respecify and call the function again.")
 }                    # check that the number of significant digits
                      # is a positive integer
 if (!((digits%%1==0)&&(digits>0))){
   message("Error: the argument digits should be a positive integer.")
   stop("Please respecify and call the function again.")
 }
}                     # end function checkinputvartypeprint

checkinputvartypepredict <- function(x,xnew,estimator,n.models,cumul.prob){
# Takes the arguments of predict.pep and checks whether they are of the 
# correct type or not. If one of them is not, the predict.pep call
# is interrupted and a message is output.
#
                      # check if the new data is a matrix of numeric (or NULL)
 if(!(((is.matrix(xnew))&&(is.numeric(xnew)))||is.null(xnew))){
   message("Error: xnew should be a matrix or a data frame of numeric.\n") 
   stop("Please respecify and call the function again.")
 }
 if (sum(is.na(xnew))>0){ # check that xnew does not contain NAs
  message("Error: the matrix xnew should not contain missing values.\n")
  stop("Please respecify and call the function again.")
 }
                      # check that the new and original data matrix are
                      # compatible in terms of dimension (i.e. they have the
                      # same number of input variables)
 if(!(is.null(xnew)||is.null(x))){
  p1 <- ncol(x)
  p2 <- ncol(xnew)
  if (p1!=p2){
   message("Error: the matrix xnew does not have the correct number 
       of explanatory variables.\n")
   stop("Please respecify and call the function again.")
  }
 }                    # check if the number of models 
                      # is a positive integer or NULL
 if ((!((n.models%%1==0)&&(n.models>0)))&&(!is.null(n.models))){
   message("Error: the argument n.models should be either a positive 
        integer or NULL.")
   stop("Please respecify and call the function again.")
 }                    # check if the argument estimator is one of 
                      # BMA, MAP or MPM
 if (!(is.character(estimator)&&estimator%in%c("BMA","MAP","MPM"))){
   message("Error: estimator should be one of BMA, MAP or MPM.\n")
   stop("Please respecify and call the function again.")
 }                    # check if the cumulative probability is a
                      # number between 0 and 1
 if (!((is.numeric(cumul.prob))&&(cumul.prob>=0)&&(cumul.prob<=1))){
   message("Error: cumul.prob should be a number between 0 and 1.\n")
   stop("Please respecify and call the function again.")
 }
}                     # end function checkinputvartypepredict

