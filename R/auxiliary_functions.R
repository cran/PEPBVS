integer.base.b <- function(x, b=2){
# Takes a vector of integers x and computes and returns
# the representation in base b (default to 2) for each integer.
# 
 xi <- as.integer(x)
 if(any(is.na(xi) | ((x-xi)!=0))){      # check whether given vector is an 
   message("Error: x should be a vector of integer") # integer vector
   stop("Please respecify and call the function again.")
 }
 N <- length(x)                        # vector length
 xMax <- max(x)                        # maximal element
                                       # number of digits needed 
 ndigits <- floor(logb(xMax, base=2))+1 # for the representation
 Base.b <- array(NA, dim=c(N, ndigits))
 for(i in 1:ndigits){                  # compute the representation
  Base.b[, ndigits-i+1] <- (x %% b)    # through successive divisions with
  x <- (x %/% b)                       # the base
 }
 if(N ==1) 
   Base.b[1, ] else 
   Base.b                              # return result
}                                      # end function integer.base.b


constant.marginal.likelihood <- function(y, d0=0){
# Takes a response vector y and the hyperparameter d0 (default value 0)
# and computes C1 from the paper in log scale (last equation on p. 1082).
#
 n <- length(y)                        # sample size
 # compute quantities needed for the result
 k0 <- 1                               # k0
 a0 <- (n+d0-k0)/2                     # a0
 RSS0 <- var(y)*(n-1)                  # residual sum of squares (RSS)
 logdetX0TX0 <- log(n)            # log|t(X0)*X0|
                                       # compute logC1
 logC1 <- (d0/2-1)*log(2)+0.5*(k0-n)*log(pi) -
          0.5*logdetX0TX0+lgamma(a0) -a0 *log(RSS0)
 return(logC1)
}                                      # end function constant.marginal.likelihood  

