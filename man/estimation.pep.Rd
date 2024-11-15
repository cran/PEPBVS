% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_complementary_functions.R
\name{estimation.pep}
\alias{estimation.pep}
\title{Model averaged estimates}
\usage{
estimation.pep(
  object,
  ssize = 10000,
  estimator = "BMA",
  n.models = NULL,
  cumul.prob = 0.99
)
}
\arguments{
\item{object}{An object of class pep (e.g., output of \code{pep.lm}).}

\item{ssize}{Positive integer, the number of values to be simulated from
the (joint) posterior distribution of the beta coefficients.
Default value=10000.}

\item{estimator}{A character, the type of estimation. One of 
``BMA'' (Bayesian model averaging, default), 
``MAP'' (maximum a posteriori model) or ``MPM'' (median probability model).
Default value=\code{"BMA"}.}

\item{n.models}{Positive integer, the number of (top) models where
the average is based on or \code{NULL}. Relevant for \code{estimator="BMA"}.
Default value=\code{NULL}.}

\item{cumul.prob}{Numeric between zero and one, cumulative probability of
top models to be used for computing the average. Relevant for \code{estimator="BMA"}. 
Default value=0.99.}
}
\value{
\code{estimation.pep} returns a matrix (of dimension 
\code{ssize} \eqn{\times \, (p+1)}) --- 
where the rows correspond
to the simulations and the columns to the beta coefficients
(including the intercept) --- containing the 
simulated data.
}
\description{
Simulates values from the (joint) posterior distribution of the 
beta coefficients under Bayesian model averaging.
}
\details{
For the computations, Equation 10 of Garcia--Donato and Forte (2018) 
is used. That (simplified) formula arises when changing the prior on the
model parameters to the reference prior. This change of prior is
justified in Garcia--Donato and Forte (2018). The resulting formula is a mixture
distribution and the simulation is implemented as follows: firstly the 
model (component) based on its posterior probability is chosen and 
subsequently the values of the beta coefficients included in the chosen model are
drawn from the corresponding multivariate Student distribution, while the
values of the beta coefficents outside the chosen model are set to zero.

Let \eqn{k} be the number of models with cumulative posterior probability up 
to the given value of \code{cumul.prob}. Then, for Bayesian model averaging 
the summation is based on the top \eqn{(k+1)} models if they exist, otherwise
on the top \eqn{k} models.

When both \code{n.models} and \code{cumul.prob} are provided --- once 
specifying the number of models for the given cumulative probability as 
described above --- the minimum between the two numbers is used for estimation.
}
\examples{
data(UScrime_data)
res <- pep.lm(y~.,data=UScrime_data)
set.seed(123)
estM1 <- estimation.pep(res,ssize=2000)
estM2 <- estimation.pep(res,ssize=2000,estimator="MPM")
}
\references{
Garcia--Donato, G. and Forte, A. (2018) Bayesian Testing, 
Variable Selection and Model Averaging in Linear Models using R with 
BayesVarSel. The R Journal, 10(1): 155–174. 
\doi{10.32614/RJ-2018-021}
}
