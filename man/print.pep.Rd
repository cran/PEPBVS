% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_complementary_functions.R
\name{print.pep}
\alias{print.pep}
\title{Printing object of class pep}
\usage{
\method{print}{pep}(
  x,
  n.models = 5,
  actual.PO = FALSE,
  digits = max(3L, getOption("digits") - 3L),
  ...
)
}
\arguments{
\item{x}{An object of class pep (e.g. output of \code{full_enumeration_pep} 
or \code{mc3_pep}).}

\item{n.models}{Positive integer, the number of top models for which information 
is provided. Default value=5.}

\item{actual.PO}{Boolean, relevant for the MC3 algorithm. If \code{TRUE}
then apart from the estimated posterior odds, the actual posterior
odds of the top models (i.e. ratios based on the marginal likelihood 
times prior probability) are also printed - which could be used as a 
convergence indicator of the algorithm. Default value=\code{FALSE}.}

\item{digits}{Positive integer, the number of digits for printing numbers. 
Default value=\code{max(3L, getOption("digits") - 3L)}.}

\item{...}{Additional parameters to be passed to 
\code{print.default}.}
}
\value{
No return value, used for printing the results on the R console.
}
\description{
For each of the top models (shown in columns), the following information is
printed: the model representation using variable inclusion indicators, 
its marginal likelihood (in log scale), the R2, the model dimension, the Bayes 
factor, posterior odds and posterior probability. An additional 
column with the posterior inclusion probabilities of the explanatory variables is 
also printed.
}
\details{
The number of models for which information is provided, is computed as the minimum
between the number asked by the user and the number of models present in
the object \code{x}.
}
\examples{
data(UScrime_data)
y <- UScrime_data[,"y"]
X <- UScrime_data[,-15]
res <- full_enumeration_pep(X,y)
print(res)

}
