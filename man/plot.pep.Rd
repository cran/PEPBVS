% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/basic_complementary_functions.R
\name{plot.pep}
\alias{plot.pep}
\title{Plots for object of class pep}
\usage{
\method{plot}{pep}(x, ...)
}
\arguments{
\item{x}{An object of class pep (e.g. output of \code{full_enumeration_pep} 
or \code{mc3_pep}).}

\item{...}{Additional graphical parameters to be passed to plotting functions.}
}
\value{
No return value, used for generating a figure.
}
\description{
Generates four plots related to an object of class pep. In particular,
the first one is a plot of the residuals against fitted values under 
Bayesian model averaging. The second plots the cumulative posterior 
probability of the top models (those with cumulative posterior probability 
larger than 0.99). The third plot depicts the marginal likelihood 
(in log scale) of a model against its dimension while 
the fourth plot shows the posterior inclusion probabilities
of the explanatory variables (with those exceeding 0.5 marked in red).
}
\details{
Let k be the number of models with cumulative posterior probability up 
to 0.99. Then, the second plot depicts the cumulative posterior probability 
of the top (k+1) models.

In the special case of no explanatory variables, the fourth plot with the
posterior inclusion probabilities is not generated.
}
\examples{
data(UScrime_data)
y <- UScrime_data[,"y"]
X <- UScrime_data[,-15]
res <- full_enumeration_pep(X,y)
plot(res)

}
\seealso{
\code{\link{image.pep}}
}
