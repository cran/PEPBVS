% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/UScrime_data.R
\docType{data}
\name{UScrime_data}
\alias{UScrime_data}
\title{US Crime Data}
\format{
This data frame contains the following columns:
\describe{
\item{\code{M}}{percentage of males aged 14--24.}
\item{\code{Ed}}{mean years of schooling.}
\item{\code{Po1}}{police expenditure in 1960.}
\item{\code{Po2}}{police expenditure in 1959.}
\item{\code{LF}}{labour force participation rate.}
\item{\code{M.F}}{number of males per 1000 females.}
\item{\code{Pop}}{state population.}
\item{\code{NW}}{number of non-whites per 1000 people.}
\item{\code{U1}}{unemployment rate of urban males 14--24.}
\item{\code{U2}}{unemployment rate of urban males 35--39.}
\item{\code{GDP}}{gross domestic product per head.}
\item{\code{Ineq}}{income inequality.}
\item{\code{Prob}}{probability of imprisonment.}
\item{\code{Time}}{average time served in state prisons.}
\item{\code{y}}{rate of crimes in a particular category per head of population.}
}
}
\source{
Data from the R package MASS
}
\description{
The dataset has been borrowed from the MASS R package and describes 
the effect of punishment regimes on crime rates. One explanatory
variable (indicator variable for a Southern state) was removed since it was binary.
}
