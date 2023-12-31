% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_bioint.R
\name{bioint}
\alias{bioint}
\title{Computation of Bioclimatic Intensities}
\usage{
bioint(bb)
}
\arguments{
\item{bb}{Bioclimatic balance.}
}
\value{
data.frame with 10 variables. See details.
}
\description{
Computes bioclimatic intensities from bioclimatic balance.
}
\details{
The function yields 10 variables at monthly scale corresponding with the warm (w) and cold (c) variants of 5 bioclimatic intensities: PBI (Potential bioclimatic intensity), RBI (Real bioclimatic intensity), CBI (Conditioned bioclimatic intensity), FBI (Free bioclimatic intensity), and DBI (Dry bioclimatic intensity).
}
\examples{
wb <- watbal(t = rnorm(12, 18, 6), p = rnorm(12, 50, 30), lat = 35, CC = 400)
bb <- biobal(wb, 400)
bi <- bioint(bb)
}
