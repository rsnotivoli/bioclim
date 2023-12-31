% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/f_thermind.R
\name{thermind}
\alias{thermind}
\title{Function to Compensated Thermal Index}
\usage{
thermind(t)
}
\arguments{
\item{t}{Monthly average temperature data (12 nueric values).}
}
\value{
Compensated Thermal Index
}
\description{
Computes Compensated Thermal Index from monthly temperature.
}
\examples{
thermind(rnorm(12, 18, 6))
}
