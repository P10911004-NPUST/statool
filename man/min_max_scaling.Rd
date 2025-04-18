% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/feature_scaling.R
\name{min_max_scaling}
\alias{min_max_scaling}
\title{Minimum-Maximum (Min-Max) scaling}
\usage{
min_max_scaling(x, range = c(0, 1))
}
\arguments{
\item{x}{A numeric vector.}

\item{range}{A numeric vector with length of 2. The desired range for the data normalization.
The default is c(0, 1).}
}
\value{
A numeric vector, i.e., the normalized data.
}
\description{
Feature scaling method normalizes the input data to a desired range (generally from 0 to 1).
The transformation was done by:
\verb{[(x - min(x)) / (max(x) - min(x))] * (range_ceiling - range_floor) + range_floor}
}
\examples{
set.seed(1)
x <- stats::rnorm(10)
min_max_scaling(x)
#> [1] 0.086047962 0.419296552 0.000000000 1.000000000 0.479300618
#> [6] 0.006236443 0.544264487 0.647475101 0.580609856 0.218124222

min_max_scaling(x, range = c(-2, 2))
#> [1] -1.65580815 -0.32281379 -2.00000000  2.00000000 -0.08279753
#> [6] -1.97505423  0.17705795  0.58990040  0.32243942 -1.12750311
}
\author{
Joon-Keat Lai
}
