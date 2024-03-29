% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{projnorm}
\alias{projnorm}
\title{Norm of difference of projections}
\usage{
projnorm(B1, B2)
}
\arguments{
\item{B1}{first matrix}

\item{B2}{second matrix}
}
\value{
scalar value of the projection difference norm between \code{B1} and \code{B2}
}
\description{
Measures distance between two subspaces
}
\examples{
b1 <- matrix(rnorm(10 * 2), ncol = 2)
b2 <- matrix(rnorm(10 * 2), ncol = 2)
projnorm(b1, b2)

## angle here should be smalls
b1 <- matrix(rnorm(10 * 2), ncol = 2)
b2 <- b1 + matrix(rnorm(10 * 2, sd = 0.2), ncol = 2)
projnorm(b1, b2)
}
