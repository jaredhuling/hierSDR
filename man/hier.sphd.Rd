% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hier_sphd.R
\name{hier.sphd}
\alias{hier.sphd}
\title{Main hierarchical sufficient dimension reduction fitting function}
\usage{
hier.sphd(
  x,
  y,
  z,
  z.combinations,
  d,
  weights = rep(1L, NROW(y)),
  maxit = 250L,
  tol = 1e-09,
  h = NULL,
  opt.method = c("lbfgs2", "lbfgs.x", "bfgs.x", "bfgs", "lbfgs", "spg", "ucminf", "CG",
    "nlm", "nlminb", "newuoa"),
  init.method = c("random", "phd"),
  vic = TRUE,
  grassmann = TRUE,
  nn = NULL,
  nn.try = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95),
  n.random = 100L,
  optimize.nn = FALSE,
  separate.nn = FALSE,
  constrain.none.subpop = TRUE,
  verbose = TRUE,
  degree = 2,
  pooled = FALSE,
  maxk = 5000,
  ...
)
}
\arguments{
\item{x}{an n x p matrix of covariates, where each row is an observation and each column is a predictor}

\item{y}{vector of responses of length n}

\item{z}{an n x C matrix of binary indicators, where each column is a binary variable indicating the presence
of a binary variable which acts as a stratifying variable. Each combination of all columns of \code{z} pertains
to a different subpopulation. WARNING: do not use too many binary variables in \code{z} or else it will quickly
result in subpopulations with no observations}

\item{z.combinations}{a matrix of dimensions 2^C x C with each row indicating a different combination of the possible
values in \code{z}. Each combination represents a subpopulation. This is necessary because we need to specify a
different structural dimension for each subpopulation, so we need to know the ordering of the subpopulations so we
can assign each one a structural dimension}

\item{d}{an integer vector of length 2^C of structural dimensions. Specified in the same order as the rows in
\code{z.combinations}}

\item{weights}{vector of observation weights}

\item{maxit}{maximum number of iterations for optimization routines}

\item{tol}{convergence tolerance for optimization routines. Defaults to \code{1e-6}}

\item{h}{bandwidth parameter. By default, a reasonable choice is selected automatically}

\item{opt.method}{optimization method to use. Available choices are
\code{c("lbfgs2", "lbfgs.x", "bfgs.x", "bfgs", "lbfgs", "spg", "ucminf", "CG", "nlm", "nlminb", "newuoa")}}

\item{init.method}{method for parameter initialization. Either \code{"random"} for random initialization or \code{"phd"}
for a principle Hessian directions initialization approach}

\item{vic}{logical value of whether or not to compute the VIC criterion for dimension determination}

\item{grassmann}{logical value of whether or not to enforce parameters to be on the Grassmann manifold}

\item{nn}{nearest neighbor parameter for \code{\link[locfit]{locfit.raw}}}

\item{nn.try}{vector of nearest neighbor parameters for \code{\link[locfit]{locfit.raw}} to try in random initialization}

\item{n.random}{integer number of random initializations for parameters to try}

\item{optimize.nn}{should \code{nn} be optimized? Not recommended}

\item{separate.nn}{should each subpopulation have its own \code{nn}? If \code{TRUE}, optimization takes
much longer. It is rarely better, so recommended to set to \code{FALSE}}

\item{constrain.none.subpop}{should the "none" subpopulation be constrained to be contained in every other subpopulation's
dimension reduction subspace? Recommended to set to \code{TRUE}}

\item{verbose}{should results be printed along the way?}

\item{degree}{degree of kernel to use}

\item{pooled}{should the estimator be a pooled estimator?}

\item{maxk}{maxk parameter for \code{\link[locfit]{locfit.raw}}. Set to a large number if an out of vertex space error occurs.}

\item{...}{extra arguments passed to \code{\link[locfit]{locfit.raw}}}
}
\value{
A list with the following elements
\itemize{
\item beta a list of estimated sufficient dimension reduction matrices, one for each subpopulation
\item beta.init a list of the initial sufficient dimension reduction matrices, one for each subpopulation -- do not use, just for the sake of comparisons
\item directions a list of estimated sufficient dimension reduction directions (i.e. the reduced dimension predictors/variables), one for each subpopulation.
These have number of rows equal to the sample size for the subpopulation and number of columns equal to the specified dimensions of the reduced dimension spaces.
\item y.list a list of vectors of responses for each subpopulation
\item z.combinations the \code{z.combinations} specified as an input
\item cov list of variance covariance matrices for the covariates for each subpopulation
\item sqrt.inv.cov list of inverse square roots of the variance covariance matrices for the covariates for each subpopulation. These are used for scaling
\item solver.obj object returned by the solver/optimization function
\item value value of the objective function at the solution
\item value.init value of the objective function at the initial beta (\code{beta.init}) used
\item vic.est.eqn the average (unpenalized) VIC value  across the r different input values. This assesses model fit
\item vic.eqns the individual (unpenalized) VIC values across the r input values. Not used.
\item vic the penalized VIC value. This is used for dimension selection, with dimensions chosen by the set of dimensions
that minimize this penalized vic value that trades off model complexity and model fit
}
}
\description{
fits hierarchically nested sufficient dimension reduction models
}
\examples{

library(hierSDR)

set.seed(123)
dat <- simulate_data(nobs = 200, nvars = 6,
                     x.type = "some_categorical",
                     sd.y = 1, model = 2)

x <- dat$x ## covariates
z <- dat$z ## factor indicators
y <- dat$y ## response

dat$beta ## true coefficients that generate the subspaces

dat$z.combinations ## what combinations of z represent different subpops

## correct structural dimensions:
dat$d.correct

## fit hier SPHD model:

\donttest{
hiermod <- hier.sphd(x, y, z, dat$z.combinations, d = dat$d.correct,
                     verbose = FALSE, maxit = 250, maxk = 8200)

## validated inf criterion for choosing dimensions (the smaller the better)
hiermod$vic


cbind(hiermod$beta[[4]], NA, dat$beta[[4]])

## angles between estimated and true subspaces for each population:
mapply(function(x,y) angle(x,y), hiermod$beta, dat$beta)

## projection difference norm between estimated and true subspaces for each population:
mapply(function(x,y) projnorm(x,y), hiermod$beta, dat$beta)
}


}
