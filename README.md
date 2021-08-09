
<!-- README.md is generated from README.Rmd. Please edit that file -->

# hierSDR

<!-- badges: start -->
<!-- badges: end -->

The `hierSDR` package implements the methodology of Huling and Yu
(2021+) which is a semiparametric sufficient dimension reduction
approach that handles heterogeneity in the relationship between
covariates and outcome where the heterogeneity is based on all
combinations of a set of binary factors. A motivating example is
modeling health risks based on covariates where the relationship between
covariates and outcome changes depending on what combination of chronic
conditions a patient has. In this example the chronic conditions are the
factors that drive heterogeneity. Each subpopulation is the set of
patients with a different combination of chronic conditions.

## Installation

You can install the released version of hierSDR from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("hierSDR")
```

And the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("jaredhuling/hierSDR")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(hierSDR)
```

Simulate data where response mean is generated from a multiple index
model where each combination of a set of binary factors *Z* has a
potentially different set of parameters in the multiple index model.
Here, *Z* drives heterogeneity in the relationship between covariates
*X* and a response *Y*.

``` r
set.seed(123)
dat <- simulate_data(nobs = 500, nvars = 6,
                     x.type = "some_categorical",
                     sd.y = 1, model = 2)

x <- dat$x ## covariates
z <- dat$z ## factor indicators, these drive heterogeneity in the data
y <- dat$y ## response

## true coefficients that generate the subspaces
dat$beta 
#> $none
#>            [,1]
#> [1,] 1.00000000
#> [2,] 0.06586526
#> [3,] 0.17240683
#> [4,] 0.10866231
#> [5,] 0.84253756
#> [6,] 0.77444774
#> 
#> $a
#>            [,1]
#> [1,] 1.00000000
#> [2,] 0.06586526
#> [3,] 0.17240683
#> [4,] 0.10866231
#> [5,] 0.84253756
#> [6,] 0.77444774
#> 
#> $b
#>            [,1]
#> [1,] 1.00000000
#> [2,] 0.06586526
#> [3,] 0.17240683
#> [4,] 0.10866231
#> [5,] 0.84253756
#> [6,] 0.77444774
#> 
#> $ab
#>            [,1]       [,2]
#> [1,] 1.00000000  1.0000000
#> [2,] 0.06586526  0.2317807
#> [3,] 0.17240683 -0.2069714
#> [4,] 0.10866231  0.7276493
#> [5,] 0.84253756  0.8026606
#> [6,] 0.77444774 -0.5029732

## what combinations of z represent different subpopulations
dat$z.combinations 
#>   a b
#> 1 0 0
#> 2 1 0
#> 3 0 1
#> 4 1 1

## correct structural dimensions:
## (the dimensions of the subpopulation-specific
##  central mean subspaces)
dat$d.correct
#> [1] 1 0 0 1
```

Now estimate the parameters in the multiple index model where each
combination of the factors in *Z* has its own set of parameters, but the
parameters are constrained so that the central mean subspaces (which
convey the regression relationship information) are hierarchical
according to the structure of the subpopulations defined by *Z*.

Here we are specifying the correct dimensions for each of the parameter
matrices for each subpopulation. However, the correct dimensions can be
estimated with a validated information criterion (VIC).

``` r
hiermod <- hier.sphd(x, y, z, dat$z.combinations, d = dat$d.correct,
                     verbose = FALSE)


## validated inf criterion for choosing dimensions 
## (the smaller the better).
## multiple models should be fit with differing dimensions and the
## dimensions with the smallest VIC should be used
hiermod$vic
#> [1] 507.1956

## true and estimated parameters
cbind(Estimated = hiermod$beta[[4]], NA, Truth = dat$beta[[4]])
#>           [,1]      [,2] [,3]       [,4]       [,5]
#> [1,] 1.0000000 1.0000000   NA 1.00000000  1.0000000
#> [2,] 0.3552258 0.3562017   NA 0.06586526  0.2317807
#> [3,] 0.2889527 0.2616284   NA 0.17240683 -0.2069714
#> [4,] 0.1832808 0.2413655   NA 0.10866231  0.7276493
#> [5,] 0.8526007 0.8487122   NA 0.84253756  0.8026606
#> [6,] 0.8043703 0.6349214   NA 0.77444774 -0.5029732

## angles between estimated and true subspaces for each subpopulation:
## (smaller means the subspaces are better estimated)
mapply(function(x,y) angle(x,y), hiermod$beta, dat$beta)
#>      0,0      1,0      0,1      1,1 
#> 11.34241 11.34241 11.34241 11.56792

## projection difference norm between estimated and true subspaces for each population:
## (smaller means the subspaces are better estimated)
mapply(function(x,y) projnorm(x,y), hiermod$beta, dat$beta)
#>       0,0       1,0       0,1       1,1 
#> 0.2781362 0.2781362 0.2781362 0.3871345
```

## References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-huling21semiparametric" class="csl-entry">

Huling, Jared D., and Menggang Yu. 2021+. “Sufficient Dimension
Reduction for Populations with Structured Heterogeneity.” *Biometrics*,
2021+.

</div>

</div>
