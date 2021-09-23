

#' Simulate data with hierarchical subspaces
#'
#' @description Simulates data with hierarchical subspaces. Data are generated with two factors that induce heterogeneity
#'
#' @param nobs positive integer for the sample size per subpopulation
#' @param nvars positive integer for the dimension
#' @param x.type variable type for covariates, either \code{"continuous"} (where the covariates are multivariate normal with a variance-matrix
#' with AR-1 form with parameter \code{rho}) or \code{"some_categorical"} (where half covariates are continuous and
#' the other half are binary with dependencies on the continuous covariates)
#' @param sd.y standard deviation of responsee
#' @param rho correlation parameter for AR-1 covariance structure for continuous covariates
#' @param model model number used, either "1", "2", or "3", each corresponds to a different outcome model setting
#' @importFrom stats rnorm var
#' @importFrom MASS mvrnorm
#' @return A list with the following elements
#' \itemize{
#' \item x a matrix of covariates with number of rows equal to the total sample size and columns equal to the number of variables
#' \item z a matrix with number of rows equal to the total sample size and columns as dummy variables indicating presence of a stratifying factor
#' \item y a vector of all responses
#' \item beta a list of the true sufficient dimension reduction matrices, one for each subpopulation
#' \item z.combinations all possible combinations of the stratifying factors \code{z}
#' \item snr scalar the observed signal-to-noise ratio for the response
#' \item d.correct the true dimensions of the dimension reduction spaces
#' }
#' @export
#' @examples
#'
#' library(hierSDR)
#'
#' set.seed(123)
#' dat <- simulate_data(nobs = 100, nvars = 6,
#'                      x.type = "some_categorical",
#'                      sd.y = 1, model = 2)
#'
#' x <- dat$x ## covariates
#' z <- dat$z ## factor indicators
#' y <- dat$y ## response
#'
#' dat$beta ## true coefficients that generate the subspaces
#'
#' dat$snr ## signal-to-noise ratio
#'
#' str(x)
#' str(z)
#'
#' dat$z.combinations ## what combinations of z represent different subpops
#'
#' ## correct structural dimensions:
#' dat$d.correct
#'
#'
simulate_data <- function(nobs,
                          nvars,
                          x.type = c("continuous", "some_categorical"),
                          sd.y = 1,
                          rho = 0.5,
                          model = c("1", "2", "3"))
{


    # stopifnot(
    #     "need nfactors > 1" = nfactors > 1,
    #     "need nfactors < 6" = nfactors < 6,
    #     "need n > 9" = n > 9,
    #     "need p > 0" = p > 0
    # )

    nfactors <- 2

    stopifnot(
        nfactors > 1,
        nfactors < 6,
        nobs > 9,
        nvars > 1,
        nvars < 101,
        rho < 1,
        rho > -1,
        sd.y > 0
    )



    ndata <- 2 ^ nfactors

    model.num <- as.numeric(model)



    if (x.type == "continuous")
    {
        cov.mat     <- rho ^ abs(outer(1:(nvars), 1:(nvars), "-"))
        x.list      <- replicate(ndata, list( mvrnorm(n = nobs, mu = rep(0, nvars), Sigma = cov.mat) ))
    } else if (x.type == "some_categorical")
    {

        if (nvars <= 10)
        {
            subtr <- floor(nvars / 2)
        } else
        {
            subtr <- 10
        }

        cov.mat <- rho ^ abs(outer(1:(nvars - subtr), 1:(nvars - subtr), "-"))


        x.list <- vector(mode = "list", length = ndata)
        for (g in 1:ndata)
        {
            x.tmp  <- mvrnorm(n = nobs, mu = rep(0, nvars - subtr), Sigma = cov.mat)

            if (subtr == 10)
            {
                x.tmp2 <- apply(x.tmp[,1:5], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
                x.tmp3 <- apply(x.tmp[,6:10], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
            } else
            {
                x.tmp2 <- apply(x.tmp[,1:floor(subtr/2),drop=FALSE], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
                x.tmp3 <- apply(x.tmp[,(floor(subtr/2) + 1):subtr,drop=FALSE], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
            }

            x.list[[g]] <- data.matrix(cbind(x.tmp, x.tmp2, x.tmp3))
        }
    }


    # put all subpops into one tall matrix
    x.tall <- do.call(rbind, x.list)


    # factor dummy variables
    z <- cbind(c(rep(0,nrow(x.list[[1]])), rep(1,nrow(x.list[[2]])),  rep(0, nrow(x.list[[3]]) ), rep(1, nrow(x.list[[4]]) )  ),
               c(rep(0,nrow(x.list[[1]])), rep(0,nrow(x.list[[2]])),  rep(1, nrow(x.list[[3]]) ), rep(1, nrow(x.list[[4]]) )  ))

    eta.b2.models   <- c(3)
    eta.zero.models <- c(2, 3)



    ## correct structural dimensions
    if (model.num == 1)
    {
        d.correct <- c(1, 0, 0, 0)
    } else if (model.num == 2)
    {
        d.correct <- c(1, 0, 0, 1)
    } else if (model.num == 3)
    {
        d.correct <- c(1, 1, 0, 0)
    }


    #
    #  beta_none - common for all. in this sim will only have 1 column
    #
    #  beta_a    - will have dimension 1 for model _____
    #
    #  beta_b    - will always be of dimension zero (ie we won't have it)
    #
    #  beta_ab   - will have dimension 1 for model _____
    #


    ## randomly generate true coefficients
    beta.none <- matrix(c(runif(floor(nvars / 2), min = 0, max = 0.25) * (2 * rbinom(floor(nvars / 2), 1, 0.5) - 1),
                          runif(ceiling(nvars / 2), min = 0, max = 0.5) * (2 * rbinom(ceiling(nvars / 2), 1, 0.5) - 1)), ncol = 1)

    eta.a     <- matrix(c(runif(floor(nvars / 2), min = 0, max = 0.25) * (2 * rbinom(floor(nvars / 2), 1, 0.5) - 1),
                          runif(ceiling(nvars / 2), min = 0, max = 0.5) * (2 * rbinom(ceiling(nvars / 2), 1, 0.5) - 1)), ncol = 1)

    eta.ab    <- matrix(c(runif(floor(nvars / 2), min = 0, max = 0.25) * (2 * rbinom(floor(nvars / 2), 1, 0.5) - 1),
                          runif(ceiling(nvars / 2), min = 0, max = 0.5) * (2 * rbinom(ceiling(nvars / 2), 1, 0.5) - 1)), ncol = 1)


    beta.none[1:ncol(beta.none),] <- 0.5*diag(ncol(beta.none))
    eta.a[1:ncol(eta.a),]         <- 0.5*diag(ncol(eta.a))
    eta.ab[1:ncol(eta.ab),]       <- 0.5*diag(ncol(eta.ab))

    beta.none <- grassmannify(beta.none)$beta
    eta.a     <- grassmannify(eta.a)$beta
    eta.ab    <- grassmannify(eta.ab)$beta



    if (model.num == 1)
    {
        beta.a <- beta.b <- beta.ab <- beta.none
    } else if (model.num == 2)
    {
        beta.a  <- beta.b <- beta.none
        beta.ab <- cbind(beta.none, eta.ab)
    } else if (model.num == 3)
    {
        beta.a  <- cbind(beta.none, eta.a)
        beta.b  <- beta.none
        beta.ab <- cbind(beta.none, eta.a)
    }

    beta.true.list <- list(none = beta.none, a = beta.a, b = beta.b, ab = beta.ab)

    if (model.num == 1)
    {
        y.true.none <- 2 * ((x.list[[1]] %*% beta.none[,1]) ^ 2)
        y.true.a <- 0.5 * (x.list[[2]] %*% beta.a[,1]) ^ 3
        y.true.b <- 2 * (x.list[[3]] %*% beta.b[,1])
        y.true.ab <- (x.list[[4]] %*% beta.ab[,1]) / (0.5 + (x.list[[4]] %*% beta.ab[,1] + 1.5) ^ 2)

        y.true <- c(y.true.none, y.true.a, y.true.b, y.true.ab)
        y      <- y.true + rnorm(length(y.true), sd = sd.y)

    } else if (model.num == 2)
    {

        y.true.none <- 2 * ((x.list[[1]] %*% beta.none[,1]) ^ 2)
        y.true.a    <- 1 / (0.1 + 0.5 * drop(x.list[[2]] %*% beta.a[,1]) ^ 2) - 0.5 * drop(x.list[[2]] %*% beta.a[,1]) ^ 2
        y.true.b    <- 2 * ((x.list[[3]] %*% beta.b[,1]) ^ 2)# ^ 2
        y.true.ab   <- 2 * ((x.list[[4]] %*% beta.ab[,1]) ^ 2) +
            (x.list[[4]] %*% beta.ab[,1]) * (x.list[[4]] %*% beta.ab[,2])

        y.true <- c(y.true.none, y.true.a, y.true.b, y.true.ab)
        y      <- y.true + rnorm(length(y.true), sd = sd.y)

    } else if (model.num == 3)
    {
        y.true.none <- 1 * abs(x.list[[1]] %*% beta.none[,1])^2
        y.true.a    <- 0.5 * (x.list[[2]] %*% beta.a[,1]) * ((x.list[[2]] %*% beta.a[,2])) ^ 2
        y.true.b    <- 1 * apply(exp( (x.list[[3]] %*% beta.b) ) , 1, sum)
        y.true.ab   <- 5 * drop(exp(-(x.list[[4]] %*% beta.ab[,1])^2 )) * drop(x.list[[4]] %*% beta.ab[,2])

        y.true <- c(y.true.none, y.true.a, y.true.b, y.true.ab)
        y      <- y.true + rnorm(length(y.true), sd = sd.y)

    }

    snr <- sqrt(var(y.true)) / sqrt(var(y - y.true))

    z.combinations <- subpop.struct(2L, incl.none = TRUE)

    colnames(z.combinations) <- colnames(z) <- c("a", "b")


    list(x = x.tall,
         z = z,
         y = drop(y),
         beta = beta.true.list,
         snr = snr,
         z.combinations = z.combinations,
         d.correct = d.correct)
}
