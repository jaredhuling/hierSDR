
##################################################################################
## requires the following variables to be set:
##
## nobs      - number of observations per group
## ncats     - number of condition categories (not used yet)
## nvars     - number of variables
## x.type    - either "simple" (for iid standard normal x)
##             or "complex" for AR cont. x + 5 binary x + 5 cont. depend on other x
## beta.type -  either "some.zero" for half of betas to be zero
##              or "some.small" for half of beta to be very small
## model.num - either 1 for the first model type or 2 for second model
##


library(MASS)

ndata <- 2 ^ ncats - 1

if (x.type == "simple")
{
    x.list      <- replicate(ndata, list(matrix(rnorm(nobs      * nvars), ncol = nvars)))
    x.list.test <- replicate(ndata, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
} else if (x.type == "complex")
{
    rho <- 0.5
    cov.mat <- rho ^ abs(outer(1:(nvars - 10), 1:(nvars - 10), "-"))


    x.list <- x.list.test <- vector(mode = "list", length = ndata)
    for (g in 1:ndata)
    {
        x.tmp  <- mvrnorm(n = nobs, mu = rep(0, nvars - 10), Sigma = cov.mat)
        x.tmp2 <- apply(x.tmp[,1:5], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
        x.tmp3 <- apply(x.tmp[,6:10], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
        x.list[[g]] <- data.matrix(cbind(x.tmp, x.tmp2, x.tmp3))

        x.tmp  <- mvrnorm(n = nobs.test, mu = rep(0, nvars - 10), Sigma = cov.mat)
        x.tmp2 <- apply(x.tmp[,1:5], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
        x.tmp3 <- apply(x.tmp[,6:10], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
        x.list.test[[g]] <- data.matrix(cbind(x.tmp, x.tmp2, x.tmp3))
    }
}


x <- bdiag(x.list)
x.test <- as.matrix(bdiag(x.list.test))


cor.directions <- function(a, b, x)
{
    cov <- cov(x)
    R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
    R.sq
}



if (beta.type == "some.zero")
{
    beta.a <- matrix(c(numeric(nvars / 2),
                       runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
    beta.b <- matrix(c(numeric(nvars / 2),
                       runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
    eta.ab <- matrix(c(numeric(nvars / 2),
                       runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
} else
{
    beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                       rnorm(nvars / 2, sd = 0.25)), ncol = 1)
    beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                       rnorm(nvars / 2, sd = 0.25)), ncol = 1)
    eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                       rnorm(nvars / 2, sd = 0.25)), ncol = 1)
}


mult.a <- diag(runif(nvars, min = 0.5, max = 1.5))
mult.b <- diag(runif(nvars, min = 0.5, max = 1.5))

eta.zero.models <- c(2, 4)

if (model.num %in% eta.zero.models)
{
    beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)

} else
{
    beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b)
}


if (model.num == 1)
{
    y.true.a <- apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) +
        0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) + # ^ 2 +
        0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))

    y.true <- c(y.true.a, y.true.b, y.true.ab)
    y <- y.true + rnorm(nobs, sd = sd.sim)


    y.true.a <- apply( exp(x.list.test[[1]] %*% beta.a), 1, sum)
    y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list.test[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) +
        0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) + # ^ 2 +
        0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

    y.true.test <- c(y.true.a, y.true.b, y.true.ab)
    y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)




} else if (model.num == 2)
{

    y.true.a <- 0.5 * apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 3)
    y.true.ab <- x.list[[3]] %*% beta.ab[,1] / (0.5 + (x.list[[3]] %*% beta.ab[,2] + 1.5) ^ 2)

    y.true <- c(y.true.a, y.true.b, y.true.ab)
    y <- y.true + rnorm(nobs, sd = sd.sim)


    y.true.a <- 0.5 * apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 3)
    y.true.ab <- x.list.test[[3]] %*% beta.ab[,1] / (0.5 + (x.list.test[[3]] %*% beta.ab[,2] + 1.5) ^ 2)

    y.true.test <- c(y.true.a, y.true.b, y.true.ab)
    y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

} else if (model.num == 3)
{
    y.true.a <- 0.25 * apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) -
        0.5 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) + # ^ 2 +
        0.25 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))

    y.true <- c(y.true.a, y.true.b, y.true.ab)
    y <- y.true + rnorm(nobs, sd = sd.sim)


    y.true.a <- 0.25 * apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list.test[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) -
        0.5 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) + # ^ 2 +
        0.25 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

    y.true.test <- c(y.true.a, y.true.b, y.true.ab)
    y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)
} else if (model.num == 4)
{
    y.true.a <- 0.25 * apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) -
        0.5 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2)

    y.true <- c(y.true.a, y.true.b, y.true.ab)
    y <- y.true + rnorm(nobs, sd = sd.sim)


    y.true.a <- 0.25 * apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)# ^ 2
    y.true.ab <- (apply( (x.list.test[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) -
        0.5 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2)

    y.true.test <- c(y.true.a, y.true.b, y.true.ab)
    y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)
} else if (model.num == 5)
{
    y.true.a <- 5 * apply(cos(x.list[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)
    y.true.ab <- 5 * apply(sin(x.list[[3]] %*% beta.ab[,2]) , 1, sum) -
        0.5 * ((x.list[[3]] %*% beta.ab[,1]) ^ 2) *  (apply(sin(x.list[[3]] %*% beta.ab[,3]), 1, sum))

    y.true <- c(y.true.a, y.true.b, y.true.ab)
    y <- y.true + rnorm(nobs, sd = sd.sim)


    y.true.a <- 5 * apply(cos(x.list.test[[1]] %*% beta.a) , 1, sum)
    y.true.b <- + 0.5 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)
    y.true.ab <- 5 * apply(sin(x.list.test[[3]] %*% beta.ab[,2]) , 1, sum) -
        0.5 * ((x.list.test[[3]] %*% beta.ab[,1]) ^ 2) *  (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

    y.true.test <- c(y.true.a, y.true.b, y.true.ab)
    y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)
}

snr <- var(y.true.test) / var(y.test - y.true.test)



