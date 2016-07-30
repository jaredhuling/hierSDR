
Kepanechnikov  <- function(u) 0.75 * (1 - (u) ^ 2) * (abs(u) < 1)
Kepanechnikov2 <- function(u) 0.75 * (1 - (u) ^ 2)


semiDR <- function(x, y, d = 5L, maxit = 10L)
{

}

createDiffEpa <- function(x, h = 1)
{
    diffs <- outer(x, x, FUN = "-") / h
    diffs[abs(diffs) >= 1] <- 0
    diffs <- as(diffs, "sparseMatrix")
    diffs
}

nwsmooth <- function(x, y, h = 1)
{
    diffmat   <- createDiffEpa(x, h = h)
    diffmat@x <- Kepanechnikov2(diffmat@x) / h
    predicted.values <- colSums(y * diffmat) / colSums(diffmat)
    list(fitted = predicted.values)
}

nwsmoothcov <- function(x, y, h = 1)
{
    nobs <- NROW(x)
    diffmat   <- createDiffEpa(y, h = h)
    diffmat@x <- Kepanechnikov2(diffmat@x) / h
    txpy <- predicted.values <- vector(mode = "list", length = nobs)

    csums <- colSums(diffmat)
    for (i in 1:nobs) txpy[[i]] <- tcrossprod(x[i,])
    for (i in 1:nobs)
    {
        sum.cov <- txpy[[1]] * diffmat[i,1]
        for (j in 2:nobs) sum.cov <- sum.cov + txpy[[j]] * diffmat[i,j]

        predicted.values[[i]] <- sum.cov / csums[i]
    }
    list(fitted = predicted.values)
}


