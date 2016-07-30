
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
    nobs      <- NROW(x)
    diffmat   <- createDiffEpa(x, h = h)
    diffmat@x <- Kepanechnikov2(diffmat@x) / h
    predicted.values <- colSums(y * diffmat) / colSums(diffmat)
    trS <- sum(diag(diffmat))
    gcv <- (1 / nobs) * sum(( (y - predicted.values) / (1 - trS / nobs) ) ^ 2)
    list(fitted = predicted.values, gcv = gcv)
}

nwsmoothcov <- function(x, y, h = 1)
{
    nobs <- NROW(x)
    diffmat   <- createDiffEpa(y, h = h)
    diffmat@x <- Kepanechnikov2(diffmat@x) / h
    txpy <- predicted.values <- vector(mode = "list", length = nobs)
    trS <- sum(diag(diffmat))

    csums <- colSums(diffmat)
    for (i in 1:nobs) txpy[[i]] <- tcrossprod(x[i,])
    for (i in 1:nobs)
    {
        sum.cov <- txpy[[1]] * diffmat[1,i]
        for (j in 2:nobs) sum.cov <- sum.cov + txpy[[j]] * diffmat[j,i]

        predicted.values[[i]] <- sum.cov / csums[i]
    }

    normdiff <- 0
    for (j in 1:nobs) normdiff <- normdiff + norm(txpy[[j]] - predicted.values[[j]], type = "F")

    gcv <- (1 / nobs) * ((normdiff / (1 - trS / nobs)) ^ 2)

    list(fitted = predicted.values, gcv = gcv)
}


