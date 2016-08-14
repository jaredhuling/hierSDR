
Kepanechnikov  <- function(u) 0.75 * (1 - (u) ^ 2) * (abs(u) < 1)
Kepanechnikov2 <- function(u) 0.75 * (1 - (u) ^ 2)


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
    nobs        <- NROW(x)
    #diffmat    <- createDiffEpa(x, h = h)
    dist.idx    <- fields.rdist.near(x, x, 1/h, max.points = nobs * 100)
    dist.idx$ra <- dist.idx$ra / h
    diffmat     <- sparseMatrix(dist.idx$ind[,1], dist.idx$ind[,2],
                                x = dist.idx$ra,  dims = dist.idx$da)
    diffmat@x   <- Kepanechnikov2(diffmat@x) / h

    RSS <- sum(((Diagonal(nobs) - diffmat) %*% y)^2)/nobs
    predicted.values <- Matrix::colSums(y * diffmat) / Matrix::colSums(diffmat)
    trS <- sum(Matrix::diag(diffmat))
    #gcv <- (1 / nobs) * sum(( (y - predicted.values) / (1 - trS / nobs) ) ^ 2)
    gcv <- RSS / ((1 - trS / nobs) ) ^ 2
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

directional.regression <- function(x, y, d = 5L)
{

}

phd <- function(x, y, d = 5L)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x)
    eig.V <- eigen(V.hat)
    eta.hat <- eig.V$vectors[,1:d]
    beta.hat <- t(t(eta.hat) %*% sqrt.inv.cov)
    list(beta.hat = beta.hat, eta.hat = eta.hat, M = V.hat, cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}

semi.phd <- function(x, y, d = 5L, maxit = 10L, h = NULL)
{
    phd.mod <- phd(x = x, y = y, d = d)
    cov          <- phd.mod$cov
    sqrt.inv.cov <- phd.mod$sqrt.inv.cov
    beta         <- phd.mod$beta.hat
    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 100))
    }

    est.eqn <- function(beta.vec)
    {
        beta.mat   <- matrix(beta.vec, ncol = d)
        directions <- x %*% beta
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv))
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- y - Ey.given.xbeta
        lhs   <- norm(crossprod(x, resid * x), type = "F") ^ 2
        lhs
    }

    slver <- optim(par = as.vector(beta),
                   fn = est.eqn,
                   method = "BFGS",
                   control = list(maxit = maxit))

    beta.semi <- matrix(slver$par, ncol = d)
    #for (i in 1:maxit)
    #{
        #beta.prev  <- beta
        #directions <- x %*% beta
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv))
        #best.h     <- h[which.min(gcv.vals)]
        #locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)

        #Ey.given.xbeta <- fitted(locfit.mod)
    #}
    list(beta = beta.semi, beta.init = beta, solver.obj = slver)
}


sir <- function(x, y, h = 10L, d = 5L, slice.ind = NULL)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    if (is.null(slice.ind))
    {
        quantiles <- quantile(y, probs = seq(0, 1, length.out = h+1))
        quantiles[1] <- quantiles[1] - 1e-5
        quantiles[length(quantiles)] <- quantiles[length(quantiles)] + 1e-5
        y.cut <- cut(y, quantiles)
        p.hat <- sapply(levels(y.cut), function(lv) mean(y.cut == lv))

        x.h <- t(sapply(levels(y.cut), function(lv) colMeans(x.tilde[y.cut == lv,])))
    } else
    {
        lvls <- sort(unique(slice.ind))
        p.hat <- sapply(lvls, function(lv) mean(slice.ind == lv))
        x.h <- t(sapply(lvls, function(lv) colMeans(x.tilde[slice.ind == lv,])))
    }

    V.hat <- crossprod(x.h, p.hat * x.h)
    eig.V <- eigen(V.hat)
    eta.hat <- eig.V$vectors[,1:d]
    beta.hat <- t(t(eta.hat) %*% sqrt.inv.cov)
    list(beta.hat = beta.hat, eta.hat = eta.hat)
}

