
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
semi.phd2 <- function(x, y, d = 5L, maxit = 10L, h = NULL)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x)
    eig.V <- eigen(V.hat)
    beta  <- eig.V$vectors[,1:d]
    beta.init <- beta
    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    est.eqn <- function(beta.vec)
    {
        beta.mat   <- rbind(diag(d), matrix(beta.vec, ncol = d))
        directions <- x %*% beta
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv)[4])
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- y - Ey.given.xbeta
        lhs   <- norm(crossprod(x, resid * x), type = "F") ^ 2
        lhs
    }

    #  beta[(d+1):nrow(beta),]
    objective <- numeric(maxit)
    for (i in 1:maxit)
    {
        directions <- x.tilde %*% beta
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv)[4])
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)


        resid <- y - fitted(locfit.mod)

        P <- beta %*% solve(crossprod(beta), t(beta))
        V.hat <- crossprod(x.tilde, resid * x.tilde) / nrow(x)
        eig.V <- eigen(V.hat)
        beta  <- eig.V$vectors[,1:d]

        objective[i] <- norm(V.hat - P %*% V.hat %*% P, type = "F") ^ 2
    }

    beta.hat  <- t(t(beta) %*% sqrt.inv.cov)
    beta.init <- t(t(beta.init) %*% sqrt.inv.cov)

    list(beta = beta.hat, beta.init = beta.init, objective = objective,
         M = V.hat, cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}



semi.phd <- function(x, y, d = 5L, maxit = 10L, h = NULL, ...)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    nobs  <- nrow(x)
    nvars <- ncol(x)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x.tilde)
    eig.V <- eigen(V.hat)
    beta.init  <- eig.V$vectors[,1:d]

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    est.eqn <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d)
        directions <- x.tilde %*% beta.mat
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, ...)[4])
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, ...)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }

    est.eqn.grad <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d)
        directions <- x.tilde %*% beta.mat
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, ...)[4])
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, ...)
        locfit.mod.deriv <- locfit.raw(x = directions, y = y, alpha = best.h, deriv = 1, ...)


        Ey.given.xbeta       <- fitted(locfit.mod)
        Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

        resid    <- drop(y - Ey.given.xbeta)
        psi      <- crossprod(x.tilde, resid * x.tilde) / nobs
        ## psi gradient with respect to just one column of beta
        psi.grad <- crossprod(x.tilde, (Ey.given.xbeta * Ey.given.xbeta.deriv * rowSums(x.tilde ^ 2) ) ) / nobs
        gradient <- 2 * t(psi %*% psi.grad)
        ## psi gradient is essentially just grad of one column repeated d times
        rep(drop(gradient), d)
    }

    est.eqn2 <- function(beta.mat)
    {
        beta.mat2  <- rbind(beta.init[1:d,], beta.mat)
        directions <- x.tilde %*% beta.mat2
        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv)[4])
        best.h     <- h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- y - Ey.given.xbeta
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2
        lhs
    }

    #  beta[(d+1):nrow(beta),]

    init <- rep(1, length(as.vector(beta.init[(d+1):nrow(beta.init),])) )


    #slver <- BBoptim(par = init,
    #               fn = est.eqn,
    #               #method = "SANN",
    #               control = list(maxit = maxit,
    #                              maxfeval = maxit * 25))


    #beta <- as.vector(beta.init)
    #for (i in 1:maxit)
    #{
    #    prev <- beta
    #    gradient <- grad(est.eqn, x = prev, method = "simple")
    #    jacobi   <- jacobian(est.eqn, x = prev, method = "simple")
    #    beta     <- beta - (gradient / jacobi)
    #
    #    print(i)
    #    if (max(abs(beta - prev) < 1e-5))
    #    {
    #        break
    #    }
    #}

    #slver <- list(par = beta, value = est.eqn(beta))

     slver <-   optim(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
                      fn      = est.eqn,
                      method  = "CG",
                      control = list(maxit = maxit))

    #slver <-   sbplx(x0      = as.vector(beta.init), # beta.init[(d+1):nrow(beta.init),],
    #                     fn      = est.eqn,
    #                  #lower = rep(-100, length(beta.init)),
    #                  #upper = rep( 100, length(beta.init)),
    #                     control = list(maxeval = maxit * 2))

    #slver <- DEoptim(fn = est.eqn,
    #                 lower = rep(-1e5, length(init)),
    #                 upper = rep(1e5, length(init)))

    #beta.semi <- rbind(diag(d), matrix(slver$par, ncol = d))


    #beta <- beta.init[(d+1):nrow(beta.init),]
    beta <- beta.init
    for (i in 1:maxit)
    {

    }

    beta.semi <- rbind(beta.init[1:d, ], matrix(slver$par, ncol = d))
    beta.semi <- matrix(slver$par, ncol = d)


    beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- t(t(beta.init) %*% sqrt.inv.cov)

    #for (i in 1:maxit)
    #{
        #beta.prev  <- beta
        #directions <- x %*% beta
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv))
        #best.h     <- h[which.min(gcv.vals)]
        #locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)

        #Ey.given.xbeta <- fitted(locfit.mod)
    #}
    list(beta = beta.semi, beta.init = beta, solver.obj = slver,
         cov = cov, sqrt.inv.cov = sqrt.inv.cov)
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

