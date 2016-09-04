
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
    sqrt.cov <- eig.cov$vectors %*% diag(sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
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
    beta.init  <- eig.V$vectors[,1:d,drop=FALSE]

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }
    directions <- x.tilde %*% beta.init
    gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
    best.h.init     <- h[which.min(gcv.vals)]

    est.eqn <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d)
        directions <- x.tilde %*% beta.mat
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
        best.h     <- best.h.init # h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, deg = 3, ...)


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
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
        best.h     <- best.h.init # h[which.min(gcv.vals)]
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, deg = 3, ...)
        locfit.mod.deriv <- locfit.raw(x = directions, y = y, alpha = best.h, deriv = 1, deg = 3, ...)


        Ey.given.xbeta       <- fitted(locfit.mod)
        Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

        resid    <- drop(y - Ey.given.xbeta)
        psi      <- crossprod(x.tilde, resid * x.tilde) / nobs
        ## psi gradient with respect to just one column of beta
        psi.grad <- -crossprod(x.tilde, (Ey.given.xbeta.deriv * rowSums(x.tilde ^ 2) ) ) / nobs
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
                      gr      = est.eqn.grad,
                      method  = "L-BFGS",
                      control = list(maxit = maxit, reltol = 1e-8))

    #slver <-   tnewton(x0      = as.vector(beta.init), # beta.init[(d+1):nrow(beta.init),],
    #                     gr      = est.eqn.grad,
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


semi.phd.for.hier <- function(x, y, d = rep(2L, 3L), maxit = 10L, h = NULL,
                              strat.id = rep(1, nrow(x)), constraints, ...)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    nobs  <- nrow(x)
    nvars <- ncol(x)
    pp <- nvars / 3
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x.tilde)

    D <- sum(d)
    #dd <- floor(d/3)
    beta.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))
    cum.d <- c(0, cumsum(d))
    for (c in 1:length(constraints))
    {
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

            eig.c <- eigen(Proj.constr.list[[c]] %*% V.hat )
            eta.hat <- eig.c$vectors[, 1:d[c] ]
            beta.list[[c]] <- Re(eta.hat)# t(t(eta.hat) %*% sqrt.inv.cov)
        }
    }

    #eig.V <- eigen(V.hat)
    beta.init  <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]
    #eig.V <- eigen(V.hat)
    #beta.init  <- eig.V$vectors[,1:d,drop=FALSE]

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)

    #for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    ## est bandwidth

    best.h.vec <- numeric(num.strata)
    for (s in 1:num.strata)
    {
        strata.idx <- which(strat.id == unique.strata[s])
        dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
            beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
        dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

        gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
                                                 y = y[strata.idx],
                                                 alpha = hv, deg = 3, ...)[4])

        best.h.vec[s] <- h[which.min(gcv.vals)]
    }


    est.eqn.no.norm <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(as.vector(beta.vec), ncol = D)

        beta.list <- vector(mode = "list", length = length(constraints))
        for (c in 1:length(constraints))
        {
            if (d[c] > 0)
            {
                #AA <- constraints[[c]]
                #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.mat[, ((c-1) * dd + 1):(c * dd)]))
                beta.list[[c]] <- Proj.constr.list[[c]] %*% beta.mat[, (cum.d[c] + 1):cum.d[c + 1]]
            }
        }
        beta.mat <- do.call(cbind, beta.list)

        #directions <- x.tilde %*% beta.mat

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx], alpha = hv, deg = 3, ...)[4])
            best.h     <- h[which.min(gcv.vals)] # best.h.vec[s] #
            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx], alpha = best.h, deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- crossprod(x.tilde, resid * x.tilde) / (nobs)
        lhs
    }

    est.eqn <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(as.vector(beta.vec), ncol = D )

        beta.list <- vector(mode = "list", length = length(constraints))
        for (c in 1:length(constraints))
        {
            if (d[c] > 0)
            {
                #AA <- constraints[[c]]
                #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.mat[, ((c-1) * dd + 1):(c * dd)]))
                beta.list[[c]] <- Proj.constr.list[[c]] %*% beta.mat[, (cum.d[c] + 1):cum.d[c + 1]]
            }
        }
        beta.mat <- do.call(cbind, beta.list)

        #directions <- x.tilde %*% beta.mat

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx], alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)] #
            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx], alpha = best.h, deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }

    est.eqn.grad <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = D)

        beta.list <- vector(mode = "list", length = length(constraints))
        for (c in 1:length(constraints))
        {
            if (d[c] > 0)
            {
                #AA <- constraints[[c]]
                #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.mat[, ((c-1) * dd + 1):(c * dd)]))
                beta.list[[c]] <- Proj.constr.list[[c]] %*% beta.mat[, (cum.d[c] + 1):cum.d[c + 1]]# - adj.fact
            }
        }
        beta.mat <- do.call(cbind, beta.list)


        #directions <- x.tilde %*% beta.mat

        Ey.given.xbeta <- Ey.given.xbeta.deriv <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx], alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)] #
            locfit.mod <- locfit.raw(x = dir.cur,
                                     y = y[strata.idx], alpha = best.h, deg = 3, ...)
            locfit.mod.deriv <- locfit.raw(x = dir.cur,
                                           y = y[strata.idx], alpha = best.h, deriv = 1, deg = 3, ...)

            Ey.given.xbeta[strata.idx]       <- fitted(locfit.mod)
            Ey.given.xbeta.deriv[strata.idx] <- fitted(locfit.mod.deriv)
        }

        resid    <- drop(y - Ey.given.xbeta)
        psi      <- crossprod(x.tilde, resid * x.tilde) / nobs
        ## psi gradient with respect to just one column of beta
        psi.grad <- -crossprod(x.tilde, (Ey.given.xbeta.deriv * rowSums(x.tilde ^ 2) ) ) / nobs
        gradient <- 2 * t(psi %*% psi.grad)
        ## psi gradient is essentially just grad of one column repeated d times
        rep(drop(gradient), D)
    }

    #slver <- BBoptim(par = init,
    #               fn = est.eqn,
    #               #method = "SANN",
    #               control = list(maxit = maxit,
    #                              maxfeval = maxit * 25))

    slver <-   optim(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
                     fn      = est.eqn,
                     gr      = est.eqn.grad,
                     method  = "L-BFGS",
                     control = list(maxit = maxit, reltol = 1e-8))

    #slver <- nlminb(start     = as.vector(beta.init),
    #                objective = est.eqn,
    #                gradient  = est.eqn.grad,
    #                control   = list(trace = 0))

    #slver <-   ucminf(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
    #                  fn      = est.eqn,
    #                  gr      = est.eqn.grad
    #                  #method  = "L-BFGS",
    #                  #control = list(maxit = maxit, reltol = 1e-8)
    #                  )

    g <- as.vector(est.eqn.no.norm(slver$par))
    #W <- tcrossprod(g) + 0.1 * diag(length(g))
    W <- NULL

    beta <- beta.init
    for (i in 1:maxit)
    {

    }

    beta.semi <- rbind(beta.init[1:D, ], matrix(slver$par, ncol = D))
    beta.semi <- matrix(slver$par, ncol = D)


    beta.semi <- beta.semi # t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.semi, beta.init = beta, solver.obj = slver,
         cov = cov, sqrt.inv.cov = sqrt.inv.cov, W = W)
}



semi.phd.for.hier.W <- function(x, y, W, d = 5L, maxit = 10L, h = NULL,
                                init = NULL,
                                strat.id = rep(1, nrow(x)), constraints, ...)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    nobs  <- nrow(x)
    nvars <- ncol(x)
    pp <- nvars / 3
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x.tilde)

    dd <- floor(d/3)
    beta.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

        Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

        eig.c <- eigen(Proj.constr.list[[c]] %*% V.hat )
        eta.hat <- eig.c$vectors[,1:dd]
        beta.list[[c]] <- Re(eta.hat)# t(t(eta.hat) %*% sqrt.inv.cov)
    }

    #eig.V <- eigen(V.hat)
    beta.init  <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]
    #eig.V <- eigen(V.hat)
    #beta.init  <- eig.V$vectors[,1:d,drop=FALSE]

    if (!is.null(init)) beta.init <- init

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)

    #for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    ## est bandwidth

    best.h.vec <- numeric(num.strata)
    for (s in 1:num.strata)
    {
        strata.idx <- which(strat.id == unique.strata[s])
        dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
            beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
        dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

        gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
                                                 y = y[strata.idx],
                                                 alpha = hv, deg = 3, ...)[4])

        best.h.vec[s] <- h[which.min(gcv.vals)]
    }


    est.eqn <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(as.vector(beta.vec), ncol = d)

        beta.list <- vector(mode = "list", length = length(constraints))
        for (c in 1:length(constraints))
        {
            #AA <- constraints[[c]]
            #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.mat[, ((c-1) * dd + 1):(c * dd)]))
            beta.list[[c]] <- Proj.constr.list[[c]] %*% beta.mat[, ((c-1) * dd + 1):(c * dd)]
        }
        beta.mat <- do.call(cbind, beta.list)

        #directions <- x.tilde %*% beta.mat

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx], alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx], alpha = best.h, deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }

        resid <- drop(y - Ey.given.xbeta)
        #g     <- as.vector(crossprod(x.tilde, resid * x.tilde))

        g     <- as.vector(crossprod(x.tilde, resid * x.tilde))
        lhs   <- norm(crossprod(g, W %*% g), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }

    est.eqn.grad <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d)

        beta.list <- vector(mode = "list", length = length(constraints))
        for (c in 1:length(constraints))
        {
            #AA <- constraints[[c]]
            #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.mat[, ((c-1) * dd + 1):(c * dd)]))
            beta.list[[c]] <- Proj.constr.list[[c]] %*% beta.mat[, ((c-1) * dd + 1):(c * dd)]# - adj.fact
        }
        beta.mat <- do.call(cbind, beta.list)


        #directions <- x.tilde %*% beta.mat

        Ey.given.xbeta <- Ey.given.xbeta.deriv <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx], alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = dir.cur,
                                     y = y[strata.idx], alpha = best.h, ...)
            locfit.mod.deriv <- locfit.raw(x = dir.cur,
                                           y = y[strata.idx], alpha = best.h, deriv = 1, deg = 3, ...)

            Ey.given.xbeta[strata.idx]       <- fitted(locfit.mod)
            Ey.given.xbeta.deriv[strata.idx] <- fitted(locfit.mod.deriv)
        }

        resid    <- drop(y - Ey.given.xbeta)
        psi      <- as.vector(crossprod(x.tilde, resid * x.tilde) / nobs)
        ## psi gradient with respect to just one column of beta
        psi.grad <- drop(-crossprod(x.tilde, (Ey.given.xbeta.deriv * rowSums(x.tilde ^ 2) ) ) / nobs)
        print(str(psi))
        print(str(W))
        print(str(psi.grad))
        gradient <- 2 * psi %*% W %*% psi.grad
        ## psi gradient is essentially just grad of one column repeated d times
        rep(drop(gradient), d)
    }

    #slver <- BBoptim(par = init,
    #               fn = est.eqn,
    #               #method = "SANN",
    #               control = list(maxit = maxit,
    #                              maxfeval = maxit * 25))

    slver <-   optim(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
                     fn      = est.eqn,
                     gr      = est.eqn.grad,
                     method  = "L-BFGS",
                     control = list(maxit = maxit, reltol = 1e-8))

    beta <- beta.init
    for (i in 1:maxit)
    {

    }

    beta.semi <- rbind(beta.init[1:d, ], matrix(slver$par, ncol = d))
    beta.semi <- matrix(slver$par, ncol = d)


    beta.semi <- beta.semi # t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.semi, beta.init = beta, solver.obj = slver,
         cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}


semi.phd.hier <- function(x.list, y, d = rep(1L, 3L), maxit = 10L, h = NULL, ...)
{
    p <- ncol(x.list[[1]])
    x <- as.matrix(bdiag(x.list))

    nobs  <- nrow(x)
    nvars <- ncol(x)
    pp <- nvars / 3

    D <- sum(d)

    cum.d <- c(0, cumsum(d))

    cov <- cov(as.matrix(x))
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                                cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                        t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                                cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                        t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                                 cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )

    #for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    strat.id <- unlist(lapply(1:length(x.list), function(id) rep(id, nrow(x.list[[id]]))))

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x.tilde)


    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

            eig.c <- eigen(Proj.constr.list[[c]] %*% V.hat )
            eta.hat <- eig.c$vectors[,1:d[c], drop=FALSE]
            beta.list[[c]] <- Re(eta.hat)

            beta.init.list[[c]] <- Proj.constr.list[[c]] %*%
                matrix(runif(prod(dim(eta.hat)),
                             min = min(beta.list[[c]]),
                             max = max(beta.list[[c]])),
                       ncol = NCOL(eta.hat))
        }
    }

    Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]

    #eig.V <- eigen(V.hat)
    beta.init  <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]
    beta.rand.init <- do.call(cbind, beta.init.list)

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)




    #######    model fitting to determine bandwidth     ########

    best.h.vec <- numeric(num.strata)

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    beta.init.cov <- t(t(beta.init) %*% sqrt.inv.cov)
    for (s in 1:num.strata)
    {
        strata.idx <- which(strat.id == unique.strata[s])
        dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
            beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
        #print(apply(dir.cur, 2, sd))
        dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

        gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
                                                 y = y[strata.idx],
                                                 kern = "gauss",
                                                 alpha = hv, deg = 3, ...)[4])
        best.h.vec[s]     <- h[which.min(gcv.vals)]
    }

    est.eqn <- function(beta.vec)
    {
        beta.mat   <- matrix(beta.vec, ncol = D)
            #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         kern = "gauss",
            #                                         alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                     kern = "gauss",
                                     alpha = best.h, deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }

    est.eqn.grad <- function(beta.vec)
    {
        beta.mat   <- matrix(beta.vec, ncol = D)
            #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)

        Ey.given.xbeta <- Ey.given.xbeta.deriv <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
                beta.mat[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         alpha = hv, kern = "gauss", deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = dir.cur,
                                     y = y[strata.idx],
                                     kern = "gauss",
                                     alpha = best.h, ...)
            locfit.mod.deriv <- locfit.raw(x = dir.cur,
                                           y = y[strata.idx],
                                           alpha = best.h,
                                           kern = "gauss",
                                           deriv = 1, deg = 3, ...)

            Ey.given.xbeta[strata.idx]       <- fitted(locfit.mod)
            Ey.given.xbeta.deriv[strata.idx] <- fitted(locfit.mod.deriv)
        }

        resid    <- drop(y - Ey.given.xbeta)
        psi      <- crossprod(x.tilde, resid * x.tilde) / nobs
        ## psi gradient with respect to just one column of beta
        psi.grad <- -crossprod(x.tilde, (Ey.given.xbeta.deriv * rowSums(x.tilde ^ 2) ) ) / nobs
        gradient <- 2 * t(psi %*% psi.grad)
        ## psi gradient is essentially just grad of one column repeated d times
        #rep(drop(gradient), d)
        unlist(lapply(Proj.constr.list, function(xx) drop(xx %*% drop(gradient) ) ))
    }

    #slver <- BBoptim(par = as.vector(beta.init),
    #               fn = est.eqn,
    #               gr = est.eqn.grad,
    #               control = list(maxit = maxit,
    #                              maxfeval = maxit * 4))

    slver <-   optim(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
                     fn      = est.eqn,
                     gr      = est.eqn.grad,
                     method  = "L-BFGS",
                     control = list(maxit = maxit, factr = 1e-10))

    #slver <- nlminb(start     = as.vector(beta.init),
    #                objective = est.eqn,
    #                gradient  = est.eqn.grad,
    #                control   = list(trace = 0))

    #slver <-   optimx(par     = as.vector(beta.init), # beta.init[(d+1):nrow(beta.init),],
    #                  fn      = est.eqn,
    #                  gr      = est.eqn.grad,
    #                  method  = "spg",
    #                  control = list(maxit = maxit, reltol = 1e-8,
    #                                 kkt = FALSE, starttests = FALSE))
    #slver <- list(obj = slver, par = as.numeric(slver[1,1:length(as.vector(beta.init))]))

    #print(str(slver$par))
    beta <- beta.init
    for (i in 1:maxit)
    {

    }

    #beta.semi <- rbind(beta.init[1:d, ], matrix(slver$par, ncol = d))
    beta.semi <- matrix(slver$par, ncol = D)
        #t(t(matrix(slver$par, ncol = d * length(constraints))) %*% sqrt.inv.cov)


    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.semi, beta.init = beta, solver.obj = slver,
         beta.rand.init = t(t(beta.rand.init) %*% sqrt.inv.cov),
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



# from http://www4.stat.ncsu.edu/~li/software/GroupDR.R
# angle between two spaces
angles<-function(B1, B2)
{
    if(!is.matrix(B1)) B1 <- as.matrix(B1)
    if(!is.matrix(B2)) B2 <- as.matrix(B2)

    if(ncol(B1) >= ncol(B2)) {
        B <- B1; B.hat <- B2
    } else {
        B <- B2; B.hat <- B1
    }

    P1 <- B %*% solve(t(B) %*% B) %*% t(B)
    if(ncol(B.hat) == 1) {
        nume  <- as.vector(t(B.hat) %*% P1 %*% B.hat)
        deno  <- as.vector(t(B.hat) %*% B.hat)
        ratio <- nume / deno
    } else {
        BtB   <- t(B.hat) %*% B.hat
        ei    <- eigen(BtB)
        BtB2  <- ei$vectors %*% diag(1/sqrt(ei$values)) %*% t(ei$vectors)
        M     <- BtB2 %*% t(B.hat) %*% P1 %*% B.hat %*% BtB2
        ratio <- abs(eigen(M)$values[nrow(M)])
    }
    ans <- acos(sqrt(ratio)) / pi * 180
    if(ans > 90) ans <- 180 - ans
    return(ans)
}



# normalize a vector
norm2 <- function(v)
{
    sumv2<-sum(v^2)
    if(sumv2 == 0) sumv2<-1
    v/sqrt(sumv2)
}

# Gram-Schmidt orthonormalization
orthnorm<-function(X)
{
    X<-as.matrix(X)
    n<-nrow(X)
    p<-ncol(X)

    W<-NULL
    if(p > 1) {
        W<-cbind(W, X[,1])
        for(k in 2:p) {
            gw<-rep(0, n)
            for(i in 1:(k-1)) {
                gki<-as.vector((t(W[,i]) %*% X[,k])/(t(W[,i]) %*% W[,i]))
                gw<-gw + gki * W[,i]
            }
            W<-cbind(W, X[,k] - gw)
        }
    } else {
        W<-cbind(W, X[,1])
    }

    W<-apply(W, 2, norm)
    W
}

