
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

    y.scaled <- scale(y, scale = FALSE)
    miny <- min(y.scaled)

    V.hat <- crossprod(x.tilde, drop(y.scaled + miny * sign(miny) + 1e-5) * x.tilde) / nrow(x)
    eig.V <- eigen(V.hat)
    eta.hat <- eig.V$vectors[,1:d]
    beta.hat <- t(t(eta.hat) %*% sqrt.inv.cov)
    list(beta.hat = beta.hat, eta.hat = eta.hat, M = V.hat, cov = cov, sqrt.inv.cov = sqrt.inv.cov, eigenvalues = eig.V$values)
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



semi.phd <- function(x, y, d = 5L, maxit = 10L, h = NULL, vic = FALSE, B = NULL, ...)
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
        sd <- sd(directions)

        best.h <- sd * (0.75 * nrow(directions)) ^ (-1/(ncol(directions)+4) )
        locfit.mod <- locfit.raw(x = directions, y = y, alpha = c(0.75, best.h), deg = 2, ...)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }


    est.eqn.vic <- function(beta.vec)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d + 1)
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
                      #gr      = est.eqn.grad,
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

    #beta.semi <- rbind(beta.init[1:d, ], matrix(slver$par, ncol = d))
    beta.semi <- matrix(slver$par, ncol = d)


    if (vic)
    {
        beta.u <- beta.semi[1,,drop=FALSE]

        if (TRUE)
        {
            v.1 <- matrix(rep(1, nrow(beta.semi) - 1), ncol=1)
            vk.1 <- matrix(0, ncol=ncol(beta.semi)+1, nrow=nrow(beta.semi))
            vk.1[-1,-ncol(vk.1)] <- vk.1[-1,-ncol(vk.1)] - drop(v.1 %*% beta.u)
            vk.1[-1,ncol(vk.1)] <- v.1
            vk.1[1,] <- c(rep(0, ncol(vk.1) - 1), 1)
            vk.1.vec <- as.vector(vk.1)


            v.2 <- matrix(rep(0, nrow(beta.semi) - 1), ncol=1)
            vk.2 <- matrix(0, ncol=ncol(beta.semi)+1, nrow=nrow(beta.semi))
            vk.2[-1,-ncol(vk.2)] <- vk.2[-1,-ncol(vk.2)] - drop(v.2 %*% beta.u)
            vk.2[-1,ncol(vk.2)] <- v.2
            vk.2[1,] <- c(rep(0, ncol(vk.2) - 1), 1)
            vk.2.vec <- as.vector(vk.2)


            v.3 <- matrix(rep(-1, nrow(beta.semi) - 1), ncol=1)
            vk.3 <- matrix(0, ncol=ncol(beta.semi)+1, nrow=nrow(beta.semi))
            vk.3[-1,-ncol(vk.3)] <- vk.2[-1,-ncol(v.3)] - drop(v.3 %*% beta.u)
            vk.3[-1,ncol(vk.3)] <- v.3
            vk.3[1,] <- c(rep(0, ncol(vk.3) - 1), 1)
            vk.3.vec <- as.vector(vk.3)


            eqn.val.1 <- est.eqn.vic(vk.1.vec) * nobs
            eqn.val.2 <- est.eqn.vic(vk.2.vec) * nobs
            eqn.val.3 <- est.eqn.vic(vk.3.vec) * nobs

            vic <- (eqn.val.1 + eqn.val.2 + eqn.val.3) / 3 + nvars * d * log(nobs)
        } else
        {
            vic <- slver$value * nobs + nvars * d * log(nobs)
        }
    } else
    {
        vic <- NULL
    }



    beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- t(t(beta.init) %*% sqrt.inv.cov)

    rsq.mean <- NULL
    if (!is.null(B))
    {
        dir.orig <- x %*% beta.semi

        cov.u <- cov(dir.orig)
        eig.cov.u <- eigen(cov.u)
        print(dim(eig.cov.u$vectors))
        print(str(eig.cov.u))
        if (ncol(beta.semi) == 1)
        {
            sqrt.inv.cov.u <- 1 / sqrt(cov.u)
        } else
        {
            sqrt.inv.cov.u <- eig.cov.u$vectors %*% diag(1 / sqrt(eig.cov.u$values)) %*% t(eig.cov.u$vectors)
        }
        #boot.beta <- vector(mode = "list", length = B)
        rsq <- numeric(B)
        for (b in 1:B)
        {
            s.idx <- sample.int(nobs, size = nobs, replace = FALSE)
            x.b <- x[s.idx,]
            y.b <- y[s.idx]
            cov.b <- cov(x)
            eig.cov.b <- eigen(cov.b)
            sqrt.inv.cov.b <- eig.cov.b$vectors %*% diag(1 / sqrt(eig.cov.b$values)) %*% t(eig.cov.b$vectors)
            x.tilde.b <- scale(x.b, scale = FALSE) %*% sqrt.inv.cov.b

            est.eqn.b <- function(beta.vec)
            {
                beta.mat   <- matrix(beta.vec, ncol = d)
                directions <- x.tilde.b %*% beta.mat
                #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
                best.h     <- best.h.init # h[which.min(gcv.vals)]
                locfit.mod <- locfit.raw(x = directions, y = y.b, alpha = best.h, deg = 3, ...)


                Ey.given.xbeta <- fitted(locfit.mod)

                resid <- drop(y.b - Ey.given.xbeta)
                lhs   <- norm(crossprod(x.tilde.b, resid * x.tilde.b), type = "F") ^ 2 / (nobs ^ 2)
                lhs
            }

            slver.b <-   optim(par     = slver$par, # beta.init[(d+1):nrow(beta.init),],
                               fn      = est.eqn.b,
                               #gr      = est.eqn.grad,
                               method  = "L-BFGS",
                               control = list(maxit = 1L, reltol = 1e-8))

            beta.semi.b <- t(t(matrix(slver.b$par, ncol = d)) %*% sqrt.inv.cov)
            dir.boot <- x %*% beta.semi.b

            #cov.v <- cov(dir.boot)
            #eig.cov.v <- eigen(cov.v)
            #sqrt.inv.cov.v <- eig.cov.v$vectors %*% diag(1 / sqrt(eig.cov.v$values)) %*% t(eig.cov.v$vectors)
            if (ncol(beta.semi) == 1)
            {
                emat <- sqrt.inv.cov.u * drop(cov(dir.orig, dir.boot)) * (drop(cov(dir.boot, dir.orig))/drop(cov(dir.boot))) * sqrt.inv.cov.u
            } else
            {
                emat <- sqrt.inv.cov.u %*% cov(dir.orig, dir.boot) %*% solve(cov(dir.boot), cov(dir.boot, dir.orig)) %*% sqrt.inv.cov.u
            }
            rsq[b] <- mean(eigen(emat)$values)
        }
        rsq.mean <- mean(rsq)
    }






    directions <- x %*% beta.semi
    gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 2, ...)[4])
    best.h     <- h[which.min(gcv.vals)]
    locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, deg = 2, ...)

    #for (i in 1:maxit)
    #{
        #beta.prev  <- beta
        #directions <- x %*% beta
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv))
        #best.h     <- h[which.min(gcv.vals)]
        #locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h)

        #Ey.given.xbeta <- fitted(locfit.mod)
    #}
    list(beta         = beta.semi,
         beta.init    = beta,
         solver.obj   = slver,
         cov          = cov,
         sqrt.inv.cov = sqrt.inv.cov,
         final.gcv    = min(gcv.vals),
         final.model  = locfit.mod,
         all.gcvs     = gcv.vals,
         rsq          = rsq.mean,
         vic          = vic)
}

semi.phd.dim.select <- function(x, y, max.d = 5L, maxit = 10L, h = NULL, ...)
{
    best.gcv <- 1e99
    gcvs <- numeric(max.d)
    models <- vector(mode = "list", length = max.d)
    for (d in 1:max.d)
    {
        print(d)
        s.phd.cur <- semi.phd(x, y, d = d, maxit = maxit, h = h, ...)
        if (s.phd.cur$final.gcv < best.gcv)
        {
            best.gcv   <- s.phd.cur$final.gcv
            best.model <- s.phd.cur
            best.d     <- d
        }
        gcvs[d] <- s.phd.cur$final.gcv
        models[[d]] <- s.phd.cur
    }
    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcvs
    ret$models <- models
    ret
}


semi.phd.dim.select <- function(x, y, max.d = 5L, maxit = 10L, h = NULL, ...)
{
    best.gcv <- 1e99
    gcvs <- bics <- numeric(max.d)
    models <- vector(mode = "list", length = max.d)

    for (d in 1:max.d)
    {
        print(d)
        s.phd.cur <- semi.phd(x, y, d = d, maxit = maxit, h = h, ...)
        if (s.phd.cur$final.gcv < best.gcv)
        {
            best.gcv   <- s.phd.cur$final.gcv
            best.model <- s.phd.cur
            best.d     <- d
        }
        bics[d] <- s.phd.cur$solver.obj$value + 2 * log(nrow(x)) * d
        gcvs[d] <- s.phd.cur$final.gcv
        models[[d]] <- s.phd.cur
    }
    ret   <- best.model
    ret$d <- which.min(bics)
    ret$gcvs <- gcvs
    ret$models <- models
    ret$bic <- bics
    ret
}


semi.phd.dim.select.boot <- function(x, y, max.d = 5L, maxit = 10L, h = NULL, B = 50L, ...)
{
    best.gcv <- 1e99
    gcvs <- rsq <- numeric(max.d)
    models <- vector(mode = "list", length = max.d)

    for (d in 1:max.d)
    {
        print(d)
        s.phd.cur <- semi.phd(x, y, d = d, maxit = maxit, h = h, B = B, ...)
        if (s.phd.cur$final.gcv < best.gcv)
        {
            best.gcv   <- s.phd.cur$final.gcv
            best.model <- s.phd.cur
            best.d     <- d
        }
        rsq[d] <- s.phd.cur$rsq
        gcvs[d] <- s.phd.cur$final.gcv
        models[[d]] <- s.phd.cur
    }
    ret   <- best.model
    ret$d <- which.max(rsq)
    ret$gcvs <- gcvs
    ret$models <- models
    ret$rsq <- rsq
    ret
}

semi.phd.dim.select.k <- function(x, y, max.d = 5L, h = NULL, k = 10L, maxit = 100, ...)
{
    best.gcv <- 1e99
    gcvs <- array(NA, dim = c(k, max.d))
    models <- vector(mode = "list", length = max.d)
    nobs <- nrow(x)

    foldid <- sample(rep(seq(k), length = nobs))

    for (i in 1:k)
    {
        which = which(foldid == i)
        x.i     <- x[-which,]
        #phd.cur <- phd(x.i, y[-which], d = max.d)
        print(i)
        for (d in 1:max.d)
        {
            if (i == nobs) print(d)

            s.phd.cur <- semi.phd(x.i, y[-which], d = d, maxit = maxit, h = h, ...)
            directions <- x.i %*% s.phd.cur$beta

            mses <- numeric(length(h))
            for (hv in 1:length(h))
            {
                loc.mod <- locfit.raw(x = directions, y = y[-which], alpha = h[hv], deg = 2, ...)
                mses[hv] <- (y[which] - predict(loc.mod, x[which,,drop=FALSE] %*% s.phd.cur$beta)) ^ 2
            }

            gcvs[i, d]   <- min(mses)

        }
        if (i > 1) print(colMeans(gcvs, na.rm=TRUE))
    }

    gcv.res    <- colMeans(gcvs)
    best.d     <- which.min(gcv.res)
    best.model <- phd(x, y, d = best.d)

    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcv.res
    ret$gcv.mat <- gcvs
    ret
}

semi.phd.hier.separate.dim.select.k <- function(x.list, y, max.d = rep(1L, 3L), k = 5L, maxit = 10L, h = NULL, ...)
{
    best.gcv <- 1e99
    d.mat <- expand.grid(1:max.d[1], 1:max.d[2], 0:max.d[3])
    n.D <- nrow(d.mat)
    gcvs <- array(NA, dim = c(k, n.D))
    #models <- vector(mode = "list", length = max.d)
    #nobs <- nrow(x)


    nobs.vec <- unlist(lapply(x.list, nrow))
    c.nobs <- c(0, cumsum(nobs.vec))

    foldid.list <- lapply(1:length(x.list), function(i) sample(rep(seq(k), length = nobs.vec[i])))

    y.list <- vector(mode = "list", length = length(x.list))
    y.list[[1]] <- y[(c.nobs[1] + 1):c.nobs[2]]
    y.list[[2]] <- y[(c.nobs[2] + 1):c.nobs[3]]
    y.list[[3]] <- y[(c.nobs[3] + 1):c.nobs[4]]

    for (i in 1:k)
    {
        which.list <- lapply(1:length(x.list) , function(idx) which(foldid.list[[idx]] == i))
        x.i        <- lapply(1:length(x.list) , function(idx) x.list[[idx]][-which.list[[idx]],])
        y.i        <- lapply(1:length(x.list) , function(idx) y.list[[idx]][-which.list[[idx]]])
        y.i.vec <- unlist(y.i)

        x.te       <- lapply(1:length(x.list) , function(idx) x.list[[idx]][which.list[[idx]],])
        y.te       <- lapply(1:length(x.list) , function(idx) y.list[[idx]][which.list[[idx]] ])
        #phd.cur <- phd(x.i, y[-which], d = max.d)
        print(i)
        for (d in 1:n.D)
        {
            if (i == nobs) print(d.mat[d,])

            print(d.mat[d,])
            d.cur <- as.vector(data.matrix(d.mat[d,]))
            print(d.cur)
            s.phd.cur <- semi.phd.hier.separate(x.i, y.i.vec, d = d.cur, maxit = maxit, h = h, ...)


            mses <- numeric(length(h))
            for (gr in 1:length(x.list))
            {
                directions      <- x.i[[gr]]  %*% s.phd.cur$beta[[gr]]
                directions.test <- x.te[[gr]] %*% s.phd.cur$beta[[gr]]
                for (hv in 1:length(h))
                {
                    loc.mod <- locfit.raw(x = directions, y = y.i[[gr]], alpha = h[hv], deg = 2, ...)
                    mses[hv] <- mses[hv] + sum((y.te[[gr]] - predict(loc.mod, directions.test)) ^ 2)
                }
            }

            mses <- mses / sum(nobs.vec)

            gcvs[i, d]   <- min(mses)

        }
        if (i > 1) print(colMeans(gcvs, na.rm=TRUE))
    }

    gcv.res    <- colMeans(gcvs)
    best.d     <- as.vector(data.matrix(d.mat[which.min(gcv.res),]))
    #best.model <- phd(x, y, d = best.d)

    #ret   <- best.model
    ret <- list(best.d = best.d, d.mat = d.mat, gcvs = gcv.res, gcv.mat = gcvs)
    #ret$best.d <- best.d
    #ret$d.mat  <- d.mat
    #ret$gcvs   <- gcv.res
    #ret$gcv.mat <- gcvs
    ret
}


semi.phd.hier.separate.dim.cv.k <- function(x.list, y, d = rep(1L, 2L), k = 5L, maxit = 10L, h = NULL, ...)
{
    best.gcv <- 1e99

    gcvs <- rep(NA, k)

    nobs.vec <- unlist(lapply(x.list, nrow))
    c.nobs <- c(0, cumsum(nobs.vec))

    foldid.list <- lapply(1:length(x.list), function(i) sample(rep(seq(k), length = nobs.vec[i])))

    y.list <- vector(mode = "list", length = length(x.list))
    y.list[[1]] <- y[(c.nobs[1] + 1):c.nobs[2]]
    y.list[[2]] <- y[(c.nobs[2] + 1):c.nobs[3]]
    y.list[[3]] <- y[(c.nobs[3] + 1):c.nobs[4]]

    for (i in 1:k)
    {
        which.list <- lapply(1:length(x.list) , function(idx) which(foldid.list[[idx]] == i))
        x.i        <- lapply(1:length(x.list) , function(idx) x.list[[idx]][-which.list[[idx]],])
        y.i        <- lapply(1:length(x.list) , function(idx) y.list[[idx]][-which.list[[idx]]])
        y.i.vec <- unlist(y.i)

        x.te       <- lapply(1:length(x.list) , function(idx) x.list[[idx]][which.list[[idx]],])
        y.te       <- lapply(1:length(x.list) , function(idx) y.list[[idx]][which.list[[idx]] ])

        print(i)

        d.cur <- as.vector(data.matrix(d))
        print(d.cur)
        s.phd.cur <- semi.phd.hier.separate(x.i, y.i.vec, d = d.cur, maxit = maxit, h = h, ...)


        mses <- numeric(length(h))
        for (gr in 1:length(x.list))
        {
            directions      <- x.i[[gr]]  %*% s.phd.cur$beta[[gr]]
            directions.test <- x.te[[gr]] %*% s.phd.cur$beta[[gr]]
            for (hv in 1:length(h))
            {
                loc.mod <- locfit.raw(x = directions, y = y.i[[gr]], alpha = h[hv], deg = 2, ...)
                mses[hv] <- mses[hv] + sum((y.te[[gr]] - predict(loc.mod, directions.test)) ^ 2)
            }
        }

        mses <- mses / sum(nobs.vec)

        gcvs[i]   <- min(mses)


        if (i > 1) print(mean(gcvs, na.rm=TRUE))
    }

    gcv.res    <- mean(gcvs)
    #ret   <- best.model
    ret <- list(d = d, gcvs = gcv.res, gcv.mat = gcvs)
    ret
}

phd.dim.select.gcv <- function(x, y, max.d = 5L, h = NULL, ...)
{
    best.gcv <- 1e99
    gcvs <- array(NA, dim = c(nrow(x), max.d))
    models <- vector(mode = "list", length = max.d)
    nobs <- nrow(x)
    phd.cur <- phd(x, y, d = max.d)

    gcv.res <- numeric(max.d)
    for (d in 1:max.d)
    {
        directions <- x %*% phd.cur$beta.hat[,1:d]

        gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 2, ...)[4])
        gcv.res[d] <- min(gcv.vals)
        print(which.min(gcv.vals))
    }

    best.d     <- which.min(gcv.res)
    best.model <- phd(x, y, d = best.d)

    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcv.res
    ret
}

phd.dim.select <- function(x, y, max.d = 5L, h = NULL, ...)
{
    best.gcv <- 1e99
    gcvs <- array(NA, dim = c(nrow(x), max.d))
    models <- vector(mode = "list", length = max.d)
    nobs <- nrow(x)
    for (i in 1:nobs)
    {
        x.i     <- x[-i,]
        phd.cur <- phd(x.i, y[-i], d = max.d)
        print(i)
        for (d in 1:max.d)
        {
            if (i == nobs) print(d)
            directions <- x.i %*% phd.cur$beta.hat[,1:d]

            mses <- numeric(length(h))
            for (hv in 1:length(h))
            {
                loc.mod <- locfit.raw(x = directions, y = y[-i], alpha = h[hv], deg = 2, ...)
                mses[hv] <- (y[i] - predict(loc.mod, x[i,,drop=FALSE] %*% phd.cur$beta.hat[,1:d])) ^ 2
            }

            gcvs[i, d]   <- min(mses)

        }
        if (i > 1) print(colMeans(gcvs, na.rm=TRUE))
    }

    gcv.res    <- colMeans(gcvs)
    best.d     <- which.min(gcv.res)
    best.model <- phd(x, y, d = best.d)

    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcv.res
    ret$gcv.mat <- gcvs
    ret
}

phd.dim.select.k <- function(x, y, max.d = 5L, h = NULL, k = 10L, ...)
{
    best.gcv <- 1e99
    gcvs <- array(NA, dim = c(k, max.d))
    models <- vector(mode = "list", length = max.d)
    nobs <- nrow(x)

    foldid <- sample(rep(seq(k), length = nobs))

    for (i in 1:k)
    {
        which = which(foldid == i)
        x.i     <- x[-which,]
        phd.cur <- phd(x.i, y[-which], d = max.d)
        print(i)
        for (d in 1:max.d)
        {
            if (i == nobs) print(d)
            directions <- x.i %*% phd.cur$beta.hat[,1:d]

            mses <- numeric(length(h))
            for (hv in 1:length(h))
            {
                loc.mod <- locfit.raw(x = directions, y = y[-which], alpha = h[hv], deg = 2, ...)
                mses[hv] <- mean((y[which] - predict(loc.mod, x[which,,drop=FALSE] %*% phd.cur$beta.hat[,1:d])) ^ 2)
            }

            gcvs[i, d]   <- min(mses)

        }
        if (i > 1) print(colMeans(gcvs, na.rm=TRUE))
    }

    gcv.res    <- colMeans(gcvs)
    best.d     <- which.min(gcv.res)
    best.model <- phd(x, y, d = best.d)

    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcv.res
    ret$gcv.mat <- gcvs
    ret
}



hier.phd.dim.select.k <- function(x.list, y, max.d = 3L, h = NULL, k = 10L, ...)
{
    best.gcv <- 1e99

    models <- vector(mode = "list", length = max.d)
    nobs <- nrow(x.list[[1]])
    nvars <- ncol(x.list[[1]])


    nobs.vec <- unlist(lapply(x.list, nobs))

    foldid.list <- lapply(1:length(nobs.vec), function(i) sample(rep(seq(k), length = nobs.vec[i])))

    y.list <- vector(mode = "list", length = length(x.list))

    c.nobs <- c(0, cumsum(nobs.vec))

    for (i in 1:length(y.list)) y.list[[i]] <- y[(c.nobs[i] + 1):c.nobs[i]]

    d.vec <- c(1, 1, 0)
    for (g in 1:3)
    {
        for (i in 1:k)
        {
            which.list <- lapply(foldid.list, function(fid) which(fid == i))
            x.i     <- lapply(1:length(nobs.vec), function(i) x[-which.list[[i]],])
            y.i     <- lapply(1:length(y.list), function(id) y.list[[id]][-which.list[[id]]] )
            y.test.list <- lapply(1:length(y.list), function(id) y.list[[id]][which.list[[id]]] )
            x.test.list <- lapply(1:length(nobs.vec), function(i) x[which.list[[i]],])
            print(i)
            if (g < 3)
            {
                loop.d <- 1:max.d
            } else
            {
                loop.d <- 0:(max.d - 1)
            }
            ct <- 0
            gcvs <- array(NA, dim = c(k, length(loop.d)))
            for (d in loop.d)
            {
                ct <- ct + 1
                d.vec.cur <- d.vec
                d.vec.cur[g] <- d
                if (i == nobs) print(d)

                phd.cur <- hier.phd(x.i, y.i, d = d.vec.cur)

                #directions <- x.i %*% phd.cur$beta.hat

                mses <- array(NA, dim = c(length(h), 3) )

                directions.phd.1 <- x.i[[1]] %*% Re(phd.cur$beta.hat)[1:nvars,1]
                directions.phd.2 <- x.i[[2]] %*% Re(phd.cur$beta.hat)[(1 + nvars):(2 * nvars),2]
                directions.phd.3 <- x.i[[3]] %*% Re(phd.cur$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

                directions.phd.test.1 <- x.test.list[[1]] %*% Re(phd.cur$beta.hat)[1:nvars,1]
                directions.phd.test.2 <- x.test.list[[2]] %*% Re(phd.cur$beta.hat)[(1 + nvars):(2 * nvars),2]
                directions.phd.test.3 <- x.test.list[[3]] %*% Re(phd.cur$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

                for (hv in 1:length(h))
                {

                    mse.tmp <- numeric(3)

                    loc.mod1 <- locfit.raw(x = directions.phd.1, y = y.i[[1]], alpha = h[hv], deg = 2, ...)
                    loc.mod2 <- locfit.raw(x = directions.phd.2, y = y.i[[2]], alpha = h[hv], deg = 2, ...)
                    loc.mod2 <- locfit.raw(x = directions.phd.3, y = y.i[[3]], alpha = h[hv], deg = 2, ...)

                    mses[hv, 1] <- (y.test.list[[1]] - predict(loc.mod1, directions.phd.test.1)) ^ 2
                    mses[hv, 2] <- (y.test.list[[2]] - predict(loc.mod2, directions.phd.test.2)) ^ 2
                    mses[hv, 3] <- (y.test.list[[3]] - predict(loc.mod3, directions.phd.test.3)) ^ 2
                }

                gcvs[i, ct]   <- mean(apply(mses, 2, min))

            }
            if (i > 1) print(colMeans(gcvs, na.rm=TRUE))
        }
        gcv.res    <- colMeans(gcvs)
        d.vec[g]   <- loop.d[which.min(gcv.res)]
    }

    gcv.res    <- colMeans(gcvs)
    best.d     <- which.min(gcv.res)
    best.model <- phd(x, y, d = best.d)

    ret   <- best.model
    ret$d <- best.d
    ret$gcvs <- gcv.res
    ret$gcv.mat <- gcvs
    ret
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

    slver <-   optim(par     = beta.init, # beta.init[(d+1):nrow(beta.init),],
                     fn      = est.eqn,
                     #gr      = est.eqn.grad,
                     method  = "L-BFGS",
                     control = list(maxit = maxit, reltol = 1e-8))

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
                     #gr      = est.eqn.grad,
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




semi.phd.hier.separate <- function(x.list, y, d = rep(1L, 3L), maxit = 10L, h = NULL, B = NULL, vic = FALSE, ...)
{
    p <- ncol(x.list[[1]])

    d <- as.vector(d)
    names(d) <- NULL

    nobs.vec  <- unlist(lapply(x.list, nrow))
    nvars.vec <- unlist(lapply(x.list, ncol))
    pp <- nvars.vec[1]

    D <- sum(d)
    d.vic <- d + 1
    D.vic <- sum(d.vic)

    cum.d     <- c(0, cumsum(d))
    cum.d.vic <- c(0, cumsum(d.vic))

    cov <- lapply(x.list, cov)

    x <- as.matrix(bdiag(x.list))
    cov.b <- cov(as.matrix(x))
    eig.cov.b <- eigen(cov.b)
    sqrt.inv.cov.b <- eig.cov.b$vectors %*% diag(1 / sqrt(eig.cov.b$values)) %*% t(eig.cov.b$vectors)
    x.tilde.b <- scale(x, scale = FALSE) %*% sqrt.inv.cov.b

    sqrt.inv.cov <- lapply(1:length(x.list), function(i) {
        eig.cov <- eigen(cov[[i]])
        eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    })

    x.tilde <- lapply(1:length(x.list), function(i) {
        scale(x.list[[i]], scale = FALSE) %*% sqrt.inv.cov[[i]]
    })


    constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                                cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                        t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                                cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                        t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                                 cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )



    constr            <- t(cbind(diag(p), -diag(p)))
    Proj.constr       <- constr %*% solve(crossprod(constr), t(constr))
    Ortho.Proj.constr <- diag(ncol(Proj.constr)) - Proj.constr

    #for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    strat.id <- unlist(lapply(1:length(x.list), function(id) rep(id, nrow(x.list[[id]]))))

    V.hat <- crossprod(x.tilde.b, drop(scale(y, scale = FALSE)) * x.tilde.b) / nrow(x.tilde.b)


    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        print(d)
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

    beta.init.cov <- t(t(beta.init) %*% sqrt.inv.cov.b)
    for (s in 1:num.strata)
    {
        strata.idx <- which(strat.id == unique.strata[s])
        dir.cur    <- x.tilde.b[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
            beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
        #print(apply(dir.cur, 2, sd))
        dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

        gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
                                                 y = y[strata.idx],
                                                 kern = "gauss",
                                                 alpha = hv, deg = 3, ...)[4])
        best.h.vec[s]     <- h[which.min(gcv.vals)]
    }

    print(best.h.vec)

    est.eqn <- function(beta.vec)
    {
        beta.mat.list <- vector(mode = "list", length = 3L)

        beta.mat.list[[1]] <- matrix(beta.vec[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
        beta.mat.list[[2]] <- matrix(beta.vec[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
        beta.mat.list[[3]] <- matrix(beta.vec[(cum.d[3] * p + 1):length(beta.vec)], ncol = cum.d[4] )

        #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         kern = "gauss",
            #                                         alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]
            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "gauss",
            #                          alpha = best.h, deg = 3, ...)


            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                     kern = "gauss",
                                     alpha = c(best.h), deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }


        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            norm(crossprod(x.tilde[[i]], resid * x.tilde[[i]]), type = "F") ^ 2 / (nobs.vec[i] ^ 2)
            })))
        lhs
    }


    est.eqn.vic <- function(beta.vec)
    {
        beta.mat.list <- vector(mode = "list", length = 3L)

        beta.mat.list[[1]] <- matrix(beta.vec[(cum.d.vic[1] * p + 1):(cum.d.vic[2] * p)], ncol = d.vic[1])
        beta.mat.list[[2]] <- matrix(beta.vec[(cum.d.vic[2] * p + 1):(cum.d.vic[3] * p)], ncol = d.vic[2])
        beta.mat.list[[3]] <- matrix(beta.vec[(cum.d.vic[3] * p + 1):length(beta.vec)], ncol = cum.d.vic[4] )

        #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]
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


        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            norm(crossprod(x.tilde[[i]], resid * x.tilde[[i]]), type = "F") ^ 2 / (nobs.vec[i])
        })))
        lhs
    }

    est.eqn.grad <- function(beta.vec)
    {
        grad.full     <- grad(est.eqn, beta.vec, method = "simple")

        grad.mat.list <- vector(mode = "list", length = 3L)

        grad.mat.list[[1]] <- matrix(grad.full[(cum.d[1] * p + 1):(cum.d[2] * p)],   ncol = d[1])
        grad.mat.list[[2]] <- matrix(grad.full[(cum.d[2] * p + 1):(cum.d[3] * p)],   ncol = d[2])
        grad.mat.list[[3]] <- matrix(grad.full[(cum.d[3] * p + 1):length(beta.vec)], ncol = cum.d[4] )



        ## projection onto stiefel manifold
        #svd1 <- svd(grad.mat.list[[1]])
        #svd2 <- svd(grad.mat.list[[2]])
        #svd3 <- svd(grad.mat.list[[3]])

        #grad.mat.list[[1]] <- tcrossprod(svd1$u, svd1$v)
        #grad.mat.list[[2]] <- tcrossprod(svd2$u, svd2$v)
        #grad.mat.list[[3]] <- tcrossprod(svd3$u, svd3$v)

        print(str(grad.mat.list[[1]]))
        print(str(grad.mat.list[[2]]))
        print(str(grad.mat.list[[3]]))
        print(str(rbind(grad.mat.list[[1]], grad.mat.list[[3]][,1:ncol(grad.mat.list[[1]]), drop=FALSE])))

        gA <- Ortho.Proj.constr %*% rbind(grad.mat.list[[1]], grad.mat.list[[3]][,1:d[1], drop=FALSE])
        gB <- Ortho.Proj.constr %*% rbind(grad.mat.list[[2]], grad.mat.list[[3]][,(d[1] + 1):cum.d[3], drop=FALSE])

        print("proj")
        print(str(gA))
        print(str(gB))


        gA  <- gA[1:nrow(grad.mat.list[[1]]),,drop=FALSE]
        gB  <- gB[1:nrow(grad.mat.list[[2]]),,drop=FALSE]
        if (cum.d[3] < cum.d[4])
        {
            gAB <- cbind(gA, gB, grad.mat.list[[3]][, (cum.d[3] + 1):cum.d[4]])
        } else
        {
            gAB <- cbind(gA, gB)
        }

        print(str(gAB))
        print(head(gAB))

        c(as.vector(gA), as.vector(gB), as.vector(gAB))
    }

    est.eqn.vic.grad <- function(beta.vec)
    {
        grad.full     <- grad(est.eqn.vic, beta.vec, method = "simple")

        grad.mat.list <- vector(mode = "list", length = 3L)

        grad.mat.list[[1]] <- matrix(grad.full[(cum.d.vic[1] * p + 1):(cum.d.vic[2] * p)],   ncol = d.vic[1])
        grad.mat.list[[2]] <- matrix(grad.full[(cum.d.vic[2] * p + 1):(cum.d.vic[3] * p)],   ncol = d.vic[2])
        grad.mat.list[[3]] <- matrix(grad.full[(cum.d.vic[3] * p + 1):length(beta.vec)], ncol = cum.d.vic[4] )



        ## projection onto stiefel manifold
        #svd1 <- svd(grad.mat.list[[1]])
        #svd2 <- svd(grad.mat.list[[2]])
        #svd3 <- svd(grad.mat.list[[3]])

        #grad.mat.list[[1]] <- tcrossprod(svd1$u, svd1$v)
        #grad.mat.list[[2]] <- tcrossprod(svd2$u, svd2$v)
        #grad.mat.list[[3]] <- tcrossprod(svd3$u, svd3$v)

        print(str(grad.mat.list[[1]]))
        print(str(grad.mat.list[[2]]))
        print(str(grad.mat.list[[3]]))
        print(str(rbind(grad.mat.list[[1]], grad.mat.list[[3]][,1:ncol(grad.mat.list[[1]]), drop=FALSE])))

        gA <- Ortho.Proj.constr %*% rbind(grad.mat.list[[1]], grad.mat.list[[3]][,1:d.vic[1], drop=FALSE])
        gB <- Ortho.Proj.constr %*% rbind(grad.mat.list[[2]], grad.mat.list[[3]][,(d.vic[1] + 1):cum.d.vic[3], drop=FALSE])

        print("proj")
        print(str(gA))
        print(str(gB))


        gA  <- gA[1:nrow(grad.mat.list[[1]]),,drop=FALSE]
        gB  <- gB[1:nrow(grad.mat.list[[2]]),,drop=FALSE]
        if (cum.d.vic[3] < cum.d.vic[4])
        {
            gAB <- cbind(gA, gB, grad.mat.list[[3]][, (cum.d.vic[3] + 1):cum.d.vic[4]])
        } else
        {
            gAB <- cbind(gA, gB)
        }

        print(str(gAB))
        print(head(gAB))

        c(as.vector(gA), as.vector(gB), as.vector(gAB))
    }

    beta.init.vec <- numeric(nvars.vec[1] * d[1] + nvars.vec[2] * d[2] + nvars.vec[3] * sum(d))

    beta.init.vec[1:(nvars.vec[1] * d[1])] <- as.vector(beta.init[1:nvars.vec[1], 1:d[1]])
    beta.init.vec[(nvars.vec[1] * d[1] + 1):
                      (nvars.vec[1] * d[1] + nvars.vec[2] * d[2])] <- as.vector(beta.init[(nvars.vec[1] + 1):
                                                                                              (nvars.vec[1] + nvars.vec[2]), (d[1] + 1):(d[1] + d[2]) ])
    beta.init.vec[(nvars.vec[1] * d[1] + nvars.vec[2] * d[2] + 1):
                      (length(beta.init.vec))] <- as.vector(beta.init[(nvars.vec[1] + nvars.vec[2] + 1):
                                                                                              (nrow(beta.init)), ])

    print("two eq:")
    print(head(beta.init[1:nvars.vec[1], 1:d[1]]))
    print(head(beta.init[(nvars.vec[1] + 1):
                             (nvars.vec[1] + nvars.vec[2]), (d[1] + 1):(d[1] + d[2]) ]))
    print(head(beta.init[(nvars.vec[1] + nvars.vec[2] + 1):
                             (nrow(beta.init)), ]))

    slver <-   optim(par     = beta.init.vec, # beta.init[(d+1):nrow(beta.init),],
                     fn      = est.eqn,
                     gr      = est.eqn.grad,
                     method  = "L-BFGS",
                     control = list(maxit = maxit, factr = 1e-10))

    vic <- slver$value * nobs + log(nobs) * nvars * (d[1] + d[2] + d[3])

    beta <- beta.init

    beta.mat.list <- vector(mode = "list", length = 3L)

    beta.mat.list[[1]] <- matrix(slver$par[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
    beta.mat.list[[2]] <- matrix(slver$par[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
    beta.mat.list[[3]] <- matrix(slver$par[(cum.d[3] * p + 1):length(slver$par)], ncol = cum.d[4] )



    model.list <- vector(mode = "list", length = 3)


    sse.vec <- mse.vec <- numeric(3)
    for (m in 1:3)
    {
        strata.idx <- which(strat.id == unique.strata[m])
        best.h     <- best.h.vec[m]
        dir.cur <- x.list[[m]] %*% beta.mat.list[[m]]
        model.list[[m]] <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                 kern = "gauss",
                                 alpha = best.h, deg = 3, ...)

        fitted.vals <- fitted(model.list[[m]])

        sse.vec[m] <-  sum((y - fitted.vals) ^ 2)
        mse.vec[m] <- mean((y - fitted.vals) ^ 2)
    }


    #beta.semi <- matrix(slver$par, ncol = D)
    #t(t(matrix(slver$par, ncol = d * length(constraints))) %*% sqrt.inv.cov)


    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.mat.list, beta.init = beta, solver.obj = slver,
         #beta.rand.init = t(t(beta.rand.init) %*% sqrt.inv.cov),
         cov = cov, sqrt.inv.cov = sqrt.inv.cov,
         vic = vic, sse = sse.vec, mse = mse.vec)
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
    list(beta.hat = beta.hat, eta.hat = eta.hat, eigenvalues = eig.V$values)
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

    P1 <- B %*% solve(crossprod(B), t(B))
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

