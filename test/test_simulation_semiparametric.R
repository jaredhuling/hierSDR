
library(hierSDR)
library(Matrix)

fig.path <- "C:/Users/Jared/Dropbox/ACO/rehospitalization/sdr_hierarchical/figures/"

hier.sir <- function(x.list, y, h = 10L, d = rep(2L, 3L), slice.ind = NULL)
{
    p <- ncol(x.list[[1]])
    x <- bdiag(x.list)
    cov <- cov(as.matrix(x))
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    D <- sum(d)
    cum.d <- c(0, cumsum(d))

    constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                                cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                        t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                                cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                        t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                                 cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )

    for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

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

    beta.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            eig.c <- eigen((diag(ncol(Pc)) - Pc) %*% V.hat )
            eta.hat <- eig.c$vectors[,1:d[c]]
            beta.list[[c]] <- t(t(eta.hat) %*% sqrt.inv.cov)
        }
    }

    list(beta.hat = do.call(cbind, beta.list), cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}



hier.phd <- function(x.list, y, d = rep(2L, 3L))
{
    p <- ncol(x.list[[1]])
    x <- bdiag(x.list)
    cov <- cov(as.matrix(x))
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    D <- sum(d)
    cum.d <- c(0, cumsum(d))

    constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                                cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                        t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                                cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                        t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                                 cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )

    for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x)

    beta.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            eig.c <- eigen((diag(ncol(Pc)) - Pc) %*% V.hat )
            eta.hat <- eig.c$vectors[, 1:d[c] ]
            beta.list[[c]] <- t(t(eta.hat) %*% sqrt.inv.cov)
        }
    }

    list(beta.hat = do.call(cbind, beta.list), cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}

hier.semi.phd <- function(x.list, y, d = rep(2L, 3L), maxit = 10L, h = NULL, ...)
{
    p <- ncol(x.list[[1]])
    x <- as.matrix(bdiag(x.list))
    cov <- cov(as.matrix(x))
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    #x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    D <- sum(d)
    cum.d <- c(0, cumsum(d))

    constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                                cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                        t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                                cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                        t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                                 cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )

    #for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    strat.id <- unlist(lapply(1:length(x.list), function(id) rep(id, nrow(x.list[[id]]))))

    s.phd <- semi.phd.for.hier(x = x, y = y, d = d,
                               maxit = maxit, h = h,
                               strat.id = strat.id,
                               constraints = constraints,
                               ...)
    beta.semi.phd <- s.phd$beta

    beta.list <- beta.W.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        if (d[c] > 0)
        {
            AA <- constraints[[c]]
            adj.fact   <- AA %*% solve(crossprod(AA), crossprod(AA, beta.semi.phd[, (cum.d[c] + 1):cum.d[c + 1]]))
            #adj.fact.W <- AA %*% solve(crossprod(AA), crossprod(AA, beta.W.semi.phd[, ((c-1) * d + 1):(c * d)]))
            beta.list[[c]]   <- beta.semi.phd[,   (cum.d[c] + 1):cum.d[c + 1]] - adj.fact
            #beta.W.list[[c]] <- beta.W.semi.phd[, ((c-1) * d + 1):(c * d)] - adj.fact.W
        }
    }
    beta.sph <- do.call(cbind, beta.list)

    #s.phd.W <- semi.phd.for.hier.W(x = x, y = y,
    #                               W = solve(s.phd$W),
    #                               init = beta.sph,
    #                               d = d * 3,
    #                               maxit = maxit, h = h,
    #                               strat.id = strat.id,
    #                               constraints = constraints,
    #                               ...)
    #beta.W.semi.phd <- s.phd.W$beta

    #for (c in 1:length(constraints))
    #{
    #    AA <- constraints[[c]]
    #    #adj.fact   <- AA %*% solve(crossprod(AA), crossprod(AA, beta.semi.phd[, ((c-1) * d + 1):(c * d)]))
    #    adj.fact.W <- AA %*% solve(crossprod(AA), crossprod(AA, beta.W.semi.phd[, ((c-1) * d + 1):(c * d)]))
    #    #beta.list[[c]]   <- beta.semi.phd[,   ((c-1) * d + 1):(c * d)] - adj.fact
    #    beta.W.list[[c]] <- beta.W.semi.phd[, ((c-1) * d + 1):(c * d)] - adj.fact.W
    #}


    list(beta.hat = do.call(cbind, beta.list),
         #beta.hat.W = do.call(cbind, beta.W.list),
         cov = cov,
         sqrt.inv.cov = sqrt.inv.cov,
         beta.unconstrained = beta.semi.phd,
         solver.obj = s.phd$solver.obj)
}


# simulation parameters
nsims     <- 50
nobs.vec  <- c(250, 500, 1000, 2000)[1:3]
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 2

set.seed(123)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 7))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "semi hier phd", "sir separate", "phd separate", "semi phd separate", "semi hier phd 2")

sim.res.list <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list <- rep(list(sim.res.dir), length(nobs.vec))
sim.subsp.angle.res.list <- sim.direction.res.list

for (n in 3:length(nobs.vec))
{
    nobs <- nobs.vec[n]

    for (s in 6:nsims)
    {

        x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        #x.list      <- replicate(3, list(matrix(c(rnorm(nobs * nvars/2), rbinom(nobs * nvars/2, 1, 0.5)), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(c(rnorm(nobs.test * nvars/2), rbinom(nobs.test * nvars/2, 1, 0.5)), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)


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

        var(y.true) / var(y - y.true)

        hier.sdr      <- hier.sir(x.list, y,  d = rep(1, 3), h = 30L)
        hier.sdr.phd  <- hier.phd(x.list, y,  d = rep(1, 3))
        sdr.sir       <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd       <- phd(as.matrix(x), y, d = 1 * 3)

        semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = rep(1, 3), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)
        semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = rep(1, 3), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)


        sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)

        s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)

        cons <- t(cbind(diag(nvars), -diag(nvars)))
        Pc   <- cons %*% solve(crossprod(cons), t(cons))

        betaA  <- (diag(nvars * 2) - Pc) %*%  rbind(s.phd.1$beta, s.phd.3$beta[,1,drop=FALSE])
        betaB  <- (diag(nvars * 2) - Pc) %*%  rbind(s.phd.2$beta, s.phd.3$beta[,2,drop=FALSE])
        betaAB <- cbind(betaA[1:nvars,], betaB[1:nvars,], s.phd.3$beta[,3,drop=FALSE])

        est.eqn.tmp <- function(beta.vec)
        {
            h = seq(1, 50, length.out = 25)
            x.tt       <- x.list[[1]]
            cov <- cov(x.tt)
            eig.cov <- eigen(cov)
            sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
            x.tildet <- scale(x.tt, scale = FALSE) %*% sqrt.inv.cov

            #beta.mat   <- rbind(beta.init[1:1,], matrix(beta.vec, ncol = 1))
            beta.mat   <- matrix(beta.vec, ncol = 1)
            directions <- x.tildet %*% beta.mat
            gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y[1:nobs],
                                                     alpha = hv, maxk = 250, deg = 3)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = directions, y = y[1:nobs], alpha = best.h, maxk = 250, deg = 3)


            Ey.given.xbeta <- fitted(locfit.mod)

            resid <- drop(y[1:nobs] - Ey.given.xbeta)
            lhs   <- norm(crossprod(x.tildet, resid * x.tildet), type = "F") ^ 2 / (nobs ^ 2)
            lhs
        }

        est.eqn.grad.tmp <- function(beta.vec)
        {
            h = seq(1, 50, length.out = 25)
            x.tt       <- x.list[[1]]
            cov <- cov(x.tt)
            eig.cov <- eigen(cov)
            sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
            x.tildet <- scale(x.tt, scale = FALSE) %*% sqrt.inv.cov

            #beta.mat   <- rbind(beta.init[1:1,], matrix(beta.vec, ncol = 1))
            beta.mat   <- matrix(beta.vec, ncol = 1)
            directions <- x.tildet %*% beta.mat
            gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y[1:nobs], alpha = hv, maxk = 250)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod       <- locfit.raw(x = directions, y = y[1:nobs],
                                           alpha = best.h,            maxk = 250,
                                           deg   = 3)
            locfit.mod.deriv <- locfit.raw(x = directions, y = y[1:nobs],
                                           alpha = best.h, deriv = 1, maxk = 250,
                                           deg   = 3)


            Ey.given.xbeta       <- fitted(locfit.mod)
            Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

            resid <- drop(y[1:nobs] - Ey.given.xbeta)
            psi      <- crossprod(x.tildet, resid * x.tildet) / nobs
            psi.grad <- -crossprod(x.tildet, (Ey.given.xbeta.deriv * rowSums(x.tildet ^ 2) ) ) / nobs
            gradient <- 2 * t(psi) %*% psi.grad
            rep(drop(gradient), 1)
        }


        est.eqn.tmp2 <- function(beta.vec, idx = 1L)
        {
            h = exp(seq(log(0.25), log(15), length.out = 25))
            x.tt       <- x.list[[1]]
            cov <- cov(x.tt)
            eig.cov <- eigen(cov)
            sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
            x.tildet <- scale(x.tt, scale = FALSE) %*% sqrt.inv.cov

            #beta.mat   <- rbind(beta.init[1:1,], matrix(beta.vec, ncol = 1))
            beta.mat   <- matrix(beta.vec, ncol = 1)
            directions <- x.tildet %*% beta.mat
            gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y[1:nobs], alpha = hv, deg = 3, maxk = 250)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = directions, y = y[1:nobs], alpha = best.h, deg = 3, maxk = 250)


            Ey.given.xbeta <- fitted(locfit.mod)

            resid <- drop(y[1:nobs] - Ey.given.xbeta)
            #lhs   <- (crossprod(x.tildet, resid * x.tildet) / nobs)[1,1]
            lhs <- sum(resid) / nobs
            #lhs
            resid[idx]
        }

        est.eqn.grad.tmp2 <- function(beta.vec, idx = 1L)
        {
            h = exp(seq(log(0.25), log(15), length.out = 25))
            x.tt       <- x.list[[1]]
            cov <- cov(x.tt)
            eig.cov <- eigen(cov)
            sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
            x.tildet <- scale(x.tt, scale = FALSE) %*% sqrt.inv.cov

            #beta.mat   <- rbind(beta.init[1:1,], matrix(beta.vec, ncol = 1))
            beta.mat   <- matrix(beta.vec, ncol = 1)
            directions <- x.tildet %*% beta.mat
            gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y[1:nobs], alpha = hv, deg = 3, maxk = 250)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = directions, y = y[1:nobs], alpha = best.h, deg = 3, maxk = 250)
            locfit.mod.deriv <- locfit.raw(x = directions, y = y[1:nobs], alpha = best.h, deriv = 1, deg = 3, maxk = 250)


            Ey.given.xbeta       <- fitted(locfit.mod)
            Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

            resid <- drop(y[1:nobs] - Ey.given.xbeta)
            #psi      <- (crossprod(x.tildet, resid * x.tildet) / nobs)
            #psi.grad <- -crossprod(x.tildet, (Ey.given.xbeta.deriv * rowSums(x.tildet ^ 2) ) ) / nobs
            #gradient <- 2 * t(psi) %*% psi.grad
            gradient <- -colSums(Ey.given.xbeta.deriv * x.tildet) / nobs
            #gradient <- -Ey.given.xbeta.deriv[10] * x.tildet[10,]
            #gradient <- -Ey.given.xbeta.deriv[idx] * x.tildet[idx,]
            rep(drop(gradient), 1)
        }


        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }

        tmp.fun.grad <- function(x.vec)
        {
            set.seed(123)
            x.tmp <- matrix(rnorm(100), ncol = 1)
            y.tmp <- drop(sin(x.tmp * pi) * x.tmp ^ 2) + rnorm(100, sd = 0.1)

            h = seq(0.5, 15, length.out = 25)
            gcv.vals   <- sapply(h, function(hv) gcv(x = x.tmp, y = y.tmp, alpha = hv, deg = 3, maxk = 250)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = x.tmp, y = y.tmp,       alpha = best.h, deg = 3, maxk = 250)
            locfit.mod.deriv <- locfit.raw(x = x.tmp, y = y.tmp, alpha = best.h, deriv = 1, deg = 3, maxk = 250)
            #Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

            predict(locfit.mod.deriv, matrix(x.vec, ncol = 1))
        }

        tmp.fun <- function(x.vec)
        {
            set.seed(123)
            x.tmp <- matrix(rnorm(100), ncol = 1)
            y.tmp <- drop(sin(x.tmp * pi) * x.tmp ^ 2) + rnorm(100, sd = 0.1)

            h = seq(0.5, 15, length.out = 25)
            gcv.vals   <- sapply(h, function(hv) gcv(x = x.tmp, y = y.tmp, alpha = hv, deg = 3, maxk = 250)[4])
            best.h     <- h[which.min(gcv.vals)]
            locfit.mod <- locfit.raw(x = x.tmp, y = y.tmp,       alpha = best.h, deg = 3, maxk = 250)
            #locfit.mod.deriv <- locfit.raw(x = x.tmp, y = y.tmp, alpha = best.h, deriv = 1, maxk = 250)
            #Ey.given.xbeta.deriv <- fitted(locfit.mod.deriv)

            predict(locfit.mod, matrix(x.vec, ncol = 1))
        }


        #round(gr.a <- est.eqn.grad.tmp( s.phd.1$solver.obj$par), 2)
        #round(gr.n <- grad(est.eqn.tmp, s.phd.1$solver.obj$par), 2)
        #cor(gr.a, gr.n)

        ## psi grad
        #grad.analytical  <- est.eqn.grad.tmp2( s.phd.1$solver.obj$par)
        #grad.numerical   <- grad(est.eqn.tmp2, s.phd.1$solver.obj$par)

        #cor(grad.numerical, grad.analytical)


        # angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a)
        # angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a)
        # angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a)
        # angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a)
        # angles(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a)

        # angles(phd.1$beta.hat[1:nvars,1], beta.a)
        # angles(sir.1$beta.hat[1:nvars,1], beta.a)
        # angles(s.phd.1$beta[1:nvars,1], beta.a)

         sim.subsp.angle.res.list[[n]][s,1] <- (1/5) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
          + 3 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

         sim.subsp.angle.res.list[[n]][s,2] <- (1/5) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                  + 3 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

         sim.subsp.angle.res.list[[n]][s,3] <- (1/5) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                  + 3 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

         sim.subsp.angle.res.list[[n]][s,7] <- (1/5) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                  + 3 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

         sim.subsp.angle.res.list[[n]][s,5] <- (1/5) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                  + 3 * angles(phd.3$beta.hat, beta.ab))

         sim.subsp.angle.res.list[[n]][s,4] <- (1/5) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                  + 3 * angles(sir.3$beta.hat, beta.ab))

         sim.subsp.angle.res.list[[n]][s,6] <- (1/5) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                  + 3 * angles(s.phd.3$beta, beta.ab))

        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor2 <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        semi.hier.phd.cor <- semi.hier.phd.cor / 5
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5
        s.phd.cor    <- s.phd.cor / 5

        hier.sir.cor
        hier.phd.cor
        semi.hier.phd.cor
        semi.hier.phd.cor2
        semi.hier.phd.cor.u
        phd.cor
        sir.cor
        s.phd.cor

        sim.direction.res.list[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list[[n]][s,2] <- hier.phd.cor
        sim.direction.res.list[[n]][s,3] <- semi.hier.phd.cor

        sim.direction.res.list[[n]][s,4] <- sir.cor
        sim.direction.res.list[[n]][s,5] <- phd.cor
        sim.direction.res.list[[n]][s,6] <- s.phd.cor
        sim.direction.res.list[[n]][s,7] <- semi.hier.phd.cor2

#         directions.sir      <- as.matrix(x %*% Re(hier.sdr$beta.hat))
#         directions.sir.test <- as.matrix(x.test %*% Re(hier.sdr$beta.hat))
#         directions.phd      <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
#         directions.phd.test <- as.matrix(x.test %*% Re(hier.sdr.phd$beta.hat))
#         dir.phd.non.hier      <- as.matrix(x %*% Re(sdr.phd$beta.hat))
#         dir.phd.non.hier.test <- as.matrix(x.test %*% Re(sdr.phd$beta.hat))
#         dir.sir.non.hier      <- as.matrix(x %*% Re(sdr.sir$beta.hat))
#         dir.sir.non.hier.test <- as.matrix(x.test %*% Re(sdr.sir$beta.hat))
#
#
#         directions.phd.1 <- x.list[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
#         directions.phd.2 <- x.list[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
#         directions.phd.3 <- x.list[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]
#
#         directions.phd.test.1 <- x.list.test[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
#         directions.phd.test.2 <- x.list.test[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
#         directions.phd.test.3 <- x.list.test[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]
#
#         #
#         directions.sir.1 <- x.list[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
#         directions.sir.2 <- x.list[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
#         directions.sir.3 <- x.list[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]
#
#         directions.sir.test.1 <- x.list.test[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
#         directions.sir.test.2 <- x.list.test[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
#         directions.sir.test.3 <- x.list.test[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]
#
#
#
#
#
#
#
#         h <- seq(1, 15, length.out = 15)
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir, y = y, alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.sir <- smooth.lf(x = directions.sir, y = y, alpha = best.h,
#                                                  xev = directions.sir.test))
#
#         ## hier SIR
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.1, y = y[1:nobs], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.sir.a <- smooth.lf(x = directions.sir.1, y = y[1:nobs], alpha = best.h,
#                                                    xev = directions.sir.test.1))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.sir.b <- smooth.lf(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
#                                                    xev = directions.sir.test.2))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.sir.ab <- smooth.lf(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
#                                                     xev = directions.sir.test.3))
#
#         hier.sir.preds <- c(locfit.hier.sir.a$y, locfit.hier.sir.b$y, locfit.hier.sir.ab$y)
#
#         ## hier PHD
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.1, y = y[1:nobs], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.phd.a <- smooth.lf(x = directions.phd.1, y = y[1:nobs], alpha = best.h,
#                                                    xev = directions.phd.test.1))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.phd.b <- smooth.lf(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
#                                                    xev = directions.phd.test.2))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.phd.ab <- smooth.lf(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
#                                                     xev = directions.phd.test.3))
#
#         hier.phd.preds <- c(locfit.hier.phd.a$y, locfit.hier.phd.b$y, locfit.hier.phd.ab$y)
#
#
#
#
#         ###
#
#         ## non hier sir
#
#         directions.sir.nh.1      <- x.list[[1]] %*% sir.1$beta.hat
#         directions.sir.nh.test.1 <- x.list.test[[1]] %*% sir.1$beta.hat
#
#         directions.sir.nh.2      <- x.list[[2]] %*% sir.2$beta.hat
#         directions.sir.nh.test.2 <- x.list.test[[2]] %*% sir.2$beta.hat
#
#         directions.sir.nh.3      <- x.list[[3]] %*% sir.3$beta.hat
#         directions.sir.nh.test.3 <- x.list.test[[3]] %*% sir.3$beta.hat
#
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.1, y = y[1:nobs], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.sir.a <- smooth.lf(x = directions.sir.nh.1, y = y[1:nobs], alpha = best.h,
#                                                     xev = directions.sir.nh.test.1))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.sir.b <- smooth.lf(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
#                                                     xev = directions.sir.nh.test.2))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.sir.ab <- smooth.lf(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
#                                                      xev = directions.sir.nh.test.3))
#
#         non.hier.sir.preds <- c(locfit.nhier.sir.a$y, locfit.nhier.sir.b$y, locfit.nhier.sir.ab$y)
#
#
#         ## non hier phd
#
#         directions.phd.nh.1      <- x.list[[1]] %*% phd.1$beta.hat
#         directions.phd.nh.test.1 <- x.list.test[[1]] %*% phd.1$beta.hat
#
#         directions.phd.nh.2      <- x.list[[2]] %*% phd.2$beta.hat
#         directions.phd.nh.test.2 <- x.list.test[[2]] %*% phd.2$beta.hat
#
#         directions.phd.nh.3      <- x.list[[3]] %*% phd.3$beta.hat
#         directions.phd.nh.test.3 <- x.list.test[[3]] %*% phd.3$beta.hat
#
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.1, y = y[1:nobs], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.phd.a <- smooth.lf(x = directions.phd.nh.1, y = y[1:nobs], alpha = best.h,
#                                                     xev = directions.phd.nh.test.1))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.phd.b <- smooth.lf(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
#                                                     xev = directions.phd.nh.test.2))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.nhier.phd.ab <- smooth.lf(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
#                                                      xev = directions.phd.nh.test.3))
#
#         non.hier.phd.preds <- c(locfit.nhier.phd.a$y, locfit.nhier.phd.b$y, locfit.nhier.phd.ab$y)
#
#         ##
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd, y = y, alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.hier.phd <- smooth.lf(x = directions.phd, y = y, alpha = best.h,
#                                                  xev = directions.phd.test))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.phd.non.hier, y = y, alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.non.hier.phd <- smooth.lf(x = dir.phd.non.hier, y = y, alpha = best.h,
#                                                      xev = dir.phd.non.hier.test))
#
#         system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.sir.non.hier, y = y, alpha = hv)[4]))
#         best.h     <- h[which.min(gcv.vals)]
#         system.time(locfit.non.hier.sir <- smooth.lf(x = dir.sir.non.hier, y = y, alpha = best.h,
#                                                      xev = dir.sir.non.hier.test))
#
#         sim.res.list[[n]][s,1] <- 1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#         sim.res.list[[n]][s,2] <- 1 - mean((y.test - hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#
#         sim.res.list[[n]][s,3] <- 1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#         sim.res.list[[n]][s,4] <- 1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#         sim.res.list[[n]][s,5] <- 1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#         sim.res.list[[n]][s,6] <- 1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#
#         sim.res.list[[n]][s,7] <- 1 - mean((y.test - non.hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#         sim.res.list[[n]][s,8] <- 1 - mean((y.test - non.hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)


        if (s %% 1 == 0) cat("sim:", s, "complete \n")
    }
}

res <- do.call(rbind, lapply(sim.res.list, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list, melt))
res.angle <- do.call(rbind, lapply(sim.subsp.angle.res.list, melt))




library(reshape2)
library(ggplot2)
df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1, 7, 2, 8, 3, 5, 4, 6)])

pdf(paste0(fig.path, "sim_semi_1a_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()

## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list[[1]])))
df.m.angle <- data.frame(res.angle, nobs = rep(nobs.vec, each = nsims * ncol(sim.subsp.angle.res.list[[1]])))

colnames(df.m.dir)[2:3] <- colnames(df.m.angle)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,4,2,5,6,3,7)])
df.m.angle$Method <- factor(df.m.angle$Method, levels = levels(df.m.angle$Method)[c(1,4,2,5,6,3,7)])


pdf(paste0(fig.path, "sim_semi_1a_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

dev.off()

pdf(paste0(fig.path, "sim_semi_1a_angles_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.angle, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("angle") + facet_wrap(~ nobs) + theme_bw()

dev.off()






############################
#
#   SIM 2
#
############################





# simulation parameters
nsims     <- 50
nobs.vec  <- c(250, 500, 1000, 2000)[1:4]
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 1

set.seed(123)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 7))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "semi hier phd", "sir separate", "phd separate", "semi phd separate", "semi hier phd 2")

sim.res.list2 <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list2 <- rep(list(sim.res.dir), length(nobs.vec))
sim.subsp.angle.res.list2 <- sim.direction.res.list2

library(MASS)
rho <- 0.5
cov.mat <- rho ^ abs(outer(1:(nvars - 20), 1:(nvars - 20), "-"))

snr2 <- numeric(length(nobs.vec) * nsims)
ct <- 0
for (n in 1:length(nobs.vec))
{
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {
        ct <- ct + 1
        x.list <- x.list.test <- vector(mode = "list", length = 3)
        for (g in 1:3)
        {
            x.tmp  <- mvrnorm(n = nobs, mu = rep(0, nvars - 20), Sigma = cov.mat)
            x.tmp2 <- apply(x.tmp[,1:10], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
            x.tmp3 <- apply(x.tmp[,11:20], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
            x.list[[g]] <- data.matrix(cbind(x.tmp, x.tmp2, x.tmp3))

            x.tmp  <- mvrnorm(n = nobs.test, mu = rep(0, nvars - 20), Sigma = cov.mat)
            x.tmp2 <- apply(x.tmp[,1:10], 2, function(col) rbinom(length(col), 1, 1 / (1 + exp(-col))))
            x.tmp3 <- apply(x.tmp[,11:20], 2, function(col) 0.1 * col ^ 2 + rnorm(length(col), sd = 1))
            x.list.test[[g]] <- data.matrix(cbind(x.tmp, x.tmp2, x.tmp3))
        }


        #x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        #x.list      <- replicate(3, list(matrix(c(rnorm(nobs * nvars/2), rbinom(nobs * nvars/2, 1, 0.5)), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(c(rnorm(nobs.test * nvars/2), rbinom(nobs.test * nvars/2, 1, 0.5)), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(runif(nvars / 2, min = -0.01, max = 0.01),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        beta.b <- matrix(c(runif(nvars / 2, min = -0.01, max = 0.01),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        eta.ab <- matrix(c(runif(nvars / 2, min = -0.01, max = 0.01),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)


#         y.true.a <- apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
#         y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
#         y.true.ab <- (apply( (x.list[[3]] %*% beta.a) ^ 2, 1, sum)) +
#             0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
#             0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))
#
#         y.true <- c(y.true.a, y.true.b, y.true.ab)
#         y <- y.true + rnorm(nobs, sd = sd.sim)


        mult.mat <- matrix(runif(ncol(beta.ab) * 1, max = 0.5, min = -0.5), ncol = 1)
        beta.ab.reduced <- beta.ab %*% mult.mat

        y.true.a  <- 1 * sin((apply(1 * exp(x.list[[1]] %*% beta.a), 1, sum)))
        y.true.b  <- ((1 * x.list[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- -apply(exp(0.2 * x.list[[3]] %*% beta.ab[,1,drop=FALSE]) , 1, sum) *
            (1 + 1 * (apply(sin( 0.2 * pi * x.list[[3]] %*% beta.ab[,2,drop=FALSE]), 1, sum) )) +
            0.2 * (apply((x.list[[3]] %*% beta.ab[,3,drop=FALSE])^2 , 1, sum))



        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = sd.sim)


#         y.true.a <- apply( exp(x.list.test[[1]] %*% beta.a), 1, sum)
#         y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
#         y.true.ab <- (apply( (x.list.test[[3]] %*% beta.a) ^ 2, 1, sum)) +
#             0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
#             0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))
#
#         y.true.test <- c(y.true.a, y.true.b, y.true.ab)
#         y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)


        y.true.a  <- 1 * sin((apply(1 * exp(x.list.test[[1]] %*% beta.a), 1, sum)))
        y.true.b  <- ((1 * x.list.test[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- -apply(exp(0.2 * x.list.test[[3]] %*% beta.ab[,1,drop=FALSE]) , 1, sum) *
            (1 + 1 * (apply(sin( 0.2 * pi * x.list.test[[3]] %*% beta.ab[,2,drop=FALSE]), 1, sum) )) +
            0.2 * (apply((x.list.test[[3]] %*% beta.ab[,3,drop=FALSE])^2 , 1, sum))


        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

        snr2[ct] <- var(y.true) / var(y - y.true)
        print(snr2[ct])

        hier.sdr      <- hier.sir(x.list, y,  d = 1, h = 30L)
        hier.sdr.phd  <- hier.phd(x.list, y,  d = 1)
        sdr.sir       <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd       <- phd(as.matrix(x), y, d = 1 * 3)

        semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = 1, h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)
        semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = 1, h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)


        sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)

        s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)

        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        sim.subsp.angle.res.list2[[n]][s,1] <- (1/5) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                       + 3 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list2[[n]][s,2] <- (1/5) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                       + 3 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list2[[n]][s,3] <- (1/5) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                       + 3 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list2[[n]][s,7] <- (1/5) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                                                       + 3 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list2[[n]][s,5] <- (1/5) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                       + 3 * angles(phd.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list2[[n]][s,4] <- (1/5) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                       + 3 * angles(sir.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list2[[n]][s,6] <- (1/5) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                       + 3 * angles(s.phd.3$beta, beta.ab))

        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor2 <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        semi.hier.phd.cor <- semi.hier.phd.cor / 5
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5
        s.phd.cor    <- s.phd.cor / 5

        hier.sir.cor
        hier.phd.cor
        semi.hier.phd.cor
        semi.hier.phd.cor2
        semi.hier.phd.cor.u
        phd.cor
        sir.cor
        s.phd.cor

        sim.direction.res.list2[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list2[[n]][s,2] <- hier.phd.cor
        sim.direction.res.list2[[n]][s,3] <- semi.hier.phd.cor

        sim.direction.res.list2[[n]][s,4] <- sir.cor
        sim.direction.res.list2[[n]][s,5] <- phd.cor
        sim.direction.res.list2[[n]][s,6] <- s.phd.cor
        sim.direction.res.list2[[n]][s,7] <- semi.hier.phd.cor2

        if (s %% 1 == 0) cat("sim:", s, "complete \n")
        print(colMeans(sim.subsp.angle.res.list2[[n]][1:s,,drop=FALSE]))
        print(colMeans(sim.direction.res.list2[[n]][1:s,,drop=FALSE]))
    }
}


library(reshape2)
library(ggplot2)

res <- do.call(rbind, lapply(sim.res.list2, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list2, melt))
res.angle <- do.call(rbind, lapply(sim.subsp.angle.res.list2, melt))



## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list2[[1]])))
df.m.angle <- data.frame(res.angle, nobs = rep(nobs.vec, each = nsims * ncol(sim.subsp.angle.res.list2[[1]])))

colnames(df.m.dir)[2:3] <- colnames(df.m.angle)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,4,2,5,6,3,7)])
df.m.angle$Method <- factor(df.m.angle$Method, levels = levels(df.m.angle$Method)[c(1,4,2,5,6,3,7)])


pdf(paste0(fig.path, "sim_semi_2a_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir[which(df.m.dir$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

pdf(paste0(fig.path, "sim_semi_2a_angles_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.angle[which(df.m.angle$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("angle") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()





### modeling plots



df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1, 7, 2, 8, 3, 5, 4, 6)])

pdf(paste0(fig.path, "sim_semi_1a_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()


###











############################
#
#   SIM 3
#
############################





# simulation parameters
nsims     <- 50
nobs.vec  <- c(250, 500, 1000, 2000)[1:3]
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 2

set.seed(42)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 8))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "semi hier phd", "sir separate", "phd separate", "semi phd separate", "semi hier phd 2", "semi hier phd W")

sim.res.list3 <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list3 <- rep(list(sim.res.dir), length(nobs.vec))
sim.subsp.angle.res.list3 <- sim.direction.res.list3

library(MASS)
rho <- 0.5
cov.mat <- rho ^ abs(outer(1:(nvars - 10), 1:(nvars - 10), "-"))

snr2 <- numeric(length(nobs.vec) * nsims)
ct <- 0
for (n in 1:length(nobs.vec))
{
    # set.seed(123)
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {
        ct <- ct + 1
        x.list <- x.list.test <- vector(mode = "list", length = 3)
        for (g in 1:3)
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


        #x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        #x.list      <- replicate(3, list(matrix(c(rnorm(nobs * nvars/2), rbinom(nobs * nvars/2, 1, 0.5)), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(c(rnorm(nobs.test * nvars/2), rbinom(nobs.test * nvars/2, 1, 0.5)), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        #beta.a <- matrix(c(runif(nvars / 2, min = -0.05, max = 0.05),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        #beta.b <- matrix(c(runif(nvars / 2, min = -0.05, max = 0.05),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        #eta.ab <- matrix(c(runif(nvars / 2, min = -0.05, max = 0.05),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)

        beta.a <- matrix(c(numeric(nvars / 2),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        beta.b <- matrix(c(numeric(nvars / 2),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        eta.ab <- matrix(c(numeric(nvars / 2),
                           runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        #mult.a <- mult.b <- diag(nvars)

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b) #, eta.ab)


        #         y.true.a <- apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
        #         y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        #         y.true.ab <- (apply( (x.list[[3]] %*% beta.a) ^ 2, 1, sum)) +
        #             0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
        #             0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))
        #
        #         y.true <- c(y.true.a, y.true.b, y.true.ab)
        #         y <- y.true + rnorm(nobs, sd = sd.sim)


        mult.mat <- matrix(runif(ncol(beta.ab) * 1, max = 0.5, min = -0.5), ncol = 1)
        beta.ab.reduced <- beta.ab %*% mult.mat

        y.true.a <- 0.5 * apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
        y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 3)
        y.true.ab <- #-0.5 * (apply( (x.list[[3]] %*% beta.ab[,3]) ^ 2, 1, sum)) -
            x.list[[3]] %*% beta.ab[,1] / (0.5 + (x.list[[3]] %*% beta.ab[,2] + 1.5) ^ 2) #+
        #0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,1]), 1, sum))



        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = sd.sim)


        #         y.true.a <- apply( exp(x.list.test[[1]] %*% beta.a), 1, sum)
        #         y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        #         y.true.ab <- (apply( (x.list.test[[3]] %*% beta.a) ^ 2, 1, sum)) +
        #             0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
        #             0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))
        #
        #         y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        #         y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)


        y.true.a <- 0.5 * apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum)
        y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 3)
        y.true.ab <- #-0.5 * (apply( (x.list.test[[3]] %*% beta.ab[,3]) ^ 2, 1, sum)) -
            x.list.test[[3]] %*% beta.ab[,1] / (0.5 + (x.list.test[[3]] %*% beta.ab[,2] + 1.5) ^ 2) #+
            #0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,1]), 1, sum))


        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

        snr2[ct] <- var(y.true) / var(y - y.true)
        print(snr2[ct])

        hier.sdr      <- hier.sir(x.list, y,  d = c(1,1,0), h = 30L)
        hier.sdr.phd  <- hier.phd(x.list, y,  d = c(1,1,0))
        sdr.sir       <- sir(as.matrix(x), y, d = 1 * 2, h = 30L)
        sdr.phd       <- phd(as.matrix(x), y, d = 1 * 2)

        semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = c(1,1,0), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)
        semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = c(1,1,0), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)


        sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2)

        s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)


        #phd.A <- phd(data.matrix(cbind(x.list[[1]], x.list[[3]])), c(y[1:nobs], y[(1 + 2*nobs):(3*nobs)]), d = 1)
        #phd.B <- phd(data.matrix(cbind(x.list[[2]], x.list[[3]])), c(y[1:nobs], y[(1 + 2*nobs):(3*nobs)]), d = 1)

        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        sim.subsp.angle.res.list3[[n]][s,1] <- (1/4) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 2 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list3[[n]][s,2] <- (1/4) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 2 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list3[[n]][s,3] <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 2 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        #sim.subsp.angle.res.list3[[n]][s,8] <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]),beta.b) )
        #                                                + 2 * angles(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list3[[n]][s,7] <- (1/4) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 2 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list3[[n]][s,5] <- (1/4) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 2 * angles(phd.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list3[[n]][s,4] <- (1/4) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 2 * angles(sir.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list3[[n]][s,6] <- (1/4) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                        + 2 * angles(s.phd.3$beta, beta.ab))

        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor   <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor2  <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        #hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        #hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        ##semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        #sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        #s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 4
        hier.phd.cor <- hier.phd.cor / 4
        semi.hier.phd.cor <- semi.hier.phd.cor / 4
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 4
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 4
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 4
        phd.cor      <- phd.cor / 4
        sir.cor      <- sir.cor / 4
        s.phd.cor    <- s.phd.cor / 4

        hier.sir.cor
        hier.phd.cor
        semi.hier.phd.cor
        #semi.hier.phd.cor.W
        semi.hier.phd.cor2
        semi.hier.phd.cor.u
        phd.cor
        sir.cor
        s.phd.cor

        sim.direction.res.list3[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list3[[n]][s,2] <- hier.phd.cor
        sim.direction.res.list3[[n]][s,3] <- semi.hier.phd.cor
        #sim.direction.res.list3[[n]][s,8] <- semi.hier.phd.cor.W

        sim.direction.res.list3[[n]][s,4] <- sir.cor
        sim.direction.res.list3[[n]][s,5] <- phd.cor
        sim.direction.res.list3[[n]][s,6] <- s.phd.cor
        sim.direction.res.list3[[n]][s,7] <- semi.hier.phd.cor2

        if (s %% 1 == 0) cat("sim:", s, "complete \n")
        print(colMeans(sim.subsp.angle.res.list3[[n]][1:s,,drop=FALSE]))
        print(colMeans(sim.direction.res.list3[[n]][1:s,,drop=FALSE]))
    }
}

mean(snr2)
median(snr2)

library(reshape2)
library(ggplot2)

res <- do.call(rbind, lapply(sim.res.list3, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list3, melt))
res.angle <- do.call(rbind, lapply(sim.subsp.angle.res.list3, melt))



## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list3[[1]])))
df.m.angle <- data.frame(res.angle, nobs = rep(nobs.vec, each = nsims * ncol(sim.subsp.angle.res.list3[[1]])))

colnames(df.m.dir)[2:3] <- colnames(df.m.angle)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,4,2,5,6,3,7)])
df.m.angle$Method <- factor(df.m.angle$Method, levels = levels(df.m.angle$Method)[c(1,4,2,5,6,3,7)])


pdf(paste0(fig.path, "sim_semi_3a_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir[which(df.m.dir$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

pdf(paste0(fig.path, "sim_semi_3a_angles_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.angle[which(df.m.angle$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("angle") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()
















############################
#
#   SIM 4
#
############################





# simulation parameters
nsims     <- 50
nobs.vec  <- c(250, 500, 1000, 2000)[1:3]
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 2

set.seed(42)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 8))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "semi hier phd", "sir separate", "phd separate", "semi phd separate", "semi hier phd 2", "semi hier phd W")

sim.res.list4 <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list4 <- rep(list(sim.res.dir), length(nobs.vec))
sim.subsp.angle.res.list4 <- sim.direction.res.list4

library(MASS)
rho <- 0.5
cov.mat <- rho ^ abs(outer(1:(nvars - 10), 1:(nvars - 10), "-"))

snr2 <- numeric(length(nobs.vec) * nsims)
ct <- 0
for (n in 1:length(nobs.vec))
{
    # set.seed(123)
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {
        ct <- ct + 1
        x.list <- x.list.test <- vector(mode = "list", length = 3)
        for (g in 1:3)
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


        x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        #x.list      <- replicate(3, list(matrix(c(rnorm(nobs * nvars/2), rbinom(nobs * nvars/2, 1, 0.5)), ncol = nvars)))
        #x.list.test <- replicate(3, list(matrix(c(rnorm(nobs.test * nvars/2), rbinom(nobs.test * nvars/2, 1, 0.5)), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)


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


        #beta.a <- matrix(c(numeric(nvars / 2),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        #beta.b <- matrix(c(numeric(nvars / 2),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)
        #eta.ab <- matrix(c(numeric(nvars / 2),
        #                   runif(nvars / 2, min = -0.5, max = 0.5)), ncol = 1)


        #         y.true.a <- apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
        #         y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        #         y.true.ab <- (apply( (x.list[[3]] %*% beta.a) ^ 2, 1, sum)) +
        #             0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
        #             0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))
        #
        #         y.true <- c(y.true.a, y.true.b, y.true.ab)
        #         y <- y.true + rnorm(nobs, sd = sd.sim)

        snr2[ct] <- var(y.true) / var(y - y.true)
        print(snr2[ct])

        hier.sdr      <- hier.sir(x.list, y,  d = c(1,1,1), h = 30L)
        hier.sdr.phd  <- hier.phd(x.list, y,  d = c(1,1,1))
        sdr.sir       <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd       <- phd(as.matrix(x), y, d = 1 * 3)

        semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = c(1,1,1), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)
        semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = c(1,1,1), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 250, maxk = 1200)


        sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)

        s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)
        s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 250, maxk = 450)


        #phd.A <- phd(data.matrix(cbind(x.list[[1]], x.list[[3]])), c(y[1:nobs], y[(1 + 2*nobs):(3*nobs)]), d = 1)
        #phd.B <- phd(data.matrix(cbind(x.list[[2]], x.list[[3]])), c(y[1:nobs], y[(1 + 2*nobs):(3*nobs)]), d = 1)

        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        sim.subsp.angle.res.list4[[n]][s,1] <- (1/5) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 3 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list4[[n]][s,2] <- (1/5) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 3 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list4[[n]][s,3] <- (1/5) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 3 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

        #sim.subsp.angle.res.list4[[n]][s,8] <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]),beta.b) )
        #                                                + 3 * angles(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list4[[n]][s,7] <- (1/5) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                                                        + 3 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

        sim.subsp.angle.res.list4[[n]][s,5] <- (1/5) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 3 * angles(phd.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list4[[n]][s,4] <- (1/5) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 3 * angles(sir.3$beta.hat, beta.ab))

        sim.subsp.angle.res.list4[[n]][s,6] <- (1/5) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                        + 3 * angles(s.phd.3$beta, beta.ab))

        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor   <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor2  <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        semi.hier.phd.cor <- semi.hier.phd.cor / 5
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 5
        semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
        semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5
        s.phd.cor    <- s.phd.cor / 5

        hier.sir.cor
        hier.phd.cor
        semi.hier.phd.cor
        #semi.hier.phd.cor.W
        semi.hier.phd.cor2
        semi.hier.phd.cor.u
        phd.cor
        sir.cor
        s.phd.cor

        sim.direction.res.list4[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list4[[n]][s,2] <- hier.phd.cor
        sim.direction.res.list4[[n]][s,3] <- semi.hier.phd.cor
        #sim.direction.res.list4[[n]][s,8] <- semi.hier.phd.cor.W

        sim.direction.res.list4[[n]][s,4] <- sir.cor
        sim.direction.res.list4[[n]][s,5] <- phd.cor
        sim.direction.res.list4[[n]][s,6] <- s.phd.cor
        sim.direction.res.list4[[n]][s,7] <- semi.hier.phd.cor2

        if (s %% 1 == 0) cat("sim:", s, "complete \n")
        print(colMeans(sim.subsp.angle.res.list4[[n]][1:s,,drop=FALSE]))
        print(colMeans(sim.direction.res.list4[[n]][1:s,,drop=FALSE]))
    }
}

mean(snr2)
median(snr2)

library(reshape2)
library(ggplot2)

res <- do.call(rbind, lapply(sim.res.list4, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list4, melt))
res.angle <- do.call(rbind, lapply(sim.subsp.angle.res.list4, melt))



## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list4[[1]])))
df.m.angle <- data.frame(res.angle, nobs = rep(nobs.vec, each = nsims * ncol(sim.subsp.angle.res.list4[[1]])))

colnames(df.m.dir)[2:3] <- colnames(df.m.angle)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,4,2,5,6,3,7)])
df.m.angle$Method <- factor(df.m.angle$Method, levels = levels(df.m.angle$Method)[c(1,4,2,5,6,3,7)])


pdf(paste0(fig.path, "sim_semi_4a_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir[which(df.m.dir$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

pdf(paste0(fig.path, "sim_semi_4a_angles_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.angle[which(df.m.angle$nobs < 2000),], aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("angle") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()







#################################


#################################























#################################


#################################





# simulation parameters
nsims     <- 100
nobs.vec  <- c(250, 500, 1000, 2000)
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 2

set.seed(123)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 4))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "sir separate", "phd separate")

sim.res.list <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list <- rep(list(sim.res.dir), length(nobs.vec))

for (n in 1:length(nobs.vec))
{
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {

        x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)

        y.true.a <- sin(apply(x.list[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list[[3]] %*% beta.ab[,3]), 1, sum))^2


        y.true.a <- apply(exp(x.list[[1]] %*% beta.a) , 1, sum)
        y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- (apply( (x.list[[3]] %*% beta.a) ^ 2, 1, sum)) +
            0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))

        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = sd.sim)


        y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list.test[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list.test[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))^2


        y.true.a <- apply( exp(x.list.test[[1]] %*% beta.a), 1, sum)
        y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- (apply( (x.list.test[[3]] %*% beta.a) ^ 2, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

        var(y.true) / var(y - y.true)

        hier.sdr     <- hier.sir(x.list, y,  d = 1, h = 30L)
        hier.sdr.phd <- hier.phd(x.list, y,  d = 1)
        sdr.sir      <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd      <- phd(as.matrix(x), y, d = 1 * 3)

        sir.1 <- sir(x.list[[1]], y[1:nobs], d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs], d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)



        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])

        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1],    beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5

        sim.direction.res.list[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list[[n]][s,2] <- hier.phd.cor

        sim.direction.res.list[[n]][s,3] <- sir.cor
        sim.direction.res.list[[n]][s,4] <- phd.cor

        directions.sir      <- as.matrix(x %*% Re(hier.sdr$beta.hat))
        directions.sir.test <- as.matrix(x.test %*% Re(hier.sdr$beta.hat))
        directions.phd      <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
        directions.phd.test <- as.matrix(x.test %*% Re(hier.sdr.phd$beta.hat))
        dir.phd.non.hier      <- as.matrix(x %*% Re(sdr.phd$beta.hat))
        dir.phd.non.hier.test <- as.matrix(x.test %*% Re(sdr.phd$beta.hat))
        dir.sir.non.hier      <- as.matrix(x %*% Re(sdr.sir$beta.hat))
        dir.sir.non.hier.test <- as.matrix(x.test %*% Re(sdr.sir$beta.hat))


        directions.phd.1 <- x.list[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.2 <- x.list[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.3 <- x.list[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.phd.test.1 <- x.list.test[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.test.2 <- x.list.test[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.test.3 <- x.list.test[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        #
        directions.sir.1 <- x.list[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.2 <- x.list[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.3 <- x.list[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.sir.test.1 <- x.list.test[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.test.2 <- x.list.test[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.test.3 <- x.list.test[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]







        h <- seq(1, 15, length.out = 15)
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir <- smooth.lf(x = directions.sir, y = y, alpha = best.h,
                                                 xev = directions.sir.test))

        ## hier SIR
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.a <- smooth.lf(x = directions.sir.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.sir.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.b <- smooth.lf(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.sir.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.ab <- smooth.lf(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.sir.test.3))

        hier.sir.preds <- c(locfit.hier.sir.a$y, locfit.hier.sir.b$y, locfit.hier.sir.ab$y)

        ## hier PHD

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.a <- smooth.lf(x = directions.phd.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.phd.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.b <- smooth.lf(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.phd.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.ab <- smooth.lf(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.phd.test.3))

        hier.phd.preds <- c(locfit.hier.phd.a$y, locfit.hier.phd.b$y, locfit.hier.phd.ab$y)




        ###

        ## non hier sir

        directions.sir.nh.1      <- x.list[[1]] %*% sir.1$beta.hat
        directions.sir.nh.test.1 <- x.list.test[[1]] %*% sir.1$beta.hat

        directions.sir.nh.2      <- x.list[[2]] %*% sir.2$beta.hat
        directions.sir.nh.test.2 <- x.list.test[[2]] %*% sir.2$beta.hat

        directions.sir.nh.3      <- x.list[[3]] %*% sir.3$beta.hat
        directions.sir.nh.test.3 <- x.list.test[[3]] %*% sir.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.a <- smooth.lf(x = directions.sir.nh.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.sir.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.b <- smooth.lf(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.sir.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.ab <- smooth.lf(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.sir.nh.test.3))

        non.hier.sir.preds <- c(locfit.nhier.sir.a$y, locfit.nhier.sir.b$y, locfit.nhier.sir.ab$y)


        ## non hier phd

        directions.phd.nh.1      <- x.list[[1]] %*% phd.1$beta.hat
        directions.phd.nh.test.1 <- x.list.test[[1]] %*% phd.1$beta.hat

        directions.phd.nh.2      <- x.list[[2]] %*% phd.2$beta.hat
        directions.phd.nh.test.2 <- x.list.test[[2]] %*% phd.2$beta.hat

        directions.phd.nh.3      <- x.list[[3]] %*% phd.3$beta.hat
        directions.phd.nh.test.3 <- x.list.test[[3]] %*% phd.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.a <- smooth.lf(x = directions.phd.nh.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.phd.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.b <- smooth.lf(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.phd.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.ab <- smooth.lf(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.phd.nh.test.3))

        non.hier.phd.preds <- c(locfit.nhier.phd.a$y, locfit.nhier.phd.b$y, locfit.nhier.phd.ab$y)

        ##

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd <- smooth.lf(x = directions.phd, y = y, alpha = best.h,
                                                 xev = directions.phd.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.phd.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.phd <- smooth.lf(x = dir.phd.non.hier, y = y, alpha = best.h,
                                                     xev = dir.phd.non.hier.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.sir.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.sir <- smooth.lf(x = dir.sir.non.hier, y = y, alpha = best.h,
                                                     xev = dir.sir.non.hier.test))

        sim.res.list[[n]][s,1] <- 1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,2] <- 1 - mean((y.test - hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,3] <- 1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,4] <- 1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,5] <- 1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,6] <- 1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,7] <- 1 - mean((y.test - non.hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,8] <- 1 - mean((y.test - non.hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)


        if (s %% 5 == 0) cat("sim:", s, "complete \n")
    }
}

res <- do.call(rbind, lapply(sim.res.list, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list, melt))




library(reshape2)
library(ggplot2)
df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1, 7, 2, 8, 3, 5, 4, 6)])

pdf(paste0(fig.path, "sim_1b_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()

## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list[[1]])))

colnames(df.m.dir)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,3,2,4)])

pdf(paste0(fig.path, "sim_1b_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

dev.off()
#############################



# simulation parameters
nsims     <- 100
nobs.vec  <- c(250, 500, 1000, 2000)
nobs.test <- 1e4
nvars     <- 100
sd.sim    <- 2

set.seed(123)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 4))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "sir separate", "phd separate")

sim.res.list <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list <- rep(list(sim.res.dir), length(nobs.vec))

for (n in 1:length(nobs.vec))
{
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {

        x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)

        y.true.a <- sin(apply(x.list[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list[[3]] %*% beta.ab[,3]), 1, sum))^2



        y.true.a <- 0.5 * sin(apply(exp(x.list[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 1 * (apply( (x.list[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) /
            (-1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2 + 1)) +
            0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))

        ################

        y.true.a <- 1 * sin(apply(exp(x.list[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 0.5 * (apply( exp( 0.2 * x.list[[3]] %*% beta.ab[,1]), 1, sum)) *
            (-1 * ( sin( 0.2 * pi * x.list[[3]] %*% beta.ab[,2]) ) ) +
            (apply(exp(0.2 * x.list[[3]] %*% beta.ab[,3]), 1, sum))

        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = sd.sim)

        var(y.true) / var(y - y.true)



        y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list.test[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list.test[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))^2


        y.true.a <- 0.5 * sin(apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 1 * (apply( (x.list.test[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) /
            (-1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2 + 1)) +
            0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))


        #########

        y.true.a <- 1 * sin(apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 0.5 * (apply( exp( 0.2 * x.list.test[[3]] %*% beta.ab[,1]), 1, sum)) *
            (-1 * ( sin( 0.2 * pi * x.list.test[[3]] %*% beta.ab[,2]) ) ) +
            (apply(exp(0.2 * x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

        var(y.true) / var(y - y.true)

        hier.sdr     <- hier.sir(x.list, y,  d = 1, h = 30L)
        hier.sdr.phd <- hier.phd(x.list, y,  d = 1)
        sdr.sir      <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd      <- phd(as.matrix(x), y, d = 1 * 3)

        sir.1 <- sir(x.list[[1]], y[1:nobs], d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs], d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)



        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])

        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1],    beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5

        sim.direction.res.list[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list[[n]][s,2] <- hier.phd.cor

        sim.direction.res.list[[n]][s,3] <- sir.cor
        sim.direction.res.list[[n]][s,4] <- phd.cor

        directions.sir      <- as.matrix(x %*% Re(hier.sdr$beta.hat))
        directions.sir.test <- as.matrix(x.test %*% Re(hier.sdr$beta.hat))
        directions.phd      <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
        directions.phd.test <- as.matrix(x.test %*% Re(hier.sdr.phd$beta.hat))
        dir.phd.non.hier      <- as.matrix(x %*% Re(sdr.phd$beta.hat))
        dir.phd.non.hier.test <- as.matrix(x.test %*% Re(sdr.phd$beta.hat))
        dir.sir.non.hier      <- as.matrix(x %*% Re(sdr.sir$beta.hat))
        dir.sir.non.hier.test <- as.matrix(x.test %*% Re(sdr.sir$beta.hat))


        directions.phd.1 <- x.list[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.2 <- x.list[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.3 <- x.list[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.phd.test.1 <- x.list.test[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.test.2 <- x.list.test[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.test.3 <- x.list.test[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        #
        directions.sir.1 <- x.list[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.2 <- x.list[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.3 <- x.list[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.sir.test.1 <- x.list.test[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.test.2 <- x.list.test[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.test.3 <- x.list.test[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]







        h <- seq(1, 15, length.out = 15)
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir <- smooth.lf(x = directions.sir, y = y, alpha = best.h,
                                                 xev = directions.sir.test))

        ## hier SIR
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.a <- smooth.lf(x = directions.sir.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.sir.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.b <- smooth.lf(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.sir.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.ab <- smooth.lf(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.sir.test.3))

        hier.sir.preds <- c(locfit.hier.sir.a$y, locfit.hier.sir.b$y, locfit.hier.sir.ab$y)

        ## hier PHD

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.a <- smooth.lf(x = directions.phd.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.phd.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.b <- smooth.lf(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.phd.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.ab <- smooth.lf(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.phd.test.3))

        hier.phd.preds <- c(locfit.hier.phd.a$y, locfit.hier.phd.b$y, locfit.hier.phd.ab$y)




        ###

        ## non hier sir

        directions.sir.nh.1      <- x.list[[1]] %*% sir.1$beta.hat
        directions.sir.nh.test.1 <- x.list.test[[1]] %*% sir.1$beta.hat

        directions.sir.nh.2      <- x.list[[2]] %*% sir.2$beta.hat
        directions.sir.nh.test.2 <- x.list.test[[2]] %*% sir.2$beta.hat

        directions.sir.nh.3      <- x.list[[3]] %*% sir.3$beta.hat
        directions.sir.nh.test.3 <- x.list.test[[3]] %*% sir.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.a <- smooth.lf(x = directions.sir.nh.1, y = y[1:nobs], alpha = best.h,
                                                    xev = directions.sir.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.b <- smooth.lf(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                    xev = directions.sir.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.ab <- smooth.lf(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                     xev = directions.sir.nh.test.3))

        non.hier.sir.preds <- c(locfit.nhier.sir.a$y, locfit.nhier.sir.b$y, locfit.nhier.sir.ab$y)


        ## non hier phd

        directions.phd.nh.1      <- x.list[[1]] %*% phd.1$beta.hat
        directions.phd.nh.test.1 <- x.list.test[[1]] %*% phd.1$beta.hat

        directions.phd.nh.2      <- x.list[[2]] %*% phd.2$beta.hat
        directions.phd.nh.test.2 <- x.list.test[[2]] %*% phd.2$beta.hat

        directions.phd.nh.3      <- x.list[[3]] %*% phd.3$beta.hat
        directions.phd.nh.test.3 <- x.list.test[[3]] %*% phd.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.a <- smooth.lf(x = directions.phd.nh.1, y = y[1:nobs], alpha = best.h,
                                                    xev = directions.phd.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.b <- smooth.lf(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                    xev = directions.phd.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.ab <- smooth.lf(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                     xev = directions.phd.nh.test.3))

        non.hier.phd.preds <- c(locfit.nhier.phd.a$y, locfit.nhier.phd.b$y, locfit.nhier.phd.ab$y)

        ##

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd <- smooth.lf(x = directions.phd, y = y, alpha = best.h,
                                                 xev = directions.phd.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.phd.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.phd <- smooth.lf(x = dir.phd.non.hier, y = y, alpha = best.h,
                                                     xev = dir.phd.non.hier.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.sir.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.sir <- smooth.lf(x = dir.sir.non.hier, y = y, alpha = best.h,
                                                     xev = dir.sir.non.hier.test))

        sim.res.list[[n]][s,1] <- 1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,2] <- 1 - mean((y.test - hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,3] <- 1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,4] <- 1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,5] <- 1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,6] <- 1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,7] <- 1 - mean((y.test - non.hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,8] <- 1 - mean((y.test - non.hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)


        if (s %% 5 == 0) cat("sim:", s, "complete \n")
    }
}

res <- do.call(rbind, lapply(sim.res.list, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list, melt))




library(reshape2)
library(ggplot2)
df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1, 7, 2, 8, 3, 5, 4, 6)])

pdf(paste0(fig.path, "sim_2a_nvars100_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()

## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list[[1]])))

colnames(df.m.dir)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,3,2,4)])

pdf(paste0(fig.path, "sim_2a_directions_nvars100_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

dev.off()

#############################






# simulation parameters
nsims     <- 100
nobs.vec  <- c(250, 500, 1000, 2000)
nobs.test <- 1e4
nvars     <- 50
sd.sim    <- 1

set.seed(123)


sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 4))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "sir separate", "phd separate")

sim.res.list <- rep(list(sim.res), length(nobs.vec))
sim.direction.res.list <- rep(list(sim.res.dir), length(nobs.vec))

for (n in 1:length(nobs.vec))
{
    nobs <- nobs.vec[n]

    for (s in 1:nsims)
    {

        x.list      <- replicate(3, list(matrix(rnorm(nobs * nvars), ncol = nvars)))
        x.list.test <- replicate(3, list(matrix(rnorm(nobs.test * nvars), ncol = nvars)))
        x <- bdiag(x.list)
        x.test <- as.matrix(bdiag(x.list.test))

        beta.a <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        beta.b <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)
        eta.ab <- matrix(c(rnorm(nvars / 2, sd = 0.025),
                           rnorm(nvars / 2, sd = 0.25)), ncol = 1)


        beta.ab <- cbind(beta.a, beta.b, eta.ab)

        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

        mult.a <- diag(runif(nvars, max = 2))
        mult.b <- diag(runif(nvars, max = 2))

        beta.ab <- cbind(mult.a %*% beta.a, mult.b %*% beta.b, eta.ab)

        y.true.a <- sin(apply(x.list[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list[[3]] %*% beta.ab[,3]), 1, sum))^2



        y.true.a <- 0.5 * sin(apply(exp(x.list[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 1 * (apply( (x.list[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) /
            (-1 * ((x.list[[3]] %*% beta.ab[,2]) ^ 2 + 1)) +
            0.1 * (apply(exp(x.list[[3]] %*% beta.ab[,3]), 1, sum))


        y.true.a <- 1 * sin(apply(exp(x.list[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 0.5 * (apply( exp( 0.2 * x.list[[3]] %*% beta.ab[,1]), 1, sum)) *
            (-1 * ( sin( 0.2 * pi * x.list[[3]] %*% beta.ab[,2]) ) ) +
            (apply(exp(0.2 * x.list[[3]] %*% beta.ab[,3]), 1, sum))


        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = sd.sim)

        var(y.true) / var(y - y.true)



        y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list.test[[1]] %*% beta.a, 1, sum))
        y.true.b <- + 0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list.test[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))^2


        y.true.a <- 0.5 * sin(apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 1 * (apply( (x.list.test[[3]] %*% beta.ab[,1]) ^ 2, 1, sum)) /
            (-1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2 + 1)) +
            0.1 * (apply(exp(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))


        y.true.a <- 1 * sin(apply(exp(x.list.test[[1]] %*% beta.a) , 1, sum))
        y.true.b <- + 1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2)
        y.true.ab <- 0.5 * (apply( exp( 0.2 * x.list.test[[3]] %*% beta.ab[,1]), 1, sum)) *
            (-1 * ( sin( 0.2 * pi * x.list.test[[3]] %*% beta.ab[,2]) ) ) +
            (apply(exp(0.2 * x.list.test[[3]] %*% beta.ab[,3]), 1, sum))

        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = sd.sim)

        var(y.true) / var(y - y.true)

        hier.sdr     <- hier.sir(x.list, y,  d = 1, h = 30L)
        hier.sdr.phd <- hier.phd(x.list, y,  d = 1)
        sdr.sir      <- sir(as.matrix(x), y, d = 1 * 3, h = 30L)
        sdr.phd      <- phd(as.matrix(x), y, d = 1 * 3)

        sir.1 <- sir(x.list[[1]], y[1:nobs], d = 1, h = 30L)
        sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1, h = 30L)
        sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

        phd.1 <- phd(x.list[[1]], y[1:nobs], d = 1)
        phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)], d = 1)
        phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)



        cor.directions <- function(a, b, x)
        {
            cov <- cov(x)
            R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
            R.sq
        }


        ## beta A
        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])

        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1],    beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5

        sim.direction.res.list[[n]][s,1] <- hier.sir.cor
        sim.direction.res.list[[n]][s,2] <- hier.phd.cor

        sim.direction.res.list[[n]][s,3] <- sir.cor
        sim.direction.res.list[[n]][s,4] <- phd.cor

        directions.sir      <- as.matrix(x %*% Re(hier.sdr$beta.hat))
        directions.sir.test <- as.matrix(x.test %*% Re(hier.sdr$beta.hat))
        directions.phd      <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
        directions.phd.test <- as.matrix(x.test %*% Re(hier.sdr.phd$beta.hat))
        dir.phd.non.hier      <- as.matrix(x %*% Re(sdr.phd$beta.hat))
        dir.phd.non.hier.test <- as.matrix(x.test %*% Re(sdr.phd$beta.hat))
        dir.sir.non.hier      <- as.matrix(x %*% Re(sdr.sir$beta.hat))
        dir.sir.non.hier.test <- as.matrix(x.test %*% Re(sdr.sir$beta.hat))


        directions.phd.1 <- x.list[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.2 <- x.list[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.3 <- x.list[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.phd.test.1 <- x.list.test[[1]] %*% Re(hier.sdr.phd$beta.hat)[1:nvars,1]
        directions.phd.test.2 <- x.list.test[[2]] %*% Re(hier.sdr.phd$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.phd.test.3 <- x.list.test[[3]] %*% Re(hier.sdr.phd$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        #
        directions.sir.1 <- x.list[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.2 <- x.list[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.3 <- x.list[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]

        directions.sir.test.1 <- x.list.test[[1]] %*% Re(hier.sdr$beta.hat)[1:nvars,1]
        directions.sir.test.2 <- x.list.test[[2]] %*% Re(hier.sdr$beta.hat)[(1 + nvars):(2 * nvars),2]
        directions.sir.test.3 <- x.list.test[[3]] %*% Re(hier.sdr$beta.hat)[(1 + 2 * nvars):(3 * nvars),]







        h <- seq(1, 15, length.out = 15)
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir <- smooth.lf(x = directions.sir, y = y, alpha = best.h,
                                                 xev = directions.sir.test))

        ## hier SIR
        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.a <- smooth.lf(x = directions.sir.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.sir.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.b <- smooth.lf(x = directions.sir.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.sir.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.sir.ab <- smooth.lf(x = directions.sir.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.sir.test.3))

        hier.sir.preds <- c(locfit.hier.sir.a$y, locfit.hier.sir.b$y, locfit.hier.sir.ab$y)

        ## hier PHD

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.a <- smooth.lf(x = directions.phd.1, y = y[1:nobs], alpha = best.h,
                                                   xev = directions.phd.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.b <- smooth.lf(x = directions.phd.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                   xev = directions.phd.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd.ab <- smooth.lf(x = directions.phd.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                    xev = directions.phd.test.3))

        hier.phd.preds <- c(locfit.hier.phd.a$y, locfit.hier.phd.b$y, locfit.hier.phd.ab$y)




        ###

        ## non hier sir

        directions.sir.nh.1      <- x.list[[1]] %*% sir.1$beta.hat
        directions.sir.nh.test.1 <- x.list.test[[1]] %*% sir.1$beta.hat

        directions.sir.nh.2      <- x.list[[2]] %*% sir.2$beta.hat
        directions.sir.nh.test.2 <- x.list.test[[2]] %*% sir.2$beta.hat

        directions.sir.nh.3      <- x.list[[3]] %*% sir.3$beta.hat
        directions.sir.nh.test.3 <- x.list.test[[3]] %*% sir.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.a <- smooth.lf(x = directions.sir.nh.1, y = y[1:nobs], alpha = best.h,
                                                    xev = directions.sir.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.b <- smooth.lf(x = directions.sir.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                    xev = directions.sir.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.sir.ab <- smooth.lf(x = directions.sir.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                     xev = directions.sir.nh.test.3))

        non.hier.sir.preds <- c(locfit.nhier.sir.a$y, locfit.nhier.sir.b$y, locfit.nhier.sir.ab$y)


        ## non hier phd

        directions.phd.nh.1      <- x.list[[1]] %*% phd.1$beta.hat
        directions.phd.nh.test.1 <- x.list.test[[1]] %*% phd.1$beta.hat

        directions.phd.nh.2      <- x.list[[2]] %*% phd.2$beta.hat
        directions.phd.nh.test.2 <- x.list.test[[2]] %*% phd.2$beta.hat

        directions.phd.nh.3      <- x.list[[3]] %*% phd.3$beta.hat
        directions.phd.nh.test.3 <- x.list.test[[3]] %*% phd.3$beta.hat


        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.1, y = y[1:nobs], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.a <- smooth.lf(x = directions.phd.nh.1, y = y[1:nobs], alpha = best.h,
                                                    xev = directions.phd.nh.test.1))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.b <- smooth.lf(x = directions.phd.nh.2, y = y[(nobs + 1):(2*nobs)], alpha = best.h,
                                                    xev = directions.phd.nh.test.2))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.nhier.phd.ab <- smooth.lf(x = directions.phd.nh.3, y = y[(2 * nobs + 1):(3*nobs)], alpha = best.h,
                                                     xev = directions.phd.nh.test.3))

        non.hier.phd.preds <- c(locfit.nhier.phd.a$y, locfit.nhier.phd.b$y, locfit.nhier.phd.ab$y)

        ##

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.phd, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.hier.phd <- smooth.lf(x = directions.phd, y = y, alpha = best.h,
                                                 xev = directions.phd.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.phd.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.phd <- smooth.lf(x = dir.phd.non.hier, y = y, alpha = best.h,
                                                     xev = dir.phd.non.hier.test))

        system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = dir.sir.non.hier, y = y, alpha = hv)[4]))
        best.h     <- h[which.min(gcv.vals)]
        system.time(locfit.non.hier.sir <- smooth.lf(x = dir.sir.non.hier, y = y, alpha = best.h,
                                                     xev = dir.sir.non.hier.test))

        sim.res.list[[n]][s,1] <- 1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,2] <- 1 - mean((y.test - hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,3] <- 1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,4] <- 1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,5] <- 1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,6] <- 1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim.res.list[[n]][s,7] <- 1 - mean((y.test - non.hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim.res.list[[n]][s,8] <- 1 - mean((y.test - non.hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)


        if (s %% 5 == 0) cat("sim:", s, "complete \n")
    }
}

res <- do.call(rbind, lapply(sim.res.list, melt))
res.dir <- do.call(rbind, lapply(sim.direction.res.list, melt))




library(reshape2)
library(ggplot2)
df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1, 7, 2, 8, 3, 5, 4, 6)])

pdf(paste0(fig.path, "sim_2b_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()

## directions R^2

df.m.dir <- data.frame(res.dir, nobs = rep(nobs.vec, each = nsims * ncol(sim.direction.res.list[[1]])))

colnames(df.m.dir)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir$Method <- factor(df.m.dir$Method, levels = levels(df.m.dir$Method)[c(1,3,2,4)])

pdf(paste0(fig.path, "sim_2b_directions_nvars50_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

dev.off()

########

plot(x = directions.phd[,1], y = y)
plot(x = directions.phd.test[,1], y = y.test)
plot(x = directions.phd.test[,2], y = y.test)
plot(x = directions.phd.test[,3], y = y.test)

plot(x = directions.sir.test[,1], y = y.test)
plot(x = directions.sir.test[,2], y = y.test)
plot(x = directions.sir.test[,3], y = y.test)



anova(lmsir3, lmsir)
anova(lmphd3, lmphd)
anova(lmphdnh3, lmphdnh)



beta.c     <- t(hier.sdr$beta.hat) %*% hier.sdr$cov %*% hier.sdr$beta.hat
beta.c.phd <- t(hier.sdr.phd$beta.hat) %*% hier.sdr.phd$cov %*% hier.sdr.phd$beta.hat

round(beta.c, 4)
round(beta.c.phd, 4)


beta.t <- matrix(rnorm(ncol(x) * 6), ncol = 6)
yy <- x %*% beta.t + matrix(rnorm(nrow(x) * 6), ncol = 6)

p <- ncol(x.list[[1]])
constraints <- list(t(rbind(cbind(diag(p), array(0, dim = c(p, p)), -diag(p)),
                            cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) ))    ), #
                    t(rbind(cbind(array(0, dim = c(p, p)), diag(p), -diag(p)),
                            cbind(diag(p), array(0, dim = c(p, p * 2)))  ) ),
                    t( rbind(cbind(diag(p), array(0, dim = c(p, p * 2))),
                             cbind(array(0, dim = c(p, p)), diag(p), array(0, dim = c(p, p)) )  ))   )


beta.ols <- solve(crossprod(x), crossprod(x, yy))


beta.constr <- array(0, dim = dim(beta.ols))
xtx <- crossprod(x)
for (i in 1:length(constraints))
{
    idx <- ((i - 1) * 2 + 1):(i * 2)
    AA <- constraints[[i]]

    f1 <- solve(xtx, AA)
    adj.fact <- f1 %*% solve(crossprod(AA, f1), crossprod(AA, beta.ols[, idx]))


    #adj.fact <- AA %*% solve(crossprod(AA), crossprod(AA, beta.ols[, idx]))

    beta.constr[, idx] <- as.matrix(beta.ols[, idx] - adj.fact)
}

round(cbind(beta.t, NA, beta.constr), 4)






###############################
#
# test dimension selection
#


cor.directions <- function(a, b, x)
{
    cov <- cov(x)
    R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
    R.sq
}


# simulation parameters
nsims     <- 50
nobs      <- c(250, 500, 1000, 2000)[2]
nobs.test <- 1e4
nvars     <- 10
sd.sim    <- 2

set.seed(12345)
d.estimates <- numeric(nsims)
for (s in 1:nsims)
{
    x <- matrix(rnorm(nobs * nvars), ncol = nvars)


    beta <- matrix(c(rep(0, nvars / 2),
                     rnorm(nvars / 2, sd = 0.5)), ncol = 1)
    beta <- data.matrix(cbind(beta, matrix(c(rep(0, nvars / 2),
                                             rnorm(nvars / 2, sd = 0.25)), ncol = 1)))

    y.true <- + 0.5 * ((x %*% beta[,1]) ^ 2) #- sin(pi * (x %*% beta[,2]) )  ^ 2# ^ 2

    y <- y.true + rnorm(nobs, sd = sd.sim)



    #phd <- semi.phd.dim.select.k(x, y, max.d = 3, k = 5, maxit = 15,
    #                             h = exp(seq(log(0.01), log(25), length.out = 15)), maxk = 450)
    phd <- phd.dim.select.k(x, y, max.d = 4, k = 5,
                            h = exp(seq(log(0.01), log(25), length.out = 15)), maxk = 450)
    phd$gcvs
    phd$d

    #s.phd <- semi.phd.dim.select(x, y, max.d = 3,
    #                             h = exp(seq(log(0.05), log(5), length.out = 15)),
    #                             maxit = 250, maxk = 450)

    d.estimates[s] <- phd$d
    print(phd$gcvs)
    print(d.estimates[1:s])
}

#phd <- phd(x, y, d = 10)


phd <- phd.dim.select(x, y, max.d = 3,
                      h = exp(seq(log(0.01), log(5), length.out = 15)), maxk = 450)
phd$gcvs
phd$d

sum(d.estimates == 0)
mean(d.estimates[d.estimates != 0] == 1)
mean(d.estimates[d.estimates != 0] == 2)
mean(d.estimates[d.estimates != 0] == 3)
mean(d.estimates[d.estimates != 0] == 4)



s.phd <- semi.phd.dim.select(x, y, max.d = 3,
                             h = exp(seq(log(0.05), log(5), length.out = 15)),
                             maxit = 250, maxk = 450)
str(s.phd)
s.phd$d
s.phd$bic

s.phd.b <- semi.phd.dim.select.boot(x, y, max.d = 3,
                                    h = exp(seq(log(0.05), log(5), length.out = 15)),
                                    maxit = 250, maxk = 450)
str(s.phd.b)
s.phd.b$d
s.phd.b$rsq

bic2 <- sapply(1:length(s.phd$models), function(i) s.phd$models[[i]]$solver.obj$value + 0.15 * log(nobs) * i )
bic2
which.min(bic2)
bic2 <- sapply(1:length(s.phd$models), function(i) s.phd$models[[i]]$solver.obj$value + 3 * (nobs^(-0.5)) * log(nobs) * i )
bic2

apply(s.phd$beta, 2, function(b) max(apply(beta, 2, function(be) cor.directions(b, be, x))))
apply(s.phd$beta.init, 2, function(b) max(apply(beta, 2, function(be) cor.directions(b, be, x))))


Proj <- function(b) b %*% solve(crossprod(b), t(b))


norm(Proj(s.phd$beta) - Proj(beta), type = "F")
norm(Proj(s.phd$beta.init) - Proj(beta), type = "F")

norm(Proj(matrix(rnorm(nvars), ncol=1)) - Proj(beta), type = "F")
