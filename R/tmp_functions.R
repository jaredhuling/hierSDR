


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
