



semi.phd.hier.newton <- function(x.list, y, d = rep(1L, 3L),
                                 maxit = 10L, h = NULL, B = NULL, vic = FALSE,
                                 weights = rep(1L, NROW(y)),
                                 opt.method = c("lbfgs.x", "bfgs", "lbfgs2",
                                                "bfgs.x",
                                                "lbfgs",
                                                "spg",
                                                "ucminf",
                                                "CG",
                                                "nlm",
                                                "nlminb",
                                                "newuoa"),
                                 init.method = c("random", "phd"),
                                 optimize.nn = FALSE,
                                 nn = NULL,
                                 separate.nn = FALSE,
                                 calc.mse = FALSE,
                                 verbose = TRUE, ...)
{
    p <- nvars <- ncol(x.list[[1]])

    opt.method  <- match.arg(opt.method)
    init.method <- match.arg(init.method)

    d <- as.vector(d)
    names(d) <- NULL

    nobs.vec  <- unlist(lapply(x.list, nrow))
    nvars.vec <- unlist(lapply(x.list, ncol))
    pp <- nvars.vec[1]

    nobs <- sum(nobs.vec)

    D <- sum(d)
    d.vic <- d + 1
    D.vic <- sum(d.vic)

    cum.d     <- c(0, cumsum(d))
    cum.d.vic <- c(0, cumsum(d.vic))

    cov <- lapply(x.list, cov)

    x <- as.matrix(bdiag(x.list))
    cov.b <- cov(as.matrix(x))
    eig.cov.b <- eigen(cov.b)
    eigs <- eig.cov.b$values
    eigs[eigs <= 0] <- 1e-5
    sqrt.inv.cov.b <- eig.cov.b$vectors %*% diag(1 / sqrt(eigs)) %*% t(eig.cov.b$vectors)
    x.tilde.b <- scale(x, scale = FALSE) %*% sqrt.inv.cov.b

    sqrt.inv.cov <- lapply(1:length(x.list), function(i) {
        eig.cov <- eigen(cov[[i]])
        eigvals <- eig.cov$values
        eigvals[eigvals <= 0] <- 1e-5
        eig.cov$vectors %*% diag(1 / sqrt(eigvals)) %*% t(eig.cov$vectors)
    })

    x.tilde <- lapply(1:length(x.list), function(i) {
        scale(x.list[[i]], scale = FALSE) %*% sqrt.inv.cov[[i]]# / sqrt(nobs.vec[i])
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

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)



    est.eqn <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        #beta.mat.list <- vector(mode = "list", length = 3L)

        #beta.mat.list[[1]] <- matrix(beta.vec[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
        #beta.mat.list[[2]] <- matrix(beta.vec[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
        #beta.mat.list[[3]] <- matrix(beta.vec[(cum.d[3] * p + 1):length(beta.vec)], ncol = cum.d[4] )

        cond.mat <- subpop.struct(2L)
        if (optimize.nn)
        {
            beta.mat.list <- vec2subpopMatsId(beta.vec[-1], p, d, cond.mat)
        } else
        {
            beta.mat.list <- vec2subpopMatsId(beta.vec, p, d, cond.mat)
        }

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }


        #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)


        Ey.given.xbeta <- numeric(nobs)


        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0),drop=FALSE]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         kern = "gauss",
            #                                         alpha = hv, deg = 3, ...)[4])
            #best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            if (optimize.nn)
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = 2, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val[s], best.h), deg = 2, ...)
            }



            #alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = 2, ...)

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod",
            #                          alpha = best.h, deg = 2, ...)


            # if (ncol(dir.cur) == 1)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur),...)
            # } else if (ncol(dir.cur) == 2)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1],dir.cur[,2]),...)
            # } else if (ncol(dir.cur) == 3)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1], dir.cur[,2], dir.cur[,3]),...)
            # } else if (ncol(dir.cur) == 4)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1], dir.cur[,2], dir.cur[,3], dir.cur[,4]),...)
            # }
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }



        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            wts.cur <- weights[strata.idx]
            norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        })))
        lhs / sum(weights)
    }





    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))

    npar <- sum(p * d - d ^ 2)
    V.hat <- crossprod(x.tilde.b, drop(scale(y, scale = FALSE)) * x.tilde.b) / nrow(x.tilde.b)
    for (c in 1:length(constraints))
    {
        #print(d)
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

            eig.c <- eigen(Proj.constr.list[[c]] %*% V.hat )
            eta.hat <- eig.c$vectors[,1:d[c], drop=FALSE]
            beta.list[[c]] <- Re(eta.hat)
        }
    }


    beta.list.init <- beta.list

    for (s in 1:length(beta.list.init))
    {
        beta.list.init[[s]] <- beta.list.init[[s]][(p * (s - 1) + 1):(p * s),,drop=FALSE]
    }

    beta.list.init <- lapply(beta.list.init, function(bb) {
        nc <- ncol(bb)
        bb[(nc + 1):nrow(bb),]
    })



    init <- unlist(beta.list.init)

    Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]

    #eig.V <- eigen(V.hat)
    beta.init  <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]


    d2 <- sapply(vec2subpopMatsId(init, p, d, cond.mat), ncol)

    est.eqn.unconstr <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        #beta.mat.list <- vector(mode = "list", length = 3L)

        #beta.mat.list[[1]] <- matrix(beta.vec[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
        #beta.mat.list[[2]] <- matrix(beta.vec[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
        #beta.mat.list[[3]] <- matrix(beta.vec[(cum.d[3] * p + 1):length(beta.vec)], ncol = cum.d[4] )

        cond.mat <- subpop.struct(2L)
        if (optimize.nn)
        {
            beta.mat.list <- vec2subpopMatsUnconstr(beta.vec[-1], p, d2, cond.mat)
        } else
        {
            beta.mat.list <- vec2subpopMatsUnconstr(beta.vec, p, d2, cond.mat)
        }

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }


        #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)


        Ey.given.xbeta <- numeric(nobs)


        for (s in 1:num.strata)
        {
            strata.idx <- which(strat.id == unique.strata[s])
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]
            #print(apply(dir.cur, 2, sd))
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0),drop=FALSE]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         kern = "gauss",
            #                                         alpha = hv, deg = 3, ...)[4])
            #best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            if (optimize.nn)
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = 2, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val[s], best.h), deg = 2, ...)
            }



            #alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = 2, ...)

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod",
            #                          alpha = best.h, deg = 2, ...)


            # if (ncol(dir.cur) == 1)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur),...)
            # } else if (ncol(dir.cur) == 2)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1],dir.cur[,2]),...)
            # } else if (ncol(dir.cur) == 3)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1], dir.cur[,2], dir.cur[,3]),...)
            # } else if (ncol(dir.cur) == 4)
            # {
            #     locfit.mod <- gam(y[strata.idx] ~ te(dir.cur[,1], dir.cur[,2], dir.cur[,3], dir.cur[,4]),...)
            # }
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }



        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            wts.cur <- weights[strata.idx]
            norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        })))
        lhs / sum(weights)
    }




    if (init.method == "random")
    {
        n.samples <- 100
        beta.list.phd <- beta.list.tmp <- beta.list
        best.value <- 1e99

        nn.vals <- c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95)
        for (tr in 1:n.samples)
        {
            par.cur <- runif(npar, min = min(init), max = max(init))

            values.cur <- numeric(length(nn.vals))
            for (i in 1:length(nn.vals) )
            {
                nh <- nn.vals[i]
                values.cur[i] <- est.eqn(par.cur, nn.val = nh)
            }

            value.cur <- min(values.cur)
            if (value.cur < best.value)
            {
                best.value <- value.cur
                best.par   <- par.cur
            }
        }
        init.rand <- best.par
    }






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
                                                 kern = "trwt", kt = "prod",
                                                 alpha = hv, deg = 3, ...)[4])
        best.h.vec[s]     <- h[which.min(gcv.vals)]
    }



    est.eqn.grad <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- grad(est.eqn, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
        if (verbose) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }

    est.eqn.grad.unconstr <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- grad(est.eqn.unconstr, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
        if (verbose) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }


    est.eqn.vic <- function(beta.vec, nn.val, v = 1)
    {
        #beta.mat.list <- vector(mode = "list", length = 3L)

        #beta.mat.list[[1]] <- matrix(beta.vec[(cum.d.vic[1] * p + 1):(cum.d.vic[2] * p)], ncol = d.vic[1])
        #beta.mat.list[[2]] <- matrix(beta.vec[(cum.d.vic[2] * p + 1):(cum.d.vic[3] * p)], ncol = d.vic[2])
        #beta.mat.list[[3]] <- matrix(beta.vec[(cum.d.vic[3] * p + 1):length(beta.vec)], ncol = cum.d.vic[4] )

        #t(t(matrix(beta.vec, ncol = d * length(constraints))) %*% sqrt.inv.cov)

        cond.mat <- subpop.struct(2L)
        beta.mat.list <- vec2subpopMatsIdVIC(beta.vec, p, d, cond.mat, v = v)

        Ey.given.xbeta <- numeric(nobs)

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }

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

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                     kern = "trwt", kt = "prod",
                                     alpha = c(nn.val[s], best.h), deg = 2, ...)

            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }


        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            wts.cur <- weights[strata.idx]
            norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        })))
        lhs / sum(weights)
    }



    beta.init.vec <- numeric(nvars.vec[1] * d[1] + nvars.vec[2] * d[2] + nvars.vec[3] * sum(d))

    beta.init.vec[1:(nvars.vec[1] * d[1])] <- as.vector(beta.init[1:nvars.vec[1], 1:d[1]])
    beta.init.vec[(nvars.vec[1] * d[1] + 1):
                      (nvars.vec[1] * d[1] + nvars.vec[2] * d[2])] <- as.vector(beta.init[(nvars.vec[1] + 1):
                                                                                              (nvars.vec[1] + nvars.vec[2]), (d[1] + 1):(d[1] + d[2]) ])
    beta.init.vec[(nvars.vec[1] * d[1] + nvars.vec[2] * d[2] + 1):
                      (length(beta.init.vec))] <- as.vector(beta.init[(nvars.vec[1] + nvars.vec[2] + 1):
                                                                          (nrow(beta.init)), ])

    # print("two eq:")
    # print(head(beta.init[1:nvars.vec[1], 1:d[1]]))
    # print(head(beta.init[(nvars.vec[1] + 1):
    #                          (nvars.vec[1] + nvars.vec[2]), (d[1] + 1):(d[1] + d[2]) ]))
    # print(head(beta.init[(nvars.vec[1] + nvars.vec[2] + 1):
    #                          (nrow(beta.init)), ]))

    # slver <-   optim(par     = beta.init.vec, # beta.init[(d+1):nrow(beta.init),],
    #                  fn      = est.eqn,
    #                  gr      = est.eqn.grad,
    #                  method  = "L-BFGS",
    #                  control = list(maxit = maxit, factr = 1e-10))





    # slver <- newton(par = unlist(beta.list.init), # subpopMats2vec(beta.list, d),
    #                 fn = est.eqn, maxit = maxit, tol = 1e-8,
    #                 verbose = verbose)



    # test which nn values minimize the most effectively
    if (is.null(nn))
    {
        tryval <- try.nn(nn.vals      = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95),
                         init         = init,
                         est.eqn      = est.eqn,
                         est.eqn.grad = est.eqn.grad,
                         opt.method   = "spg",
                         optimize.nn  = optimize.nn,
                         separate.nn  = separate.nn,
                         num.subpops  = num.strata,
                         maxit        = 20L,
                         verbose      = verbose)
        nn   <- tryval$nn
        init <- tryval$par

        if (init.method == "random")
        {
            tryval2 <- try.nn(nn.vals      = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95),
                              init         = init.rand,
                              est.eqn      = est.eqn,
                              est.eqn.grad = est.eqn.grad,
                              opt.method   = "spg",
                              optimize.nn  = optimize.nn,
                              separate.nn  = separate.nn,
                              num.subpops  = num.strata,
                              maxit        = 20L,
                              verbose      = verbose)
            nn2   <- tryval2$nn
            init2 <- tryval2$par

            if (tryval$value > tryval2$value)
            {
                init <- init2
                nn   <- nn2
            }
        }

        if (verbose) print(paste("best nn:", nn))
    } else
    {
        if (init.method == "random") init <- init.rand
    }

    slver <- opt.est.eqn(init         = init,
                         est.eqn      = est.eqn,
                         est.eqn.grad = est.eqn.grad,
                         opt.method   = opt.method,
                         nn           = nn,
                         optimize.nn  = optimize.nn,
                         maxit        = maxit,
                         verbose      = verbose)

    if (verbose)
    {
        print(slver)

        print("done optimizing :)")
    }


    beta <- beta.init

    cond.mat <- subpop.struct(2L)
    beta.mat.list <- vec2subpopMatsId(slver$par, p, d, cond.mat)

    if (verbose) print("Unconstr opt")
    slver.unconstr <- opt.est.eqn(init         = unlist(beta.mat.list),
                                  est.eqn      = est.eqn.unconstr,
                                  est.eqn.grad = est.eqn.grad.unconstr,
                                  opt.method   = opt.method,
                                  nn           = nn,
                                  optimize.nn  = optimize.nn,
                                  maxit        = 50,
                                  verbose      = verbose)

    #init2 <- Unconstrvec2Constrvec(slver.unconstr$par, p, d, cond.mat)
    #value.init2 <- est.eqn.vic(init2, nn.val = nn)




    #beta.mat.list <- vector(mode = "list", length = 3L)

    #beta.mat.list[[1]] <- matrix(slver$par[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
    #beta.mat.list[[2]] <- matrix(slver$par[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
    #beta.mat.list[[3]] <- matrix(slver$par[(cum.d[3] * p + 1):length(slver$par)], ncol = cum.d[4] )




    if (verbose) print("computing VIC's")

    vic.eqns <- c(est.eqn.vic(slver$par, nn.val = nn, v = -1),
                  est.eqn.vic(slver$par, nn.val = nn, v = -0.5),
                  est.eqn.vic(slver$par, nn.val = nn, v = 0.1),
                  est.eqn.vic(slver$par, nn.val = nn, v = 0.5),
                  est.eqn.vic(slver$par, nn.val = nn, v = 1))

    if (verbose) print("computed VIC's")

    vic.eqn <- mean(vic.eqns)

    vic  <- vic.eqn     + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    vic2 <- slver$value + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    vic3 <- slver.unconstr$value + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))

    model.list <- vector(mode = "list", length = 3)


    sse.vec <- mse.vec <- numeric(3)
    if (calc.mse)
    {
        for (m in 1:3)
        {
            strata.idx <- which(strat.id == unique.strata[m])
            best.h     <- best.h.vec[m]
            dir.cur <- x.tilde[[m]] %*% beta.mat.list[[m]]
            #model.list[[m]] <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                              kern = "gauss",
            #                              alpha = best.h, deg = 2, ...)

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            model.list[[m]] <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                          kern = "trwt", kt = "prod",
                                          alpha = c(nn, best.h), deg = degree, ...)

            fitted.vals <- fitted(model.list[[m]])

            sse.vec[m] <-  sum((y[strata.idx] - fitted.vals) ^ 2)
            mse.vec[m] <- mean((y[strata.idx] - fitted.vals) ^ 2)
        }
    }


    #beta.semi <- matrix(slver$par, ncol = D)
    #t(t(matrix(slver$par, ncol = d * length(constraints))) %*% sqrt.inv.cov)


    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.mat.list, beta.init = beta,
         solver.obj = slver,
         #beta.rand.init = t(t(beta.rand.init) %*% sqrt.inv.cov),
         cov = cov, sqrt.inv.cov = sqrt.inv.cov,
         nn  = nn,
         value = slver$value, value.init = est.eqn(init, nn.val = nn),
         value.unconstr = slver.unconstr$value,
         vic.est.eqn = vic.eqn, vic.eqns = vic.eqns,
         vic = vic, vic2 = vic2, vic3 = vic3,
         sse = sse.vec, mse = mse.vec)
}

hasAllConds <- function(checkvec, combvec)
{
    comb <- which(combvec != 0)
    all(checkvec[comb] != 0)
}


hier.sphd <- function(x, y, z, z.combinations, d,
                      maxit = 250L,
                      h = NULL,
                      weights = rep(1L, NROW(y)),
                      opt.method = c("lbfgs2", "lbfgs.x",
                                     "bfgs.x",
                                     "bfgs",
                                     "lbfgs",
                                     "spg",
                                     "ucminf",
                                     "CG",
                                     "nlm",
                                     "nlminb",
                                     "newuoa"),
                      init.method           = c("random", "phd"),
                      vic                   = TRUE,     # should the VIC be calculated?
                      vic.free.params       = FALSE,
                      grassmann             = TRUE,     # constrain parameters to the Grassmann manifold?
                      nn                    = NULL,     # nn fraction. Will be chosen automatically if nn = NULL
                      nn.try                = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95), # values to try if nn = NULL (more values takes longer)
                      optimize.nn           = FALSE,    # should nn be optimized? not recommended.
                      separate.nn           = FALSE,    # should each subpopulation have a separate nn? only used if nn = NULL
                      calc.mse              = FALSE,    # should the MSE of each subpop model be calculated after fitting?
                      constrain.none.subpop = FALSE,    # should we constrain S_{none} \subseteq S_{M} for all M?
                      verbose               = TRUE,     # should messages be printed out during optimization?
                      degree                = 2,        # degree of kernel
                      pooled                = FALSE,
                      ...)
{

    nobs <- NROW(x)

    if (nobs != NROW(z)) stop("number of observations in x doesn't match number of observations in z")

    z.combinations <- data.matrix(z.combinations)
    combinations   <- apply(z.combinations, 1, function(rr) paste(rr, collapse = ","))
    n.combinations <- length(combinations)

    if (nrow(z.combinations) != n.combinations) stop("duplicated row in z.combinations")

    z.combs   <- apply(z, 1, function(rr) paste(rr, collapse = ","))
    n.combs.z <- length(unique(z.combs))

    if (n.combinations != n.combs.z) stop("number of combinations of factors in z differs from that in z.combinations")

    x.list <- y.list <- weights.list <- strat.idx.list <- vector(mode = "list", length = n.combinations)

    for (c in 1:n.combinations)
    {
        strat.idx.list[[c]] <- which(z.combs == combinations[c])
        x.list[[c]]         <- x[strat.idx.list[[c]],]
        y.list[[c]]         <- y[strat.idx.list[[c]]]
        weights.list[[c]]   <- weights[strat.idx.list[[c]]]
    }


    p <- nvars <- ncol(x)

    opt.method  <- match.arg(opt.method)
    init.method <- match.arg(init.method)

    d <- as.vector(d)
    names(d) <- NULL

    if (length(d) != nrow(z.combinations)) stop("number of subpopulations implied by 'z.combinations' does not match that of 'd'")
    if (length(d) != length(x.list) ) stop("number of subpopulations implied by 'x.list' does not match that of 'd'")

    nobs.vec  <- unlist(lapply(x.list, nrow))
    nvars.vec <- unlist(lapply(x.list, ncol))
    pp        <- nvars.vec[1]

    nobs <- sum(nobs.vec)

    if (nobs != NROW(y)) stop("unequal number of observations in x.list and y")

    D <- sum(d)
    d.vic <- d + 1
    D.vic <- sum(d.vic)

    cum.d     <- c(0, cumsum(d))
    cum.d.vic <- c(0, cumsum(d.vic))

    cov <- lapply(x.list, cov)

    x.big <- as.matrix(bdiag(x.list))
    cov.b <- cov(as.matrix(x.big))
    eig.cov.b <- eigen(cov.b)
    eigs <- eig.cov.b$values
    eigs[eigs <= 0] <- 1e-5
    sqrt.inv.cov.b <- eig.cov.b$vectors %*% diag(1 / sqrt(eigs)) %*% t(eig.cov.b$vectors)
    x.tilde.b <- scale(x.big, scale = FALSE) %*% sqrt.inv.cov.b

    sqrt.inv.cov <- lapply(1:length(x.list), function(i) {
        eig.cov <- eigen(cov[[i]])
        eigvals <- eig.cov$values
        eigvals[eigvals <= 0] <- 1e-5
        eig.cov$vectors %*% diag(1 / sqrt(eigvals)) %*% t(eig.cov$vectors)
    })

    x.tilde <- lapply(1:length(x.list), function(i) {
        scale(x.list[[i]], scale = FALSE) %*% sqrt.inv.cov[[i]]# / sqrt(nobs.vec[i])
    })

    x.tilde.tall <- do.call(rbind, x.tilde)


    # construct linear constraint matrices
    constraints <- vector(mode = "list", length = nrow(z.combinations))

    names(constraints) <- combinations


    for (cr in 1:nrow(z.combinations))
    {
        cur.comb <- combinations[cr]

        other.idx <- (1:nrow(z.combinations))[-cr]

        first.mat <- TRUE
        for (i in other.idx)
        {
            if (hasAllConds(z.combinations[i,], z.combinations[cr,]) & !(!constrain.none.subpop & all(z.combinations[cr,] == 0)))
            { # if the ith subpopulation is nested within the one indexed by 'cr'
                # then we apply equality constraints (but only if it's not the none subpop
                # (unless we want to constrain the none subpop))

                constr.idx <- rep(0, nrow(z.combinations))
                constr.idx[c(cr, i)] <- c(1, -1)
                constr.mat.eq  <- do.call(rbind, lapply(constr.idx, function(ii) diag(rep(ii, p))))


                if (first.mat)
                {
                    first.mat  <- FALSE
                    constr.mat <- constr.mat.eq # cbind(constr.mat.eq, constr.mat.zero)
                } else
                {
                    constr.mat <- cbind(constr.mat, constr.mat.eq) # constr.mat.zero
                }
            } else
            {
                # else we apply constraints to be zero

                constr.idx.zero <- rep(0, nrow(z.combinations))
                constr.idx.zero[i] <- 1
                constr.mat.zero <- do.call(rbind, lapply(constr.idx.zero, function(ii) diag(rep(ii, p))))

                if (first.mat)
                {
                    first.mat  <- FALSE
                    constr.mat <- constr.mat.zero # cbind(constr.mat.eq, constr.mat.zero)
                } else
                {
                    constr.mat <- cbind(constr.mat, constr.mat.zero) # constr.mat.zero
                }
            }
        }

        constraints[[cr]] <- constr.mat

    }

    strat.id <- unlist(lapply(1:length(x.list), function(id) rep(id, nrow(x.list[[id]]))))

    V.hat    <- crossprod(x.tilde.b, drop(scale(y, scale = FALSE)) * x.tilde.b) / nrow(x.tilde.b)

    #V.hat.us    <- crossprod(x.big, drop(scale(y, scale = FALSE)) * x.big) / nrow(x.big)


    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))

    #qrlist <- lapply(1:length(x.list), function(i) {
    #    qr(scale(x.list[[i]], scale = FALSE))
    #})

    #qr.b <- qr(scale(x.big, scale = FALSE))

    #qrz <- qr(scale(x, center=TRUE, scale=FALSE))


    for (c in 1:length(constraints))
    {
        #print(d)
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

            #eig.c          <- eigen(Proj.constr.list[[c]] %*% V.hat )
            #eta.hat        <- eig.c$vectors[,1:d[c], drop=FALSE]

            #real.eta <- Re(eta.hat)


            D          <- eigen(Proj.constr.list[[c]] %*% V.hat )

            or <- rev(order(abs(D$values)))
            evalues <- D$values[or]
            raw.evectors <- Re(D$vectors[,or])

            qr_R <- function(object)
            {
                qr.R(object)[1:object$rank,1:object$rank]
            }

            # evectors <- backsolve(sqrt(nrow(x.big))*qr_R(qr.b),raw.evectors)
            # evectors <- if (is.matrix(evectors)) evectors else matrix(evectors,ncol=1)
            # evectors <- apply(evectors,2,function(x) x/sqrt(sum(x^2)))

            #real.eta <- Re(raw.evectors[,1:d[c], drop=FALSE])
            real.eta <- Re(raw.evectors[,1:d[c], drop=FALSE])

            whichzero <- apply(real.eta, 2, function(cl) all(abs(cl) <= 1e-12))

            #print((real.eta[,!whichzero,drop=FALSE]))
            #real.eta[,!whichzero] <- grassmannify(real.eta[,!whichzero,drop=FALSE])$beta

            beta.list[[c]] <- real.eta
        }
    }


    Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]

    #eig.V <- eigen(V.hat)
    beta.init     <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)

    beta.component.init <- beta.list
    ct <- 0

    for (s in 1:length(beta.list))
    {
        if (d[s] > 0)
        {
            beta.component.init[[s]] <- beta.list[[s]][(p * (s - 1) + 1):(p * s),,drop = FALSE]
            #beta.component.init[[s]] <- beta.list[[s]][(p * (s - 1) + 1):(p * s),,drop = FALSE]
        }
    }

    beta.component.init <- lapply(beta.component.init, function(bb) {
        if (!is.null(bb))
        {
            nc <- ncol(bb)

            whichzero <- apply(bb, 2, function(cl) all(abs(cl) <= 1e-12))

            #print((bb[,!whichzero,drop=FALSE]))
            #bb[,!whichzero] <- grassmannify(bb[,!whichzero,drop=FALSE])$beta


            bb <- grassmannify(bb)$beta
            bb[(nc + 1):nrow(bb),]
        } else
        {
            NULL
        }
    })

    init <- unlist(beta.component.init)

    txpy.list <- vector(mode = "list", length = n.combinations)
    for (s in 1:n.combinations)
    {
        txpy.list[[s]] <- vector(mode = "list", length = NROW(x.tilde[[s]]))
        for (j in 1:NROW(x.tilde[[s]])) txpy.list[[s]][[j]] <- tcrossprod(x.tilde[[s]][j,])
    }


    ncuts <- 3

    est.eqn <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        if (optimize.nn)
        {
            beta.mat.list  <- vec2subpopMatsId(beta.vec[-1], p, d, z.combinations)
        } else
        {
            beta.mat.list  <- vec2subpopMatsId(beta.vec, p, d, z.combinations)
        }

        Ey.given.xbeta <- numeric(nobs)

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }


        resid.Exxt.given.x.beta <- vector(mode = "list", length = nobs)


        Ey.given.xbeta.list <- vector(mode = "list", length = n.combinations)

        for (s in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[s]]
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]

            # remove any directions with no variation.
            # probably not needed, but just in case
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0), drop = FALSE]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            if (optimize.nn)
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = degree, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val[s], best.h), deg = degree, ...)
            }

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod",
            #                          alpha = best.h, deg = 2, ...)


            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
            Ey.given.xbeta.list[[s]]   <- Ey.given.xbeta[strata.idx]
            #nws <- nwcov(dir.cur, x = x.tilde[[s]], txpy.list = txpy.list[[s]], h = 0.25, max.points = 500, ncuts = ncuts)
            #resid.Exxt.given.x.beta[strata.idx] <- nws$resid
        }



        if (pooled)
        {
            resid      <- drop(unlist(y.list)) - drop(unlist(Ey.given.xbeta.list))
            return(norm(crossprod(x.tilde.tall, (weights * resid) * x.tilde.tall), type = "F") ^ 2 / sum(weights))
        } else
        {
            lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
                strata.idx <- strat.idx.list[[i]]
                resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
                wts.cur    <- weights.list[[i]]

                # lhss <- 0
                # for (j in 1:length(strata.idx))
                # {
                #     lhss <- lhss + sum((resid.Exxt.given.x.beta[[strata.idx[j]]] * resid[j]) ^ 2)
                # }
                # lhss
                norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
            })))
            return(lhs / sum(weights))
        }
    }


    # not used because it really doesn't work well
    est.eqn.by.component <- function(component.vec, nn.val, s.idx, beta.vec, optimize.nn = FALSE)
    {
        optimize.nn  <- FALSE
        cumul.params <- c(0, cumsum(d * p - d ^ 2))

        if (d[s.idx] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            beta.vec[vec.idx] <- component.vec
        }

        if (optimize.nn)
        {
            beta.mat.list  <- vec2subpopMatsId(beta.vec[-1], p, d, z.combinations)
        } else
        {
            beta.mat.list  <- vec2subpopMatsId(beta.vec, p, d, z.combinations)
        }

        Ey.given.xbeta <- numeric(nobs)

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }


        Ey.given.xbeta.list <- vector(mode = "list", length = n.combinations)

        for (s in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[s]]
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]

            # remove any directions with no variation.
            # probably not needed, but just in case
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0), drop = FALSE]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            if (optimize.nn)
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = degree, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val[s], best.h), deg = degree, ...)
            }

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod",
            #                          alpha = best.h, deg = 2, ...)

            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
            Ey.given.xbeta.list[[s]]   <- Ey.given.xbeta[strata.idx]
        }



        if (pooled)
        {
            resid      <- drop(y) - drop(unlist(Ey.given.xbeta.list))
            return(norm(crossprod(x.tilde.tall, (weights * resid) * x.tilde.tall), type = "F") ^ 2 / sum(weights))
        } else
        {
            lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
                strata.idx <- strat.idx.list[[i]]
                resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
                wts.cur    <- weights.list[[i]]

                # lhss <- 0
                # for (j in 1:length(strata.idx))
                # {
                #     lhss <- lhss + sum((resid.Exxt.given.x.beta[[strata.idx[j]]] * resid[j]) ^ 2)
                # }
                # lhss
                norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
            })))
            return(lhs / sum(weights))
        }
    }

    est.eqn.grad <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- numDeriv::grad(est.eqn, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
        if (verbose > 1) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }

    est.eqn.by.component.grad <- function(component.vec, nn.val, s.idx, beta.vec, optimize.nn = FALSE)
    {
        grad.full     <- numDeriv::grad(est.eqn.by.component, component.vec, method = "simple", nn.val = nn.val,
                              s.idx = s.idx, beta.vec = beta.vec, optimize.nn = optimize.nn)
        if (verbose > 1) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }


    est.eqn.vic <- function(beta.vec, nn.val, v = 1)
    {

        beta.mat.list <- vec2subpopMatsIdVIC(beta.vec, p, d, z.combinations, v = v)

        Ey.given.xbeta <- numeric(nobs)

        if (length(nn.val) == 1)
        {
            nn.val <- rep(nn.val, length(beta.mat.list))
        }

        Ey.given.xbeta.list <- vector(mode = "list", length = n.combinations)

        for (s in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[s]]
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]

            # remove a direction if it has no variation.
            # this should never happen, but just in case
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0)]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                     kern = "trwt", kt = "prod",
                                     alpha = c(nn.val[s], best.h), deg = degree, ...)



            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
            Ey.given.xbeta.list[[s]]   <- Ey.given.xbeta[strata.idx]
        }


        # lhs   <- sum(unlist(lapply(1:n.combinations, function(i) {
        #     strata.idx <- strat.idx.list[[i]]
        #     resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
        #     wts.cur    <- weights.list[[i]]
        #     norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        # })))
        # lhs / sum(weights)

        if (pooled)
        {
            resid      <- drop(y) - drop(unlist(Ey.given.xbeta.list))
            return(norm(crossprod(x.tilde.tall, (weights * resid) * x.tilde.tall), type = "F") ^ 2 / sum(weights))
        } else
        {
            lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
                strata.idx <- strat.idx.list[[i]]
                resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
                wts.cur    <- weights.list[[i]]

                norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
            })))
            return(lhs / sum(weights))
        }
    }



    #######    model fitting to determine bandwidth     ########

    best.h.vec <- numeric(num.strata)

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    beta.init.cov <- t(t(beta.init) %*% sqrt.inv.cov.b)
    # for (s in 1:n.combinations)
    # {
    #     strata.idx <- strat.idx.list[[s]]
    #     dir.cur    <- x.tilde.b[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
    #         beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]
    #
    #     # remove irrelevant dimensions
    #     dir.cur    <- dir.cur[, which(apply(dir.cur, 2, sd) != 0)]
    #
    #     gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
    #                                              y = y[strata.idx],
    #                                              kern = "trwt", kt = "prod",
    #                                              alpha = hv, deg = 3, ...)[4])
    #     best.h.vec[s]     <- h[which.min(gcv.vals)]
    # }



    d2   <- sapply(vec2subpopMatsId(init, p, d, z.combinations), ncol)
    npar <- sum(p * d - d ^ 2)

    value.init <- est.eqn(init, nn.val = nn)

    if (init.method == "random")
    {
        n.samples <- 100
        beta.list.phd <- beta.list.tmp <- beta.list
        best.value <- 1e99

        nn.vals <- c(0.05, 0.15, 0.25, 0.5, 0.75, 0.9, 0.95)
        for (tr in 1:n.samples)
        {
            par.cur <- runif(npar, min = min(init), max = max(init))

            values.cur <- numeric(length(nn.vals))
            for (i in 1:length(nn.vals) )
            {
                nh <- nn.vals[i]
                values.cur[i] <- est.eqn(par.cur, nn.val = nh)
            }

            value.cur <- min(values.cur)
            if (value.cur < best.value)
            {
                best.value <- value.cur
                best.par   <- par.cur
            }
        }
        init.rand <- best.par
    }





    #init <- rnorm(length(init))



    # test which nn values minimize the most effectively
    if (is.null(nn))
    {
        tryval <- try.nn(nn.vals      = nn.try,
                         init         = init,
                         est.eqn      = est.eqn,
                         est.eqn.grad = est.eqn.grad,
                         opt.method   = "spg",
                         optimize.nn  = optimize.nn,
                         separate.nn  = separate.nn,
                         num.subpops  = n.combinations,
                         maxit        = 20L,
                         verbose      = verbose)
        nn   <- tryval$nn
        init <- tryval$par


        if (init.method == "random")
        {
            tryval2 <- try.nn(nn.vals      = nn.try,
                              init         = init.rand,
                              est.eqn      = est.eqn,
                              est.eqn.grad = est.eqn.grad,
                              opt.method   = "spg",
                              optimize.nn  = optimize.nn,
                              separate.nn  = separate.nn,
                              num.subpops  = n.combinations,
                              maxit        = 20L,
                              verbose      = verbose)
            nn2   <- tryval2$nn
            init2 <- tryval2$par

            if (tryval$value > tryval2$value)
            {
                init <- init2
                nn   <- nn2
            }
        }

        if (verbose) print(paste("best nn:", nn))
    } else
    {
        if (init.method == "random") init <- init.rand
    }


    slver <- opt.est.eqn(init         = init,
                         est.eqn      = est.eqn,
                         est.eqn.grad = est.eqn.grad,
                         opt.method   = opt.method,
                         nn           = nn,
                         optimize.nn  = optimize.nn,
                         maxit        = maxit,
                         verbose      = verbose)


    beta.mat.list <- vec2subpopMatsId(slver$par, p, d, z.combinations)

    if (vic)
    {
        vic.eqns <- c(est.eqn.vic(slver$par, nn.val = nn, v = -1),
                      est.eqn.vic(slver$par, nn.val = nn, v = -0.5),
                      est.eqn.vic(slver$par, nn.val = nn, v = 0.1),
                      est.eqn.vic(slver$par, nn.val = nn, v = 0.5),
                      est.eqn.vic(slver$par, nn.val = nn, v = 1))

        vic.eqn <- mean(vic.eqns)
        if (vic.free.params)
        {
            vic1    <- vic.eqn       + log(sum(nobs.vec)) * length(slver$par)
            vic3    <- min(vic.eqns) + log(sum(nobs.vec)) * length(slver$par)
        } else
        {
            vic1    <- vic.eqn       + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
            vic3    <- min(vic.eqns) + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
        }
    } else
    {
        vic.eqns <- vic.eqn <- vic3 <- NULL
    }

    if (vic.free.params)
    {
        vic2 <- slver$value   + log(sum(nobs.vec)) * length(slver$par)
    } else
    {
        vic2 <- slver$value   + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    }



    model.list <- vector(mode = "list", length = 3)

    directions.list <- vector(mode = "list", length = n.combinations)

    for (m in 1:n.combinations)
    {
        directions.list[[m]] <- x.tilde[[m]] %*% beta.mat.list[[m]]
    }

    sse.vec <- mse.vec <- numeric(n.combinations)
    if (calc.mse)
    {
        for (m in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[m]]
            #best.h     <- best.h.vec[m]
            dir.cur    <- directions.list[[m]]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            model.list[[m]] <- locfit.raw(x = dir.cur, y = y.list[[m]],
                                          kern = "trwt", kt = "prod",
                                          alpha = c(nn, best.h), deg = 2, ...)

            fitted.vals <- fitted(model.list[[m]])

            sse.vec[m] <-  sum((y[strata.idx] - fitted.vals) ^ 2)
            mse.vec[m] <- mean((y[strata.idx] - fitted.vals) ^ 2)
        }
    }

    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    #beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)

    names(beta.mat.list) <- names(directions.list) <- combinations


    ret <- list(beta           = beta.mat.list,
                beta.init      = beta.init,
                directions     = directions.list,
                y.list         = y.list,
                z.combinations = z.combinations,
                cov            = cov,
                sqrt.inv.cov   = sqrt.inv.cov,
                solver.obj     = slver,
                value          = slver$value,
                value.init     = value.init,
                vic.est.eqn    = vic.eqn,
                vic.eqns       = vic.eqns,
                vic = vic1, vic2 = vic2, vic3 = vic3,
                sse = sse.vec, mse = mse.vec)
    class(ret) <- "hier_sdr_fit"
    ret
}

plot.hier_sdr_fit <- function(x)
{
    n.subpops    <- length(x$beta)
    subpop.names <- names(x$beta)
    dimensions   <- sapply(x$beta, ncol)

    maxd <- max(dimensions)
    par(mfrow = c(n.subpops, maxd))

    rbPal <- colorRampPalette(c('red','blue'))

    for (s in 1:n.subpops)
    {
        for (d in 1:maxd)
        {
            if (d <= dimensions[s])
            {
                col.use <- "#000000BF"
                if (dimensions[s] == 2)
                {
                    nc <- 12
                    if (d == 1)
                    {
                        col.use <- rbPal(nc)[as.numeric(cut(x$directions[[s]][,2],
                                                            breaks = nc))]
                    } else
                    {
                        col.use <- rbPal(nc)[as.numeric(cut(x$directions[[s]][,1],
                                                            breaks = nc))]
                    }

                }
                plot(x = x$directions[[s]][,d],
                     y = x$y.list[[s]], xlab = paste0("x * beta[", s, ",", d, "]"),
                     ylab = paste0("Y[", subpop.names[s], "]"),
                     pch = 19, col = col.use)
            } else
            {
                plot.new()
            }
        }
    }

}




