


newton <- function(par, fn, constr, maxit = 100L, tol = 1e-5, verbose = TRUE, ...)
{
    func.hist <- numeric(maxit)
    grad.hist <- numeric(maxit)

    beta <- par
    niter <- maxit
    for (i in 1:maxit)
    {
        g     <- grad(fn, x = beta, ...)
        H     <- hessian(fn, x = beta, ...)
        delta <- -solve(H, g)


        beta <- beta + delta

        func.hist[i] <- fn(beta, ...)
        grad.hist[i] <- sqrt(sum((g^2)))
        if (verbose) cat("grad norm: ", grad.hist[i], "\n")
        if (sum(grad.hist[i] <= tol) ||
            all(abs(delta) <= tol))
        {
            break
        }
    }
    func.hist <- func.hist[1:i]
    grad.hist <- grad.hist[1:i]
    niter <- i
    list(par = beta, niter = niter, value = func.hist[i],
         func.hist = func.hist, grad.norm = grad.hist)
}

newton.constr <- function(par, fn, constr.mat, maxit = 100L, tol = 1e-5, ...)
{
    func.hist <- numeric(maxit)
    grad.hist <- numeric(maxit)

    beta <- par
    niter <- maxit
    for (i in 1:maxit)
    {
        g      <- grad(fn, x = beta, ...)
        H      <- hessian(fn, x = beta, ...)

        constr <- constr.mat %*% beta

        nc <- nrow(constr)

        rhs <- -rbind(matrix(g, ncol = 1), constr)
        lhs <- cbind(rbind(H, constr.mat),
                     rbind(t(constr.mat),
                           matrix(0, ncol = nrow(constr.mat), nrow = nrow(constr.mat))))
        soln <- solve(lhs, rhs)

        delta <- soln[1:length(beta)]

        beta <- beta + delta

        func.hist[i] <- fn(beta, ...)
        grad.hist[i] <- sqrt(sum((g^2)))
        if (sum(grad.hist[i] <= tol) ||
            all(abs(delta) <= tol))
        {
            break
        }
    }
    func.hist <- func.hist[1:i]
    grad.hist <- grad.hist[1:i]
    niter <- i
    list(beta = beta, niter = niter, value = func.hist[i],
         func.hist = func.hist, grad.norm = grad.hist)
}


createBelowList <- function(combin.mat)
{

    M <- nrow(combin.mat)

    below.idx.list <- vector(mode = "list", length = M)

    for (c in 1:M)
    {
        indi <- which(combin.mat[c,] == 1)

        # get all indices of g terms which are directly below the current var
        # in the hierarchy, ie. for three categories A, B, C, for the A group,
        # this will return the indices for AB and AC. For AB, it will return
        # the index for ABC. for none, it will return the indices for
        # A, B, and C, etc
        inner.loop <- (1:(M))[-c]
        for (j in inner.loop)
        {
            diffs.tmp <- combin.mat[j,] - combin.mat[c,]
            if (all( diffs.tmp <= 0 ))
            {
                below.idx.list[[c]] <- c(below.idx.list[[c]], j)
            }
        }
        below.idx.list[[c]] <- c(below.idx.list[[c]], c)
    }
    rsc <- rowSums(combin.mat)
    if (any(rsc == 0))
    {
        below.idx.list[[which(rsc == 0)]] <- which(rsc == 0)
    }
    below.idx.list
}


subpop.struct <- function(n.conditions = 2, incl.none = FALSE)
{
    cond.mat <- expand.grid(rep(list(c(0,1)), n.conditions))

    rsm <- rowSums(cond.mat)
    cond.mat <- cond.mat[order(rsm),]

    if (!incl.none)
    {
        cond.mat <- cond.mat[-which(rowSums(cond.mat) == 0),]
    }
    cond.mat
}

vec2subpopMats <- function(vec, p, d, cond.mat, incl.none = FALSE)
{
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)
    component.list  <- vector(mode = "list", length = n.subpops)

    createBelowList <- createBelowList(cond.mat)

    cumul.params    <- c(0, cumsum(d * p))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            component.list[[s]] <- matrix(vec[vec.idx], ncol = d[s])
        }
    }

    for (s in 1:n.subpops)
    {
        mat.list[[s]] <- do.call(cbind, component.list[createBelowList[[s]]])
    }
    mat.list
}

vec2subpopMatsId <- function(vec, p, d, cond.mat, incl.none = FALSE)
{
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)
    component.list  <- vector(mode = "list", length = n.subpops)

    createBelowList <- createBelowList(cond.mat)

    cumul.params    <- c(0, cumsum(d * p - d ^ 2))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            component.list[[s]] <- rbind(diag(d[s]),
                                         matrix(vec[vec.idx], ncol = d[s]))
        }
    }

    for (s in 1:n.subpops)
    {
        mat.list[[s]] <- do.call(cbind, component.list[createBelowList[[s]]])
    }
    mat.list
}


vec2subpopMatsIdVIC <- function(vec, p, d, cond.mat, incl.none = FALSE, v = 1)
{
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)
    component.list  <- vector(mode = "list", length = n.subpops)

    createBelowList <- createBelowList(cond.mat)

    cumul.params    <- c(0, cumsum(d * p - d ^ 2))
    dim.added <- FALSE
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx   <- (cumul.params[s] + 1):cumul.params[s + 1]
            component.list[[s]] <- rbind(diag(d[s]),
                                             matrix(vec[vec.idx], ncol = d[s]))
        }
    }

    for (s in 1:n.subpops)
    {
        mat.list[[s]] <- do.call(cbind, component.list[createBelowList[[s]]])
    }

    change.idx <- unname(which.min(rowSums(cond.mat)))

    beta.orig <- mat.list[[change.idx]][-(1:ncol(mat.list[[change.idx]])),,drop = FALSE]
    beta.U    <- beta.orig[1,,drop = FALSE]
    beta.L    <- beta.orig[-1,,drop = FALSE]
    V         <- matrix(v,  nrow = nrow(beta.L), ncol = 1)
    beta.new  <- cbind(beta.L - V %*% beta.U, V)
    mat.list[[change.idx]] <- rbind(diag(ncol(beta.orig) + 1),
                                 beta.new)

    mat.list
}


vec <- rnorm(20)

p <- 10
d <- c(1, 1, 0)

cond.mat <- subpop.struct(2)
vec2subpopMats(vec, p, d, cond.mat)

subpopMats2vec <- function(param.list, d, incl.none = FALSE)
{
    n.conditions <- as.integer(log2(length(d) + !incl.none))
    n.subpops    <- length(d)
    p            <- nrow(param.list[[1]])

    vec.list <- vector(mode = "list", length = length(param.list))

    for (s in 1:length(param.list))
    {
        if (d[s] > 0)
        {
            print(str(param.list[[s]]))
            param.list[[s]] <- as.matrix(param.list[[s]])
            nc <- NCOL(param.list[[s]])
            print((nc - d[s] + 1):nc)
            vec.list[[s]] <- param.list[[s]][,(nc - d[s] + 1):nc,drop=FALSE]
        }
    }
    unlist(vec.list)
}

createAboveList <- function(combin.mat)
{

    M <- nrow(combin.mat)

    direct.above.idx.list <- above.idx.list <- vector(mode = "list", length = M)

    for (c in 1:M)
    {
        indi <- which(combin.mat[c,] == 1)

        # get all indices of g terms which are directly above the current var
        # in the hierarchy, ie. for three categories A, B, C, for the A group,
        # this will return the indices for AB and AC. For AB, it will return
        # the index for ABC. for none, it will return the indices for
        # A, B, and C, etc
        inner.loop <- (1:(M))[-c]
        for (j in inner.loop)
        {
            diffs.tmp <- combin.mat[j,] - combin.mat[c,]
            if (all( diffs.tmp >= 0 ))
            {
                if (sum( diffs.tmp == 1 ) == 1)
                {
                    direct.above.idx.list[[c]] <- c(direct.above.idx.list[[c]], j)
                }
                above.idx.list[[c]] <- c(above.idx.list[[c]], j)
            }
        }
        above.idx.list[[c]] <- c(above.idx.list[[c]], c)
    }
    rsc <- rowSums(combin.mat)
    if (any(rsc == 0))
    {
        above.idx.list[[which(rsc == 0)]] <- which(rsc == 0)
    }
    above.idx.list
}




semi.phd.hier.newton <- function(x.list, y, d = rep(1L, 3L),
                                 maxit = 10L, h = NULL, B = NULL, vic = FALSE,
                                 weights = rep(1L, NROW(y)),
                                 opt.method = c("lbfgs.x", "bfgs", "lbfgs2",
                                                "bfgs.x",
                                                "lbfgs",
                                                "spg", "ucminf"),
                                 init.method = c("phd", "random"),
                                 optimize.nn = TRUE,
                                 nn = NULL, calc.mse = FALSE,
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




    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))

    if (init.method == "phd")
    {
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

        Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]

        #eig.V <- eigen(V.hat)
        beta.init  <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]
    } else
    {
        for (c in 1:length(constraints))
        {
            #print(d)
            if (d[c] > 0)
            {
                Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

                Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

                beta.list[[c]] <- Proj.constr.list[[c]] %*%
                    matrix(runif(prod(dim(eta.hat)),
                                 min = min(beta.list[[c]]),
                                 max = max(beta.list[[c]])),
                           ncol = NCOL(eta.hat))
            }
        }

        Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]
        beta.init <- do.call(cbind, beta.list)
    }



    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)

    beta.list.init <- beta.list

    for (s in 1:length(beta.list.init))
    {
        beta.list.init[[s]] <- beta.list.init[[s]][(p * (s - 1) + 1):(p * s),,drop=FALSE]
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
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]

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
                                         alpha = c(nn.val, best.h), deg = 2, ...)
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

    est.eqn.grad <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- grad(est.eqn, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
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

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1/(ncol(dir.cur)+4) )

            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                     kern = "trwt", kt = "prod",
                                     alpha = c(nn.val, best.h), deg = 2, ...)

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



    beta.list.init <- lapply(beta.list.init, function(bb) {
        nc <- ncol(bb)
        bb[(nc + 1):nrow(bb),]
    })

    # slver <- newton(par = unlist(beta.list.init), # subpopMats2vec(beta.list, d),
    #                 fn = est.eqn, maxit = maxit, tol = 1e-8,
    #                 verbose = verbose)

    init <- unlist(beta.list.init)

    # test which nn values minimize the most effectively
    if (is.null(nn))
    {
        nn <- try.nn(nn.vals      = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95),
                     init         = init,
                     est.eqn      = est.eqn,
                     est.eqn.grad = est.eqn.grad,
                     opt.method   = opt.method,
                     optimize.nn  = optimize.nn,
                     maxit        = 10L,
                     verbose      = verbose)

        if (verbose) print(paste("best nn:", nn))
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

    #beta.mat.list <- vector(mode = "list", length = 3L)

    #beta.mat.list[[1]] <- matrix(slver$par[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
    #beta.mat.list[[2]] <- matrix(slver$par[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
    #beta.mat.list[[3]] <- matrix(slver$par[(cum.d[3] * p + 1):length(slver$par)], ncol = cum.d[4] )


    cond.mat <- subpop.struct(2L)
    #beta.mat.list <- vec2subpopMats(slver$par, p, d, cond.mat)
    beta.mat.list <- vec2subpopMatsId(slver$par, p, d, cond.mat)

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
                                          alpha = c(nn, best.h), deg = 2, ...)

            fitted.vals <- fitted(model.list[[m]])

            sse.vec[m] <-  sum((y[strata.idx] - fitted.vals) ^ 2)
            mse.vec[m] <- mean((y[strata.idx] - fitted.vals) ^ 2)
        }
    }


    #beta.semi <- matrix(slver$par, ncol = D)
    #t(t(matrix(slver$par, ncol = d * length(constraints))) %*% sqrt.inv.cov)


    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.mat.list, beta.init = beta, solver.obj = slver,
         #beta.rand.init = t(t(beta.rand.init) %*% sqrt.inv.cov),
         cov = cov, sqrt.inv.cov = sqrt.inv.cov,
         nn  = nn,
         vic.est.eqn = vic.eqn, vic.eqns = vic.eqns,
         vic = vic, vic2 = vic2,
         sse = sse.vec, mse = mse.vec)
}

hasAllConds <- function(checkvec, combvec)
{
    comb <- which(combvec != 0)
    all(checkvec[comb] != 0)
}


hier.s.phd <- function(x, y, z, z.combinations, d,
                       maxit = 250L, h = NULL, B = NULL, vic = FALSE,
                       weights = rep(1L, NROW(y)),
                       opt.method = c("lbfgs.x", "bfgs", "lbfgs2",
                                      "bfgs.x",
                                      "lbfgs",
                                      "spg", "ucminf"),
                       init.method = c("phd", "random"),
                       nn = NULL,
                       optimize.nn = TRUE,
                       calc.mse = FALSE,
                       constrain.none.subpop = FALSE,
                       verbose = TRUE, ...)
{

    nobs <- NROW(x)

    if (nobs != NROW(z)) stop("number of observations in x doesn't match number of observations in z")

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
                # then we apply equality constraints (but only if it's not the none subpop (unless we want to constrain the none subpop))

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


    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))


    for (c in 1:length(constraints))
    {
        #print(d)
        if (d[c] > 0)
        {
            Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

            Proj.constr.list[[c]] <- diag(ncol(Pc)) - Pc

            eig.c          <- eigen(Proj.constr.list[[c]] %*% V.hat )
            eta.hat        <- eig.c$vectors[,1:d[c], drop=FALSE]
            beta.list[[c]] <- Re(eta.hat)

            if (init.method == "random")
            {
                beta.list[[c]] <- Proj.constr.list[[c]] %*%
                    matrix(runif(prod(dim(eta.hat)),
                                 min = min(beta.list[[c]]),
                                 max = max(beta.list[[c]])),
                           ncol = NCOL(eta.hat))
            }
        }
    }


    Proj.constr.list <- Proj.constr.list[!sapply(Proj.constr.list, is.null)]

    #eig.V <- eigen(V.hat)
    beta.init      <- do.call(cbind, beta.list) #eig.V$vectors[,1:d,drop=FALSE]

    unique.strata <- unique(strat.id)
    num.strata    <- length(unique.strata)

    beta.component.init <- beta.list
    for (s in 1:length(beta.list))
    {
        beta.component.init[[s]] <- beta.list[[s]][(p * (s - 1) + 1):(p * s),,drop = FALSE]
    }



    #######    model fitting to determine bandwidth     ########

    best.h.vec <- numeric(num.strata)

    if (is.null(h))
    {
        h <- exp(seq(log(0.1), log(25), length.out = 25))
    }

    beta.init.cov <- t(t(beta.init) %*% sqrt.inv.cov.b)
    for (s in 1:n.combinations)
    {
        strata.idx <- strat.idx.list[[s]]
        dir.cur    <- x.tilde.b[strata.idx, ((s - 1) * pp + 1):(s * pp)] %*%
            beta.init[((s - 1) * pp + 1):(s * pp),,drop=FALSE]

        # remove irrelevant dimensions
        dir.cur    <- dir.cur[, which(apply(dir.cur, 2, sd) != 0)]

        gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
                                                 y = y[strata.idx],
                                                 kern = "trwt", kt = "prod",
                                                 alpha = hv, deg = 3, ...)[4])
        best.h.vec[s]     <- h[which.min(gcv.vals)]
    }

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

        for (s in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[s]]
            dir.cur    <- x.tilde[[s]] %*% beta.mat.list[[s]]

            # probably not needed, but just in case
            dir.cur    <- dir.cur[,which(apply(dir.cur, 2, sd) != 0), drop = FALSE]


            #gcv.vals   <- sapply(h, function(hv) gcv(x = dir.cur,
            #                                         y = y[strata.idx],
            #                                         kern = "gauss",
            #                                         alpha = hv, deg = 3, ...)[4])
            best.h     <- best.h.vec[s] # h[which.min(gcv.vals)]

            sd <- sd(dir.cur)

            best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )

            if (optimize.nn)
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = 2, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val, best.h), deg = 2, ...)
            }

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod",
            #                          alpha = best.h, deg = 2, ...)

            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }



        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- strat.idx.list[[i]]
            resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
            wts.cur    <- weights.list[[i]]
            norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        })))
        lhs / sum(weights)
    }

    est.eqn.grad <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- grad(est.eqn, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
        if (verbose) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }


    est.eqn.vic <- function(beta.vec, nn.val, v = 1)
    {

        beta.mat.list <- vec2subpopMatsIdVIC(beta.vec, p, d, z.combinations, v = v)

        Ey.given.xbeta <- numeric(nobs)

        for (s in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[s]]
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

            locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                     kern = "trwt", kt = "prod",
                                     alpha = c(nn.val, best.h), deg = 1, ...)

            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }


        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- strat.idx.list[[i]]
            resid      <- drop(y.list[[i]] - Ey.given.xbeta[strata.idx])
            wts.cur    <- weights.list[[i]]
            norm(crossprod(x.tilde[[i]], (wts.cur * resid) * x.tilde[[i]]), type = "F") ^ 2
        })))
        lhs / sum(weights)
    }

    beta.component.init <- lapply(beta.component.init, function(bb) {
        if (!is.null(bb))
        {
            nc <- ncol(bb)
            bb[(nc + 1):nrow(bb),]
        } else
        {
            NULL
        }
    })

    init <- unlist(beta.list.init)

    # test which nn values minimize the most effectively
    if (is.null(nn))
    {
        nn <- try.nn(nn.vals      = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                     init         = init,
                     est.eqn      = est.eqn,
                     est.eqn.grad = est.eqn.grad,
                     opt.method   = opt.method,
                     optimize.nn  = optimize.nn,
                     maxit        = 10L,
                     verbose      = verbose)

        if (verbose) print(paste("best nn:", nn))
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

    vic.eqns <- c(est.eqn.vic(slver$par, nn.val = nn, v = -1),
                  est.eqn.vic(slver$par, nn.val = nn, v = -0.5),
                  est.eqn.vic(slver$par, nn.val = nn, v = 0.1),
                  est.eqn.vic(slver$par, nn.val = nn, v = 0.5),
                  est.eqn.vic(slver$par, nn.val = nn, v = 1))

    vic.eqn <- mean(vic.eqns)

    vic  <- vic.eqn       + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    vic2 <- slver$value   + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    vic3 <- min(vic.eqns) + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))

    model.list <- vector(mode = "list", length = 3)


    sse.vec <- mse.vec <- numeric(n.combinations)
    if (calc.mse)
    {
        for (m in 1:n.combinations)
        {
            strata.idx <- strat.idx.list[[m]]
            best.h     <- best.h.vec[m]
            dir.cur    <- x.tilde[[m]] %*% beta.mat.list[[m]]

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


    list(beta = beta.mat.list,
         beta.init = beta.init,
         cov = cov, sqrt.inv.cov = sqrt.inv.cov,
         solver.obj = slver,
         vic.est.eqn = vic.eqn, vic.eqns = vic.eqns,
         vic = vic, vic2 = vic2, vic3 = vic3,
         sse = sse.vec, mse = mse.vec)
}

