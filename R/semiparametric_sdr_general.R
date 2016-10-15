

semi.phd.hier.separate.general <- function(x,
                                           y,
                                           conditions,
                                           d = NULL,
                                           maxit = 10L,
                                           h = NULL, ...)
{

    nvars <- ncol(x)
    nobs  <- nrow(x)
    y     <- drop(y)
    dimy  <- dim(y)
    leny  <- ifelse(is.null(dimy), length(y), dimy[1])
    stopifnot(leny == nobs)

    if (any(apply(conditions, 2, function(cl) length(unique(cl))) > 2))
        stop("Each column of the conditions matrix should only
             represent an indicator of one condition. Continuous
             values or non-binary values not allowed.")

    if (!is.matrix(conditions))
        stop("conditions must be a matrix. Please read the help file for details")

    if (nobs != nrow(conditions))
        stop("conditions matrix must have same number of rows as x matrix")


    ## starting and ending index of each dataset
    n.conditions <- ncol(conditions)
    M            <- 2 ^ n.conditions
    p            <- ncol(x)              # number of covariates
    N            <- nrow(x)              # total number of observations
    vnames       <- colnames(x)
    if(is.null(vnames)) vnames <- paste0("V", seq(v))

    # all possible combinations of conditions
    combin.mat   <- array(NA, dim = c(M, n.conditions))
    data.indices <- vector(mode = "list", length = M)
    group.vec    <- apply(conditions, 1, FUN=function(x) {paste(x,collapse = ",")})
    ind.2.remove <- NULL
    for (c in 1:M)
    {
        combin.mat[c,]    <- as.integer(intToBits(c)[1:n.conditions])
        data.indices[[c]] <- which(group.vec == paste(combin.mat[c,], collapse = ","))
        if (length(data.indices[[c]]) == 0)
        {
            ind.2.remove <- c(c, ind.2.remove)
        }
    }
    # remove combinations that have no observations
    if (!is.null(ind.2.remove))
    {
        data.indices[ind.2.remove] <- NULL
        M <- length(data.indices)
        combin.mat <- combin.mat[-ind.2.remove,]
    }
    # sort in order that allows fitting to work
    #ord.idx <- order(rowSums(combin.mat), decreasing = FALSE)
    #order combin.mat by each column
    ord.idx      <- do.call(order, split(combin.mat, rep(1:ncol(combin.mat), each = nrow(combin.mat))))
    combin.mat   <- combin.mat[ord.idx,]
    data.indices <- data.indices[ord.idx]

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

    ## this list will be used to create constraint matrices
    above.idx.noself.list <- above.idx.list

    for (l in 1:length(above.idx.list))
    {
        above.idx.noself.list[[l]] <- above.idx.noself.list[[l]][which(above.idx.noself.list[[l]] != l)]
    }


    ## set up data
    x.list <- y.list <- vector(mode = "list", length = M)
    for (c in 1:M)
    {
        y.list[[c]] <- y[data.indices[[c]]  ]
        x.list[[c]] <- x[data.indices[[c]], ]
    }


    if (is.null(d))
    {
        ## default to smallest
        ## possible set of dimensions
        d <- rep(1L, M)
        d[which(rowSums(combin.mat) > 1)] <- 0
    }

    d <- as.vector(data.matrix(d))
    names(d) <- NULL

    nobs.vec  <- unlist(lapply(x.list, nrow))
    nvars.vec <- unlist(lapply(x.list, ncol))
    pp <- nvars.vec[1]

    ## total number of dimensions
    D <- sum(d)

    cum.d <- c(0, cumsum(d))

    cov.list <- lapply(x.list, cov)

    x <- as.matrix(bdiag(x.list))
    cov.b <- cov(as.matrix(x))
    eig.cov.b <- eigen(cov.b)
    sqrt.inv.cov.b <- eig.cov.b$vectors %*% diag(1 / sqrt(eig.cov.b$values)) %*% t(eig.cov.b$vectors)
    x.tilde.b <- scale(x, scale = FALSE) %*% sqrt.inv.cov.b

    sqrt.inv.cov <- lapply(1:length(x.list), function(i) {
        eig.cov <- eigen(cov.list[[i]])
        eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    })

    x.tilde <- lapply(1:length(x.list), function(i) {
        scale(x.list[[i]], scale = FALSE) %*% sqrt.inv.cov[[i]]
    })



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
            locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
                                     kern = "gauss",
                                     alpha = best.h, deg = 3, ...)
            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
        }


        lhs   <- sum(unlist(lapply(1:length(x.tilde), function(i) {
            strata.idx <- which(strat.id == unique.strata[i])
            resid <- drop(y[strata.idx] - Ey.given.xbeta[strata.idx])
            norm(crossprod(x.tilde[[i]], resid * x.tilde[[i]]), type = "F") ^ 2 / (nobs.vec[i] ^ 2)
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

    beta <- beta.init

    beta.mat.list <- vector(mode = "list", length = 3L)

    beta.mat.list[[1]] <- matrix(slver$par[(cum.d[1] * p + 1):(cum.d[2] * p)], ncol = d[1])
    beta.mat.list[[2]] <- matrix(slver$par[(cum.d[2] * p + 1):(cum.d[3] * p)], ncol = d[2])
    beta.mat.list[[3]] <- matrix(slver$par[(cum.d[3] * p + 1):length(slver$par)], ncol = cum.d[4] )

    #beta.semi <- matrix(slver$par, ncol = D)
    #t(t(matrix(slver$par, ncol = d * length(constraints))) %*% sqrt.inv.cov)


    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- beta.init # t(t(beta.init) %*% sqrt.inv.cov)


    list(beta = beta.mat.list, beta.init = beta, solver.obj = slver,
         #beta.rand.init = t(t(beta.rand.init) %*% sqrt.inv.cov),
         cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}
