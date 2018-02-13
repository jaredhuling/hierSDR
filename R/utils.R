


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

    BelowList <- createBelowList(cond.mat)

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
        mat.list[[s]] <- do.call(cbind, component.list[BelowList[[s]]])
    }
    mat.list
}

subpopMats2components <- function(mat.list, p, d, cond.mat, incl.none = FALSE)
{
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    nnz             <- sum(d > 0)
    component.list  <- vector(mode = "list", length = n.subpops)

    BelowList      <- createBelowList(cond.mat)

    cumul.params    <- c(0, cumsum(d * p))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            component.list[[s]] <- matrix(vec[vec.idx], ncol = d[s])
        }
    }

    BelowListCols <- lapply(BelowList, function(mt) {
        vl <- NULL
        for (i in 1:length(mt))
        {
            vl <- c(vl, rep(mt[i], d[mt[i]]))
        }
        vl
    })


    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            subpops.related <- rep(FALSE, n.subpops)
            for (si in 1:n.subpops)
            {
                if (s %in% BelowList[[si]])
                {
                    subpops.related[si] <- TRUE
                }
            }
            subpops.related <- which(subpops.related)
            if (length(subpops.related))
            {
                component.list[[s]] <- vector(mode = "list", length = length(subpops.related))
                for (si in 1:length(subpops.related))
                {
                    idx.cur <- subpops.related[si]
                    component.list[[s]][[si]] <- mat.list[[idx.cur]][,which(BelowListCols[[idx.cur]] == s),drop = FALSE]
                }
            }
        }
    }

    component.list
}


vec2subpopMatsUnconstr <- function(vec, p, d, cond.mat, incl.none = FALSE)
{
    # d here needs to be the total ncols for each beta matrix
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)

    createBelowList <- createBelowList(cond.mat)


    cumul.params    <- c(0, cumsum(d * p))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            mat.list[[s]] <- matrix(vec[vec.idx], ncol = d[s])
        }
    }
    mat.list
}

vec2subpopMatsUnconstrConstr <- function(vec, p, d, cond.mat, incl.none = FALSE)
{
    # d here needs to be the total ncols for each beta matrix
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)

    BelowList <- createBelowList(cond.mat)


    cumul.params    <- c(0, cumsum(d * p))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            mat.list[[s]] <- matrix(vec[vec.idx], ncol = d[s])
        }
    }

    component.list <- subpopMats2components(mat.list, p, d, cond.mat, incl.none)

    for (i in 1:length(component.list))
    {
        if (!is.null(component.list[[i]]))
        {
            for (s in 1:length(component.list[[i]]))
            {
                if (s == 1)
                {
                    tmp <- component.list[[i]][[s]]
                } else
                {
                    tmp <- tmp + component.list[[i]][[s]]
                }
            }
            tmp <- tmp / length(component.list[[i]])
            component.list[[i]] <- tmp
        }
    }

    for (s in 1:n.subpops)
    {
        mat.list[[s]] <- do.call(cbind, component.list[BelowList[[s]]])
    }
    mat.list
}


Unconstrvec2Constrvec <- function(vec, p, d, cond.mat, incl.none = FALSE)
{
    # d here needs to be the total ncols for each beta matrix
    n.conditions    <- as.integer(log2(length(d) + !incl.none))
    n.subpops       <- length(d)
    #cond.mat       <- subpop.struct(n.conditions, incl.none)
    mat.list        <- vector(mode = "list", length = n.subpops)
    nnz             <- sum(d > 0)

    BelowList <- createBelowList(cond.mat)


    cumul.params    <- c(0, cumsum(d * p))
    for (s in 1:n.subpops)
    {
        if (d[s] > 0)
        {
            vec.idx <- (cumul.params[s] + 1):cumul.params[s + 1]
            mat.list[[s]] <- matrix(vec[vec.idx], ncol = d[s])
        }
    }

    component.list <- subpopMats2components(mat.list, p, d, cond.mat, incl.none)

    for (i in 1:length(component.list))
    {
        if (!is.null(component.list[[i]]))
        {
            for (s in 1:length(component.list[[i]]))
            {
                if (s == 1)
                {
                    tmp <- component.list[[i]][[s]]
                } else
                {
                    tmp <- tmp + component.list[[i]][[s]]
                }
            }
            tmp <- tmp / length(component.list[[i]])
            component.list[[i]] <- tmp[-c(1:ncol(tmp)),]
        }
    }

    unlist(component.list)
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

    change.idx <- unname(which.min(rowSums(cond.mat)))

    beta.orig <- component.list[[change.idx]][-(1:ncol(component.list[[change.idx]])),,drop = FALSE]
    beta.U    <- beta.orig[1,,drop = FALSE]
    beta.L    <- beta.orig[-1,,drop = FALSE]
    V         <- matrix(v,  nrow = nrow(beta.L), ncol = 1)
    beta.new  <- cbind(beta.L - V %*% beta.U, V)
    component.list[[change.idx]] <- rbind(diag(ncol(beta.orig) + 1),
                                          beta.new)

    for (s in 1:n.subpops)
    {
        mat.list[[s]] <- do.call(cbind, component.list[createBelowList[[s]]])
    }

    # change.idx <- unname(which.min(rowSums(cond.mat)))
    #
    # beta.orig <- mat.list[[change.idx]][-(1:ncol(mat.list[[change.idx]])),,drop = FALSE]
    # beta.U    <- beta.orig[1,,drop = FALSE]
    # beta.L    <- beta.orig[-1,,drop = FALSE]
    # V         <- matrix(v,  nrow = nrow(beta.L), ncol = 1)
    # beta.new  <- cbind(beta.L - V %*% beta.U, V)
    # mat.list[[change.idx]] <- rbind(diag(ncol(beta.orig) + 1),
    #                              beta.new)

    mat.list
}



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


