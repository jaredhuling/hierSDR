
Kepanechnikov  <- function(u) 0.75 * (1 - (u) ^ 2) * (abs(u) < 1)
Kepanechnikov2 <- function(u) 0.75 * (1 - (u) ^ 2)


# nwsmooth <- function(x, y, h = 1, max.points = 100)
# {
#     nobs        <- NROW(x)
#     dist.idx    <- fields.rdist.near(x, x, h, max.points = nobs * max.points)
#     dist.idx$ra <- Kepanechnikov2(dist.idx$ra / h)
#     diffmat     <- sparseMatrix(dist.idx$ind[,1], dist.idx$ind[,2],
#                                 x = dist.idx$ra,  dims = dist.idx$da)
#     #diffmat@x   <- Kepanechnikov2(diffmat@x)# / h
#
#     RSS <- sum(((Diagonal(nobs) - diffmat) %*% y)^2)/nobs
#     predicted.values <- Matrix::colSums(y * diffmat) / Matrix::colSums(diffmat)
#     trS <- sum(Matrix::diag(diffmat))
#     #gcv <- (1 / nobs) * sum(( (y - predicted.values) / (1 - trS / nobs) ) ^ 2)
#     gcv <- RSS / ((1 - trS / nobs) ) ^ 2
#     list(fitted = predicted.values, gcv = gcv)
# }
#
#
# nwsmoothcov <- function(directions, x, h = 1, max.points = 100, ret.xtx = FALSE)
# {
#     nobs <- NROW(directions)
#     dist.idx    <- fields.rdist.near(directions, directions, h,
#                                      max.points = nobs * max.points)
#     dist.idx$ra <- Kepanechnikov2(dist.idx$ra / h)
#     diffmat     <- sparseMatrix(dist.idx$ind[,1], dist.idx$ind[,2],
#                                 x = dist.idx$ra,  dims = dist.idx$da)
#
#     txpy <- predicted.values <- vector(mode = "list", length = nobs)
#     trS <- sum(Matrix::diag(diffmat))
#
#     csums <- Matrix::colSums(diffmat)
#     for (i in 1:nobs) txpy[[i]] <- tcrossprod(x[i,])
#     for (i in 1:nobs)
#     {
#         sum.cov <- txpy[[1]] * as.numeric(diffmat[1,i])
#         for (j in 2:nobs) sum.cov <- sum.cov + txpy[[j]] * as.numeric(diffmat[j,i])
#
#         predicted.values[[i]] <- as.matrix(sum.cov / csums[i])
#     }
#
#     normdiff <- 0
#     diff.cov <- vector(mode = "list", length = length(predicted.values))
#     for (j in 1:nobs)
#     {
#         diff.cov[[j]] <- txpy[[j]] - predicted.values[[j]]
#         normdiff <- normdiff + norm(diff.cov[[j]], type = "F")^2
#     }
#
#     gcv <- (1 / nobs) * ((sqrt(normdiff) / (1 - trS / nobs)) ^ 2)
#
#     if (!ret.xtx) txpy <- NULL
#
#     list(fitted = predicted.values,
#          resid  = diff.cov,
#          gcv = gcv, xtx.list = txpy)
# }
#
# nwcov <- function(directions, x, txpy.list, h = 1, max.points = 100, ret.xtx = FALSE, ncuts = 3)
# {
#     nobs <- NROW(directions)
#     diff.cov <- vector(mode = "list", length = nobs)
#
#
#     pr.vals <- seq(0, 1, length.out = 3 + ncuts)
#     pr.vals <- pr.vals[-c(1, length(pr.vals))]
#     cuts <- apply(apply(directions, 2, function(x)
#         cut(x, breaks = c(min(x) - 0.001, quantile(x, probs = pr.vals), max(x) + 0.001 )  )), 1,
#         function(r) paste(r, collapse = "|"))
#
#     unique.cuts <- unique(cuts)
#
#     for (i in 1:length(unique.cuts))
#     {
#         in.idx <- which(cuts == unique.cuts[i])
#         n.curr <- length(in.idx)
#         mean.cov <- crossprod(x[in.idx,,drop=FALSE]) / n.curr
#         for (j in in.idx)
#         {
#             diff.cov[[j]] <- txpy.list[[j]] - mean.cov
#         }
#     }
#
#
#     list(resid  = diff.cov)
# }


#' PHD SDR fitting function
#'
#' @description fits SDR models (PHD approach)
#'
#' @param x an n x p matrix of covariates, where each row is an observation and each column is a predictor
#' @param y vector of responses of length n
#' @param d an integer representing the structural dimension
#' @return A list with the following elements
#' \itemize{
#' \item beta.hat estimated sufficient dimension reduction matrix
#' \item eta.hat coefficients on the scale of the scaled covariates
#' \item cov variance covariance matric for the covariates
#' \item sqrt.inv.cov inverse square root of the variance covariance matrix for the covariates. Used for scaling
#' \item M matrix from principal Hessian directions
#' \item eigenvalues eigenvalues of the M matrix
#' }
#' @export
phd <- function(x, y, d = 5L)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    eigvals <- eig.cov$values
    eigvals[eigvals <= 0] <- 1e-5

    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eigvals)) %*% t(eig.cov$vectors)

    sqrt.cov <- eig.cov$vectors %*% diag(sqrt(eigvals)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    y.scaled <- scale(y, scale = FALSE)
    miny     <- min(y.scaled)

    V.hat    <- crossprod(x.tilde, drop(y.scaled + miny * sign(miny) + 1e-5) * x.tilde) / nrow(x)
    eig.V    <- eigen(V.hat)
    eta.hat  <- eig.V$vectors[,1:d]
    beta.hat <- t(t(eta.hat) %*% sqrt.inv.cov)
    list(beta.hat = beta.hat,
         eta.hat = eta.hat,
         M = V.hat,
         cov = cov,
         sqrt.inv.cov = sqrt.inv.cov,
         eigenvalues = eig.V$values)
}


#' Semiparametric PHD SDR fitting function
#'
#' @description fits semiparametric SDR models (PHD approach)
#'
#' @param x an n x p matrix of covariates, where each row is an observation and each column is a predictor
#' @param y vector of responses of length n
#' @param d an integer representing the structural dimension
#' @param maxit maximum number of iterations
#' @param h bandwidth parameter. By default, a reasonable choice is selected automatically
#' @param opt.method optimization method to use. Available choices are
#' \code{c("lbfgs2", "lbfgs.x", "bfgs.x", "bfgs", "lbfgs", "spg", "ucminf", "CG", "nlm", "nlminb", "newuoa")}
#' @param init.method method for parameter initialization. Either \code{"random"} for random initialization or \code{"phd"}
#' for a principle Hessian directions initialization approach
#' @param vic logical value of whether or not to compute the VIC criterion for dimension determination
#' @param nn nearest neighbor parameter for \code{\link[locfit]{locfit.raw}}
#' @param optimize.nn should \code{nn} be optimized? Not recommended
#' @param n.samples number of samples for the random initialization method
#' @param verbose should results be printed along the way?
#' @param degree degree of kernel to use
#' @param ... extra arguments passed to \code{\link[locfit]{locfit.raw}}
#' @return A list with the following elements
#' \itemize{
#' \item beta estimated sufficient dimension reduction matrix
#' \item beta.init initial sufficient dimension reduction matrix -- do not use, just for the sake of comparisons
#' \item cov variance covariance matric for the covariates
#' \item sqrt.inv.cov inverse square root of the variance covariance matrix for the covariates. Used for scaling
#' \item solver.obj object returned by the solver/optimization function
#' \item vic the penalized VIC value. This is used for dimension selection, with dimension chosen to
#' minimize this penalized vic value that trades off model complexity and model fit
#' }
#' @export
semi.phd <- function(x, y, d = 5L, maxit = 100L, h = NULL,
                     opt.method = c("lbfgs.x", "bfgs", "lbfgs2",
                                    "bfgs.x",
                                    "lbfgs",
                                    "spg",
                                    "ucminf",
                                    "CG",
                                    "nlm",
                                    "nlminb",
                                    "newuoa"),
                     nn = 0.95,
                     init.method = c("random", "phd"),
                     optimize.nn = FALSE, verbose = TRUE,
                     n.samples = 100,
                     degree = 2,
                     vic = TRUE, ...)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    nobs  <- nrow(x)
    nvars <- ncol(x)
    eigvals <- eig.cov$values
    eigvals[eigvals <= 0] <- 1e-5

    init.method <- match.arg(init.method)

    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eigvals)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov



    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x.tilde)
    eig.V <- eigen(V.hat)
    beta.init <- eig.V$vectors[,1:d,drop=FALSE]

    beta.init <- grassmannify(beta.init)$beta

    if (is.null(h))
    {
        h <- exp(seq(log(0.5), log(25), length.out = 25))
    }
    directions <- x.tilde %*% beta.init
    gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = degree, ...)[4])
    best.h.init     <- h[which.min(gcv.vals)]

    est.eqn <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {

        if (optimize.nn)
        {
            beta.mat   <- rbind(diag(d), matrix(beta.vec[-1], ncol = d))
        } else
        {
            beta.mat   <- rbind(diag(d), matrix(beta.vec, ncol = d))
        }

        directions <- x.tilde %*% beta.mat
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
        #best.h     <- best.h.init # h[which.min(gcv.vals)]
        sd <- sd(directions)

        best.h <- sd * (0.75 * nrow(directions)) ^ (-1 / (ncol(directions) + 4) )

        if (optimize.nn)
        {
            locfit.mod <- locfit.raw(x = directions, y = y,
                                     kern = "trwt", kt = "prod",
                                     alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = degree, ...)
        } else
        {
            locfit.mod <- locfit.raw(x = directions, y = y,
                                     kern = "trwt", kt = "prod",
                                     alpha = c(nn.val, best.h), deg = degree, ...)
        }


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }


    est.eqn.vic <- function(beta.mat, nn.val)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        #beta.mat   <- matrix(beta.vec, ncol = d + 1)
        directions <- x.tilde %*% beta.mat
        #gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 3, ...)[4])
        best.h     <- best.h.init # h[which.min(gcv.vals)]

        sd <- sd(directions)

        best.h <- sd * (0.75 * nrow(directions)) ^ (-1 / (ncol(directions) + 4) )

        locfit.mod <- locfit.raw(x = directions, y = y,
                                 kern = "trwt", kt = "prod",
                                 alpha = c(nn.val, best.h), deg = degree, ...)


        Ey.given.xbeta <- fitted(locfit.mod)

        resid <- drop(y - Ey.given.xbeta)
        lhs   <- norm(crossprod(x.tilde, resid * x.tilde), type = "F") ^ 2 / (nobs ^ 2)
        lhs
    }

    est.eqn.grad <- function(beta.vec, nn.val, optimize.nn = FALSE)
    {
        grad.full     <- numDeriv::grad(est.eqn, beta.vec, method = "simple", nn.val = nn.val, optimize.nn = optimize.nn)
        if (verbose) cat("grad norm: ", sqrt(sum(grad.full ^ 2)), "\n")
        grad.full
    }

    init <- as.vector(beta.init[-(1:ncol(beta.init)),])

    npar <- length(init)

    if (init.method == "random")
    {

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



    # test which nn values minimize the most effectively
    if (is.null(nn))
    {
        tryval <- try.nn(nn.vals      = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                         init         = init,
                         est.eqn      = est.eqn,
                         est.eqn.grad = est.eqn.grad,
                         opt.method   = opt.method,
                         optimize.nn  = optimize.nn,
                         maxit        = 10L,
                         verbose      = verbose)
        nn   <- tryval$nn
        init <- tryval$par


        if (init.method == "random")
        {
            tryval2 <- try.nn(nn.vals      = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                              init         = init.rand,
                              est.eqn      = est.eqn,
                              est.eqn.grad = est.eqn.grad,
                              opt.method   = opt.method,
                              optimize.nn  = optimize.nn,
                              maxit        = 10L,
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

    beta <- beta.init
    for (i in 1:maxit)
    {

    }

    beta.semi.lower <- matrix(slver$par, ncol = d)
    beta.semi <- rbind(diag(d), beta.semi.lower)

    ## calculate VIC
    if (vic)
    {
        beta.u <- beta.semi.lower[1,,drop=FALSE]
        beta.l <- beta.semi.lower[-1,,drop=FALSE]

        v_vec <- c(-2, -1, 0, 1, 2)

        eqn.vals.vic <- numeric(length(v_vec))

        for (v in 1:length(v_vec))
        {
            v.1  <- matrix(rep(v_vec[v], nrow(beta.semi.lower) - 1), ncol=1)
            vk.1 <- matrix(0, ncol=ncol(beta.semi.lower) + 1, nrow=nrow(beta.semi.lower)-1)
            vk.1[,-ncol(vk.1)] <- beta.l - drop(v.1 %*% beta.u)
            vk.1[,ncol(vk.1)]  <- v.1
            vk.1.vec <- as.vector(vk.1)

            beta.mat.vic <- rbind(diag(ncol(vk.1)), vk.1)

            eqn.vals.vic[v] <- est.eqn.vic(beta.mat.vic, nn) * nobs
        }

        vic <- mean(eqn.vals.vic) + nvars * d * log(nobs)


        # vic <- slver$value * nobs + nvars * d * log(nobs)
    } else
    {
        vic <- NULL
    }



    #beta.semi <- t(t(beta.semi) %*% sqrt.inv.cov)
    beta      <- t(t(beta.init) %*% sqrt.inv.cov)


    directions <- x %*% beta.semi
    gcv.vals   <- sapply(h, function(hv) gcv(x = directions, y = y, alpha = hv, deg = 2, ...)[4])
    best.h     <- h[which.min(gcv.vals)]
    locfit.mod <- locfit.raw(x = directions, y = y, alpha = best.h, deg = 2, ...)

    list(beta         = beta.semi,
         beta.init    = beta,
         solver.obj   = slver,
         cov          = cov,
         sqrt.inv.cov = sqrt.inv.cov,
         final.gcv    = min(gcv.vals),
         final.model  = locfit.mod,
         all.gcvs     = gcv.vals,
         #rsq          = rsq.mean,
         nn           = nn,
         vic          = vic)
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





#' Main hierarchical SDR fitting function
#'
#' @description fits hierarchical SDR models
#'
#' @param x an n x p matrix of covariates, where each row is an observation and each column is a predictor
#' @param y vector of responses of length n
#' @param z an n x C matrix of binary indicators, where each column is a binary variable indicating the presence
#' of a binary variable which acts as a stratifying variable. Each combination of all columns of \code{z} pertains
#' to a different subpopulation. WARNING: do not use too many binary variables in \code{z} or else it will quickly
#' result in subpopulations with no observations
#' @param z.combinations a matrix of dimensions 2^C x C with each row indicating a different combination of the possible
#' values in \code{z}. Each combination represents a subpopulation. This is necessary because we need to specify a
#' different structural dimension for each subpopulation, so we need to know the ordering of the subpopulations so we
#' can assign each one a structural dimension
#' @param d an integer vector of length 2^C of structural dimensions. Specified in the same order as the rows in
#' \code{z.combinations}
#' @param weights vector of observation weights
#' @param constrain.none.subpop should the "none" subpopulation be constrained to be contained in every other subpopulation's
#' dimension reduction subspace? Recommended to set to \code{TRUE}
#' @param pooled should the estimator be a pooled estimator?
#' @param ... not used
#' @return A list with the following elements
#' \itemize{
#' \item beta a list of estimated sufficient dimension reduction matrices, one for each subpopulation
#' \item directions a list of estimated sufficient dimension reduction directions (i.e. the reduced dimension predictors/variables), one for each subpopulation.
#' These have number of rows equal to the sample size for the subpopulation and number of columns equal to the specified dimensions of the reduced dimension spaces.
#' \item y.list a list of vectors of responses for each subpopulation
#' \item z.combinations the \code{z.combinations} specified as an input
#' \item cov list of variance covariance matrices for the covariates for each subpopulation
#' \item sqrt.inv.cov list of inverse square roots of the variance covariance matrices for the covariates for each subpopulation. These are used for scaling
#' }
#' @export
#' @examples
#'
#' library(hierSDR)
#'
hier.phd.nt <- function(x, y, z, z.combinations, d,
                        weights = rep(1L, NROW(y)),
                        constrain.none.subpop = TRUE,     # should we constrain S_{none} \subseteq S_{M} for all M?
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

    strat.id  <- unlist(lapply(1:length(x.list), function(id) rep(id, nrow(x.list[[id]]))))

    V.hat     <- crossprod(x.tilde.b, drop(scale(y, scale = FALSE)) * x.tilde.b) / nrow(x.tilde.b)

    beta.list <- beta.init.list <- Proj.constr.list <- vector(mode = "list", length = length(constraints))



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
        }
    }

    beta.component.init <- lapply(beta.component.init, function(bb) {
        if (!is.null(bb))
        {
            nc <- ncol(bb)

            whichzero <- apply(bb, 2, function(cl) all(abs(cl) <= 1e-12))

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


    beta.mat.list  <- vec2subpopMatsId(init, p, d, z.combinations)

    ncuts <- 3




    #######    model fitting to determine bandwidth     ########






    model.list <- vector(mode = "list", length = 3)

    directions.list <- vector(mode = "list", length = n.combinations)

    for (m in 1:n.combinations)
    {
        directions.list[[m]] <- x.tilde[[m]] %*% beta.mat.list[[m]]
    }

    names(beta.mat.list) <- names(directions.list) <- combinations


    ret <- list(beta           = beta.mat.list,
                directions     = directions.list,
                y.list         = y.list,
                z.combinations = z.combinations,
                cov            = cov,
                sqrt.inv.cov   = sqrt.inv.cov)
    ret
}

