
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
    list(beta.hat = beta.hat, eta.hat = eta.hat, M = V.hat, cov = cov, sqrt.inv.cov = sqrt.inv.cov, eigenvalues = eig.V$values)
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
                     nn = NULL,
                     init.method = c("random", "phd"),
                     optimize.nn = FALSE, verbose = TRUE,
                     n.samples = 100,
                     degree = 2,
                     vic = FALSE, ...)
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


    est.eqn.vic <- function(beta.vec, nn.val)
    {
        #beta.mat   <- rbind(beta.init[1:d,], matrix(beta.vec, ncol = d))
        beta.mat   <- matrix(beta.vec, ncol = d + 1)
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

    beta.semi <- rbind(diag(d), matrix(slver$par, ncol = d))

    ## calculate VIC
    if (vic)
    {
        beta.u <- beta.semi[1,,drop=FALSE]

        v_vec <- c(-2, -1, 0, 1, 2)

        eqn.vals.vic <- numeric(length(v_vec))

        for (v in 1:v_vec)
        {
            v.1  <- matrix(rep(v_vec[v], nrow(beta.semi) - 1), ncol=1)
            vk.1 <- matrix(0, ncol=ncol(beta.semi) + 1, nrow=nrow(beta.semi))
            vk.1[-1,-ncol(vk.1)] <- vk.1[-1,-ncol(vk.1)] - drop(v.1 %*% beta.u)
            vk.1[-1,ncol(vk.1)]  <- v.1
            vk.1[1,] <- c(rep(0, ncol(vk.1) - 1), 1)
            vk.1.vec <- as.vector(vk.1)

            eqn.vals.vic[v] <- est.eqn.vic(vk.1.vec) * nobs
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
