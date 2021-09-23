


#' Main hierarchical sufficient dimension reduction fitting function
#'
#' @description fits hierarchically nested sufficient dimension reduction models
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
#' @param maxit maximum number of iterations for optimization routines
#' @param tol convergence tolerance for optimization routines. Defaults to \code{1e-6}
#' @param h bandwidth parameter. By default, a reasonable choice is selected automatically
#' @param opt.method optimization method to use. Available choices are
#' \code{c("lbfgs2", "lbfgs.x", "bfgs.x", "bfgs", "lbfgs", "spg", "ucminf", "CG", "nlm", "nlminb", "newuoa")}
#' @param init.method method for parameter initialization. Either \code{"random"} for random initialization or \code{"phd"}
#' for a principle Hessian directions initialization approach
#' @param vic logical value of whether or not to compute the VIC criterion for dimension determination
#' @param grassmann logical value of whether or not to enforce parameters to be on the Grassmann manifold
#' @param nn nearest neighbor parameter for \code{\link[locfit]{locfit.raw}}
#' @param maxk maxk parameter for \code{\link[locfit]{locfit.raw}}. Set to a large number if an out of vertex space error occurs.
#' @param n.random integer number of random initializations for parameters to try
#' @param nn.try vector of nearest neighbor parameters for \code{\link[locfit]{locfit.raw}} to try in random initialization
#' @param optimize.nn should \code{nn} be optimized? Not recommended
#' @param separate.nn should each subpopulation have its own \code{nn}? If \code{TRUE}, optimization takes
#' much longer. It is rarely better, so recommended to set to \code{FALSE}
#' @param constrain.none.subpop should the "none" subpopulation be constrained to be contained in every other subpopulation's
#' dimension reduction subspace? Recommended to set to \code{TRUE}
#' @param verbose should results be printed along the way?
#' @param degree degree of kernel to use
#' @param pooled should the estimator be a pooled estimator?
#' @param ... extra arguments passed to \code{\link[locfit]{locfit.raw}}
#' @return A list with the following elements
#' \itemize{
#' \item beta a list of estimated sufficient dimension reduction matrices, one for each subpopulation
#' \item beta.init a list of the initial sufficient dimension reduction matrices, one for each subpopulation -- do not use, just for the sake of comparisons
#' \item directions a list of estimated sufficient dimension reduction directions (i.e. the reduced dimension predictors/variables), one for each subpopulation.
#' These have number of rows equal to the sample size for the subpopulation and number of columns equal to the specified dimensions of the reduced dimension spaces.
#' \item y.list a list of vectors of responses for each subpopulation
#' \item z.combinations the \code{z.combinations} specified as an input
#' \item cov list of variance covariance matrices for the covariates for each subpopulation
#' \item sqrt.inv.cov list of inverse square roots of the variance covariance matrices for the covariates for each subpopulation. These are used for scaling
#' \item solver.obj object returned by the solver/optimization function
#' \item value value of the objective function at the solution
#' \item value.init value of the objective function at the initial beta (\code{beta.init}) used
#' \item vic.est.eqn the average (unpenalized) VIC value  across the r different input values. This assesses model fit
#' \item vic.eqns the individual (unpenalized) VIC values across the r input values. Not used.
#' \item vic the penalized VIC value. This is used for dimension selection, with dimensions chosen by the set of dimensions
#' that minimize this penalized vic value that trades off model complexity and model fit
#' }
#' @export
#' @examples
#'
#' library(hierSDR)
#'
#' set.seed(123)
#' dat <- simulate_data(nobs = 200, nvars = 6,
#'                      x.type = "some_categorical",
#'                      sd.y = 1, model = 2)
#'
#' x <- dat$x ## covariates
#' z <- dat$z ## factor indicators
#' y <- dat$y ## response
#'
#' dat$beta ## true coefficients that generate the subspaces
#'
#' dat$z.combinations ## what combinations of z represent different subpops
#'
#' ## correct structural dimensions:
#' dat$d.correct
#'
#' ## fit hier SPHD model:
#'
#' \donttest{
#' hiermod <- hier.sphd(x, y, z, dat$z.combinations, d = dat$d.correct,
#'                      verbose = FALSE, maxit = 250, maxk = 8200)
#'
#' ## validated inf criterion for choosing dimensions (the smaller the better)
#' hiermod$vic
#'
#'
#' cbind(hiermod$beta[[4]], NA, dat$beta[[4]])
#'
#' ## angles between estimated and true subspaces for each population:
#' mapply(function(x,y) angle(x,y), hiermod$beta, dat$beta)
#'
#' ## projection difference norm between estimated and true subspaces for each population:
#' mapply(function(x,y) projnorm(x,y), hiermod$beta, dat$beta)
#' }
#'
#'
hier.sphd <- function(x, y, z, z.combinations, d,
                      weights = rep(1L, NROW(y)),
                      maxit = 250L, tol = 1e-9,
                      h = NULL,
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
                      grassmann             = TRUE,     # constrain parameters to the Grassmann manifold?
                      nn                    = NULL,     # nn fraction. Will be chosen automatically if nn = NULL
                      nn.try                = c(0.15, 0.25, 0.5, 0.75, 0.9, 0.95), # values to try if nn = NULL (more values takes longer)
                      n.random              = 100L,
                      optimize.nn           = FALSE,    # should nn be optimized? not recommended.
                      separate.nn           = FALSE,    # should each subpopulation have a separate nn? only used if nn = NULL
                      constrain.none.subpop = TRUE,     # should we constrain S_{none} \subseteq S_{M} for all M?
                      verbose               = TRUE,     # should messages be printed out during optimization?
                      degree                = 2,        # degree of kernel
                      pooled                = FALSE,
                      maxk = 5000,
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
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h),
                                         maxk = maxk,
                                         deg = degree, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod",
                                         alpha = c(nn.val[s], best.h),
                                         maxk = maxk,
                                         deg = degree, ...)
            }


            Ey.given.xbeta[strata.idx] <- fitted(locfit.mod)
            Ey.given.xbeta.list[[s]]   <- Ey.given.xbeta[strata.idx]
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
                                         kern = "trwt", kt = "prod", maxk = maxk,
                                         alpha = c(exp(beta.vec[1]) / (1 + exp(beta.vec[1])), best.h), deg = degree, ...)
            } else
            {
                locfit.mod <- locfit.raw(x = dir.cur, y = y.list[[s]],
                                         kern = "trwt", kt = "prod", maxk = maxk,
                                         alpha = c(nn.val[s], best.h), deg = degree, ...)
            }

            # locfit.mod <- locfit.raw(x = dir.cur, y = y[strata.idx],
            #                          kern = "trwt", kt = "prod", maxk = maxk,
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
                                     kern = "trwt", kt = "prod", maxk = maxk,
                                     alpha = c(nn.val[s], best.h), deg = degree, ...)



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

    d2   <- sapply(vec2subpopMatsId(init, p, d, z.combinations), ncol)
    npar <- sum(p * d - d ^ 2)

    value.init <- est.eqn(init, nn.val = nn)

    if (init.method == "random")
    {
        n.samples <- n.random
        beta.list.phd <- beta.list.tmp <- beta.list
        best.value <- 1e99

        nn.vals <- nn.try
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
                         tol          = tol,
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
        vic1    <- vic.eqn       + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
        vic3    <- min(vic.eqns) + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))
    } else
    {
        vic.eqns <- vic.eqn <- vic3 <- NULL
    }

    vic2 <- slver$value + log(sum(nobs.vec)) * nvars * (sum(sapply(beta.mat.list, ncol)))



    model.list <- vector(mode = "list", length = 3)

    directions.list <- vector(mode = "list", length = n.combinations)

    for (m in 1:n.combinations)
    {
        directions.list[[m]] <- x.tilde[[m]] %*% beta.mat.list[[m]]
    }

    sse.vec <- mse.vec <- numeric(n.combinations)
    # if (calc.mse)
    # {
    #     for (m in 1:n.combinations)
    #     {
    #         strata.idx <- strat.idx.list[[m]]
    #         #best.h     <- best.h.vec[m]
    #         dir.cur    <- directions.list[[m]]
    #
    #         sd <- sd(dir.cur)
    #
    #         best.h <- sd * (0.75 * nrow(dir.cur)) ^ (-1 / (ncol(dir.cur) + 4) )
    #
    #         model.list[[m]] <- locfit.raw(x = dir.cur, y = y.list[[m]],
    #                                       kern = "trwt", kt = "prod",
    #                                       alpha = c(nn, best.h), deg = 2, ...)
    #
    #         fitted.vals <- fitted(model.list[[m]])
    #
    #         sse.vec[m] <-  sum((y[strata.idx] - fitted.vals) ^ 2)
    #         mse.vec[m] <- mean((y[strata.idx] - fitted.vals) ^ 2)
    #     }
    # }

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


#' Plotting hierarchical SDR models
#'
#' @description Plots hier.sdr objects
#'
#' @param x fitted object returned by \code{\link[hierSDR]{hier.sphd}}
#' @param ... not used
#' @seealso \code{\link[hierSDR]{hier.sphd}} for function which fits hierarchical SDR model
#' @return No return value, called for side effects
#' @rdname plot
#' @export
#' @examples
#'
#' library(hierSDR)
#'
plot.hier_sdr_fit <- function(x, ...)
{
    n.subpops    <- length(x$beta)
    subpop.names <- names(x$beta)
    dimensions   <- sapply(x$beta, ncol)

    maxd <- max(dimensions)

    oldpar <- par(no.readonly = TRUE)    # code line i
    on.exit(par(oldpar))            # code line i + 1

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




