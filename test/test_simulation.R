

fig.path <- "C:/Users/Jared/Dropbox/ACO/rehospitalization/sdr_hierarchical/figures/"

hier.sir <- function(x.list, y, h = 10L, d = 2L, slice.ind = NULL)
{
    p <- ncol(x.list[[1]])
    x <- bdiag(x.list)
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
        Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

        eig.c <- eigen((diag(ncol(Pc)) - Pc) %*% V.hat )
        eta.hat <- eig.c$vectors[,1:d]
        beta.list[[c]] <- t(t(eta.hat) %*% sqrt.inv.cov)
    }

    list(beta.hat = do.call(cbind, beta.list), cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}



hier.phd <- function(x.list, y, d = 2L)
{
    p <- ncol(x.list[[1]])
    x <- bdiag(x.list)
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

    for (c in 1:length(constraints)) constraints[[c]] <- sqrt.inv.cov %*% constraints[[c]]

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x)

    beta.list <- vector(mode = "list", length = length(constraints))
    for (c in 1:length(constraints))
    {
        Pc <- constraints[[c]] %*% solve(crossprod(constraints[[c]]), t(constraints[[c]]))

        eig.c <- eigen((diag(ncol(Pc)) - Pc) %*% V.hat )
        eta.hat <- eig.c$vectors[,1:d]
        beta.list[[c]] <- t(t(eta.hat) %*% sqrt.inv.cov)
    }

    list(beta.hat = do.call(cbind, beta.list), cov = cov, sqrt.inv.cov = sqrt.inv.cov)
}


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
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     mult.a %*% beta.a, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), mult.a %*% beta.a, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     mult.b %*% beta.b, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), mult.b %*% beta.b, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     eta.ab, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), eta.ab, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])

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

pdf(paste0(fig.path, "sim_1a_nvars100_boxplots.pdf"), height = 8, width = 10)

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

pdf(paste0(fig.path, "sim_1a_directions_nvars100_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m.dir, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

dev.off()

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
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     mult.a %*% beta.a, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), mult.a %*% beta.a, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     mult.b %*% beta.b, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), mult.b %*% beta.b, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     eta.ab, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), eta.ab, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])

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
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     mult.a %*% beta.a, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), mult.a %*% beta.a, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     mult.b %*% beta.b, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), mult.b %*% beta.b, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     eta.ab, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), eta.ab, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])

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
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     mult.a %*% beta.a, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), mult.a %*% beta.a, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     mult.b %*% beta.b, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), mult.b %*% beta.b, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])

        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     eta.ab, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), eta.ab, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])

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
