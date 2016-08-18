

fig.path <- "C:/Users/Jared/Dropbox/ACO/rehospitalization/sdr_hierarchical/figures/"

set.seed(123)
y <- runif(200, -pi, pi)
x <- sin(y) + rnorm(200, sd = 0.25)

bandwidths <- exp(seq(log(0.25), log(3), length.out = 25))
nw.fits <- vector(mode = "list", length = length(bandwidths))
for (b in 1:length(bandwidths)) nw.fits[[b]] <- nwsmooth(x = y, y = x, h = bandwidths[b])
gcvs <- sapply(nw.fits, function(f) f$gcv)
bandwidths[which.min(gcvs)]

nw.fit <- nw.fits[[which.min(gcvs)]]
plot(x = y, y = x)
ordy <- order(y)
lines(x = y[ordy], y = nw.fit$fitted[ordy], col = "blue", lwd = 2)

set.seed(123)
nvars <- 10
nobs <- 250
nobs.test <- 1000
dd <- 3
x <- matrix(runif(nobs * nvars, min = -2, max = 2), ncol = nvars)
y <- (apply(x, 1, function(xr) sum( (xr[1] + xr[2])^2 ))) +
    (apply(x, 1, function(xr) ( (xr[3] + xr[4])^2 ))) +
    (apply(x, 1, function(xr) ( sum(exp(xr[5:min(length(xr), 10)])) ))) +
    rnorm(nobs, sd = 0.25)

beta <- matrix(runif(nvars * dd, min = -1, max = 1), ncol = dd)


y.true <- 1 * (x %*% beta[,1])^2 + 1 * sin(x %*% beta[,2]) * (x %*% beta[,2])^2 + 0.05 * exp(x %*% beta[,3])
y <- y.true + rnorm(nobs, sd = 2)

var(y.true) / var(y - y.true)


x.test <- matrix(runif(nobs.test * nvars, min = -2, max = 2), ncol = nvars)
y.test <- (apply(x.test, 1, function(xr) sum( (xr[1] + xr[2])^2 ))) +
    (apply(x.test, 1, function(xr) ( (xr[3] + xr[4])^2 ))) +
    (apply(x.test, 1, function(xr) ( sum(exp(xr[5:min(length(xr), 10)])) ))) +
    rnorm(nobs.test, sd = 0.25)

y.true.test <- 1 * (x.test %*% beta[,1])^2 + 1 * sin(x.test %*% beta[,2]) * (x.test %*% beta[,2])^2 + 0.05 * exp(x.test %*% beta[,3])
y.test <- y.true.test + rnorm(nobs.test, sd = 2)

ordy <- order(y)


h <- seq(1, 15, length.out = 15)
system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = x[,1:5], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.mod <- smooth.lf(x = x[,1:5], y = y, alpha = best.h, xev = x.test[,1:5]))

df.np <- data.frame(y = y, x[,1:5])

system.time(bw.obj <- npregbw(xdat = x[,1:5], ydat = y, bwmethod = "cv.ls",
                              bandwidth.compute = FALSE,
                              bws = h))

npr <- npreg(bws = bw.obj, txdat = x[,1:5], tydat = y, exdat = x.test[,1:5])

plot(x = (locfit.mod$y), y = y.test); mean((locfit.mod$y - y.test) ^ 2)
plot(x = npr$mean, y = y.test); mean((npr$mean - y.test) ^ 2)

bandwidths <- exp(seq(log(0.25), log(30), length.out = 25))[10:25]
nw.cov.fits <- vector(mode = "list", length = length(bandwidths))
for (b in 1:length(bandwidths)) nw.cov.fits[[b]] <- nwsmoothcov(x = x, y = y, h = bandwidths[b])
gcvs.cov <- sapply(nw.cov.fits, function(f) f$gcv)
bandwidths[which.min(gcvs.cov)]

nw.cov.fit <- nw.cov.fits[[which.min(gcvs.cov)]]

#    nw.cov.fit$fitted[ordy]

cov2plot <- sapply(nw.cov.fit$fitted[ordy], function(x) x[9,3])

plot(x = y[ordy], y = cov2plot, type = "b")





cor.directions <- function(a, b, x)
{
    cov <- cov(x)
    R.sq <- as.vector(crossprod(a, cov %*% b) ^ 2 / (crossprod(a, cov %*% a) * crossprod(b, cov %*% b)) )
    R.sq
}


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

beta.hat       <- sir(x, y, h = 20, d = 3)
beta.hat.phd   <- phd(x, y, d = 3)
beta.hat.sphd  <- semi.phd(x,  drop(y), d = 3, h = seq(1, 25, length.out = 25), maxit = 151)
beta.hat.sphd2 <- semi.phd2(x, drop(y), d = 3, h = seq(1, 25, length.out = 25), maxit = 101)


beta.hat.sphd$solver.obj$value

mean(sapply(1:dd, function(didx) cor.directions(beta.hat$beta.hat[,didx],     beta[,didx], x.test )))
mean(sapply(1:dd, function(didx) cor.directions(beta.hat.phd$beta.hat[,didx], beta[,didx], x.test )))
mean(sapply(1:dd, function(didx) cor.directions(beta.hat.sphd$beta[,didx],    beta[,didx], x.test )))
mean(sapply(1:dd, function(didx) cor.directions(beta.hat.sphd2$beta[,didx],    beta[,didx], x.test )))


res.vec.sir <- res.vec.phd <- res.vec.sphd <- numeric(dd)
for (i in 1:(dd))
{
    res.vec.sir.tmp <- res.vec.phd.tmp <- res.vec.sphd.tmp <- numeric(dd)
    for (j in 1:dd)
    {
        res.vec.sir.tmp[j]  <- cor.directions(beta.hat$beta.hat[,i],     beta[,j], x.test )
        res.vec.phd.tmp[j]  <- cor.directions(beta.hat.phd$beta.hat[,i], beta[,j], x.test )
        res.vec.sphd.tmp[j] <- cor.directions(beta.hat.sphd$beta[,i],    beta[,j], x.test )
    }
    res.vec.sir[i]  <- max(res.vec.sir.tmp)
    res.vec.phd[i]  <- max(res.vec.phd.tmp)
    res.vec.sphd[i] <- max(res.vec.sphd.tmp)
}
mean(res.vec.sir)
mean(res.vec.phd)
mean(res.vec.sphd)


t(beta.hat.sphd$beta) %*% beta.hat.sphd$cov %*% beta.hat.sphd$beta
t(beta.hat.sphd$beta.init) %*% beta.hat.sphd$cov %*% beta.hat.sphd$beta.init
t(beta.hat.sphd2$beta) %*% beta.hat.sphd2$cov %*% beta.hat.sphd2$beta


plot(beta.hat.sphd2$objective)

plot(x = (x %*% beta.hat$beta.hat)[,1], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,2], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,3], y = y)

plot(x = (x %*% beta.hat.phd$beta.hat)[,1], y = y)
plot(x = (x %*% beta.hat.phd$beta.hat)[,2], y = y)
plot(x = (x %*% beta.hat.phd$beta.hat)[,3], y = y)

plot(x = (x %*% beta.hat.sphd$beta)[,1], y = y)
plot(x = (x %*% beta.hat.sphd$beta)[,2], y = y)
plot(x = (x %*% beta.hat.sphd$beta)[,3], y = y)


plot(x = (x %*% beta)[,1], y = y)
plot(x = (x %*% beta)[,2], y = y)
plot(x = (x %*% beta)[,3], y = y)

library(dr)


df <- data.frame(y = y, x)

dr.mod <- dr(y ~ ., data = df, numdir = 4, nslices = 20, method = "sir")
dr.phd <- dr(y ~ ., data = df, numdir = 5, method = "phdy")

apply(x, 2, function(xx) cor(xx, y))
apply(x %*% dr.mod$evectors, 2, function(xx) cor(xx, y))
apply(x %*% beta.hat$beta.hat, 2, function(xx) cor(xx, y))


apply(x %*% dr.phd$evectors,       2, function(xx) cor(xx, y))[1:ncol(beta.hat.phd$beta.hat)]
apply(x %*% beta.hat.phd$beta.hat, 2, function(xx) cor(xx, y))
apply(x %*% beta.hat.sphd$beta,    2, function(xx) cor(xx, y))
apply(x %*% beta.hat.sphd2$beta,   2, function(xx) cor(xx, y))

apply(x.test %*% dr.phd$evectors,       2, function(xx) round(cor(xx, y.test), 3) )[1:ncol(beta.hat.phd$beta.hat)]
apply(x.test %*% beta.hat.phd$beta.hat, 2, function(xx) round(cor(xx, y.test), 3))
apply(x.test %*% beta.hat.sphd$beta,    2, function(xx) round(cor(xx, y.test), 3))
apply(x.test %*% beta.hat.sphd2$beta,   2, function(xx) round(cor(xx, y.test), 3))
apply(x.test %*% beta,   2, function(xx) round(cor(xx, y.test), 3))



system.time(bw.obj <- npregbw(xdat = x %*% beta.hat.sphd$beta[,1:3],
                              ydat = y, bwmethod = "cv.ls",
                              bandwidth.compute = FALSE,
                              bws = h))

npr.sphd <- npreg(bws = bw.obj, txdat = x %*% beta.hat.sphd$beta[,1:3],
                  tydat = y, exdat = x.test %*% beta.hat.sphd$beta[,1:3])


h <- seq(1, 15, length.out = 15)
system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = x %*% beta.hat.sphd$beta[,1:3], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.sphd <- smooth.lf(x = x %*% beta.hat.sphd$beta[,1:3], y = y, alpha = best.h,
                                     xev = x.test %*% beta.hat.sphd$beta[,1:3]))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = x %*% beta.hat.sphd2$beta[,1:3], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.sphd2 <- smooth.lf(x = x %*% beta.hat.sphd2$beta[,1:3], y = y, alpha = best.h,
                                     xev = x.test %*% beta.hat.sphd2$beta[,1:3]))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = x %*% beta.hat.phd$beta.hat[,1:3], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.phd <- smooth.lf(x = x %*% beta.hat.phd$beta.hat[,1:3], y = y, alpha = best.h,
                                      xev = x.test %*% beta.hat.phd$beta.hat[,1:3]))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = x %*% beta.hat$beta.hat[,1:3], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.sir <- smooth.lf(x = x %*% beta.hat$beta.hat[,1:3], y = y, alpha = best.h,
                                    xev = x.test %*% beta.hat$beta.hat[,1:3]))

mean((y.test - locfit.sir$y) ^ 2)
mean((y.test - locfit.phd$y) ^ 2)
mean((y.test - locfit.sphd$y) ^ 2)
mean((y.test - locfit.sphd2$y) ^ 2)

round(beta.hat.phd$beta.hat, 4)[1:7,1:5]
round(dr.phd$evectors, 4)[1:7,1:5]
round(beta.hat.sphd$beta, 4)[1:7,]


dr.phd$M[1:5, 1:5]
beta.hat.phd$M[1:5, 1:5]

plot(abs(apply(x %*% beta.hat$beta.hat, 2, function(xx) cor(xx, y))))
plot(abs(apply(x %*% dr.mod$evectors, 2, function(xx) cor(xx, y))))

plot(abs(apply(x %*% beta.hat.phd$beta.hat, 2, function(xx) cor(xx, y))))
plot(abs(apply(x %*% dr.phd$evectors, 2, function(xx) cor(xx, y))))

summary(lm1 <- lm(y ~ x %*% dr.mod$evectors[,1:5] + I((x %*% dr.mod$evectors[,1:5])^2 )
                  + I((x %*% dr.mod$evectors[,1:2])^3 )) )

summary(lm2 <- lm(y ~ x %*% dr.mod$evectors[,1:5] + I((x %*% dr.mod$evectors[,1:5])^2 )
                  + I((x %*% dr.mod$evectors[,1:3])^3)
                  + I((x %*% dr.mod$evectors[,1:3])^4)) )


summary(lm1p <- lm(y ~ x %*% dr.phd$evectors[,1:5] + I((x %*% dr.phd$evectors[,1:5])^2 )
                  + I((x %*% dr.phd$evectors[,1:2])^3 )) )

summary(lm2p <- lm(y ~ x %*% dr.phd$evectors[,1:5] + I((x %*% dr.phd$evectors[,1:5])^2 )
                  + I((x %*% dr.phd$evectors[,1:3])^3)
                  + I((x %*% dr.phd$evectors[,1:3])^4)) )

anova(lm1, lm2)
anova(lm1p, lm2p)

summary(lm1a <- lm(y ~ x %*% beta.hat$beta.hat[,1:5] + I((x %*% beta.hat$beta.hat[,1:5])^2 )
                  + I((x %*% beta.hat$beta.hat[,1:5])^3 )) )

summary(lm2a <- lm(y ~ x %*% beta.hat$beta.hat[,1:5] + I((x %*% beta.hat$beta.hat[,1:5])^2 )
                  + I((x %*% beta.hat$beta.hat[,1:3])^3)
                  + I((x %*% beta.hat$beta.hat[,1:3])^4)) )

summary(lm1ap <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:3] + I((x %*% beta.hat.phd$beta.hat[,1:3])^2 )
                   + I((x %*% beta.hat.phd$beta.hat[,1:3])^3 )) )


summary(lm1ap1 <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:3]) )

summary(lm1ap2 <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:3] + I((x %*% beta.hat.phd$beta.hat[,1:3])^2 )
                     ) )

summary(lm2ap <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:3] + I((x %*% beta.hat.phd$beta.hat[,1:3])^2 )
                   + I((x %*% beta.hat.phd$beta.hat[,1:3])^3)
                   + I((x %*% beta.hat.phd$beta.hat[,1:3])^4)) )


summary(lm1asp1 <- lm(y ~ x %*% beta.hat.sphd$beta[,1:3] ) )

summary(lm1asp2 <- lm(y ~ x %*% beta.hat.sphd$beta[,1:3] + I((x %*% beta.hat.sphd$beta[,1:3])^2 )
                     ) )

summary(lm1asp <- lm(y ~ x %*% beta.hat.sphd$beta[,1:3] + I((x %*% beta.hat.sphd$beta[,1:3])^2 )
                    + I((x %*% beta.hat.sphd$beta[,1:3])^3 )) )

summary(lm2asp <- lm(y ~ x %*% beta.hat.sphd$beta[,1:3] + I((x %*% beta.hat.sphd$beta[,1:3])^2 )
                    + I((x %*% beta.hat.sphd$beta[,1:3])^3)
                    + I((x %*% beta.hat.sphd$beta[,1:3])^4)) )



summary(lm1asp <- lm(y ~ x %*% beta.hat.sphd2$beta[,1:3] + I((x %*% beta.hat.sphd2$beta[,1:3])^2 )
                     + I((x %*% beta.hat.sphd2$beta[,1:3])^3 )) )
anova(lm1a, lm2a)


pred1sphd1 <- cbind(1, x.test %*% beta.hat.sphd$beta[,1:3]) %*% coef(lm1asp1)

pred1phd1 <- cbind(1, x.test %*% beta.hat.phd$beta.hat[,1:3]) %*% coef(lm1ap1)

pred1sphd2 <- cbind(1, x.test %*% beta.hat.sphd$beta[,1:3],
                   (x.test %*% beta.hat.sphd$beta[,1:3])^2) %*% coef(lm1asp2)

pred1phd2 <- cbind(1, x.test %*% beta.hat.phd$beta.hat[,1:3],
                  (x.test %*% beta.hat.phd$beta.hat[,1:3])^2) %*% coef(lm1ap2)

pred1sphd <- cbind(1, x.test %*% beta.hat.sphd$beta[,1:3],
                   (x.test %*% beta.hat.sphd$beta[,1:3])^2,
                   (x.test %*% beta.hat.sphd$beta[,1:3])^3) %*% coef(lm1asp)

pred1phd <- cbind(1, x.test %*% beta.hat.phd$beta.hat[,1:3],
                  (x.test %*% beta.hat.phd$beta.hat[,1:3])^2,
                  (x.test %*% beta.hat.phd$beta.hat[,1:3])^3) %*% coef(lm1ap)

mean((y.test - pred1sphd) ^ 2)
mean((y.test - pred1phd) ^ 2)
mean((y.test - pred1sphd2) ^ 2)
mean((y.test - pred1phd2) ^ 2)
mean((y.test - pred1sphd1) ^ 2)
mean((y.test - pred1phd1) ^ 2)

plot(x = (x %*% dr.mod$evectors[,1:5])[,1], y = y)
plot(x = (x %*% dr.mod$evectors[,1:5])[,2], y = y)
plot(x = (x %*% dr.mod$evectors[,1:5])[,3], y = y)
plot(x = (x %*% dr.mod$evectors[,1:5])[,4], y = y)
plot(x = (x %*% dr.mod$evectors[,1:5])[,5], y = y)

dr.mod$evectors[,1:3]
beta.hat$eta.hat

beta.hat2 <- sir(x, y, d = 3, slice.ind = dr.mod$slice.info$slice.indicator)

dr.mod$evectors[,1:3]
beta.hat2$beta.hat




set.seed(12345)
nobs <- 200
nobs.test <- 1e4
nvars <- 40
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

y.true.a <- sin(apply(x.list[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
y.true.b <- cos((x.list[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
y.true.b <- sin(apply(x.list[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
y.true.ab <- (apply(x.list[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))

y.true <- c(y.true.a, y.true.b, y.true.ab)
y <- y.true + rnorm(nobs, sd = 1)


y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))


y.true.test <- c(y.true.a, y.true.b, y.true.ab)
y.test <- y.true.test + rnorm(nobs.test, sd = 1)

var(y.true) / 1 ^ 2

ordy <- order(y)

hier.sdr     <- hier.sir(x.list, y,  d = 1)
hier.sdr.phd <- hier.phd(x.list, y,  d = 1)
sdr.sir      <- sir(as.matrix(x), y, d = 1 * 3, h = 20L)
sdr.phd      <- phd(as.matrix(x), y, d = 1 * 3)


directions.sir      <- as.matrix(x %*% Re(hier.sdr$beta.hat))
directions.sir.test <- as.matrix(x.test %*% Re(hier.sdr$beta.hat))
directions.phd      <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
directions.phd.test <- as.matrix(x.test %*% Re(hier.sdr.phd$beta.hat))
dir.phd.non.hier      <- as.matrix(x %*% Re(sdr.phd$beta.hat))
dir.phd.non.hier.test <- as.matrix(x.test %*% Re(sdr.phd$beta.hat))
dir.sir.non.hier      <- as.matrix(x %*% Re(sdr.sir$beta.hat))
dir.sir.non.hier.test <- as.matrix(x.test %*% Re(sdr.sir$beta.hat))

library(randomForest)

rf.all          <- randomForest(x = as.matrix(x), y = y, ntree = 1000)
rf.phd.hier     <- randomForest(x = directions.phd, y = y, ntree = 1000)
rf.sir.hier     <- randomForest(x = directions.sir, y = y, ntree = 1000)
rf.sir.non.hier <- randomForest(x = dir.sir.non.hier, y = y, ntree = 1000)

summary(lmsir1  <- lm(y ~ directions.sir ))
summary(lmsir   <- lm(y ~ directions.sir + I(directions.sir ^ 2)))
summary(lmphd1  <- lm(y ~ directions.phd ))
summary(lmphd   <- lm(y ~ directions.phd + I(directions.phd ^ 2)))
summary(lmphdnh <- lm(y ~ dir.phd.non.hier + I(dir.phd.non.hier ^ 2)))
summary(lmshdnh <- lm(y ~ dir.sir.non.hier + I(dir.sir.non.hier ^ 2)))

summary(lmsir3   <- lm(y ~ directions.sir + I(directions.sir ^ 2) + I(directions.sir ^ 3)))
summary(lmphd3   <- lm(y ~ directions.phd + I(directions.phd ^ 2) + I(directions.phd ^ 3)))
summary(lmphdnh3 <- lm(y ~ dir.phd.non.hier + I(dir.phd.non.hier ^ 2) + I(dir.phd.non.hier ^ 3)))
summary(lmshdnh3 <- lm(y ~ dir.sir.non.hier + I(dir.sir.non.hier ^ 2) + I(dir.sir.non.hier ^ 3)))

preds.sir1 <- cbind(1, directions.sir.test) %*% coef(lmsir1)
preds.sir  <- cbind(1, directions.sir.test, directions.sir.test ^ 2) %*% coef(lmsir)
preds.sir3 <- cbind(1, directions.sir.test, directions.sir.test ^ 2, directions.sir.test ^ 3) %*% coef(lmsir3)
preds.phd1 <- cbind(1, directions.phd.test) %*% coef(lmphd1)
preds.phd  <- cbind(1, directions.phd.test, directions.phd.test ^ 2) %*% coef(lmphd)
preds.phd3 <- cbind(1, directions.phd.test, directions.phd.test ^ 2, directions.phd.test ^ 3) %*% coef(lmphd3)
preds.phd.nh  <- cbind(1, dir.phd.non.hier.test, dir.phd.non.hier.test ^ 2) %*% coef(lmphdnh)
preds.phd.nh3 <- cbind(1, dir.phd.non.hier.test, dir.phd.non.hier.test ^ 2, dir.phd.non.hier.test ^ 3) %*% coef(lmphdnh3)
preds.sir.nh  <- cbind(1, dir.sir.non.hier.test, dir.sir.non.hier.test ^ 2) %*% coef(lmshdnh)
preds.sir.nh3 <- cbind(1, dir.sir.non.hier.test, dir.sir.non.hier.test ^ 2, dir.sir.non.hier.test ^ 3) %*% coef(lmshdnh3)


datx   <- as.matrix(cbind(x, x ^ 2))
datxte <- as.matrix(cbind(x.test, x.test ^ 2))
library(glmnet)

glmn <- cv.glmnet(y = y, x = datx)

lm.all <- lm(y ~ as.matrix(x))

lm.all.preds <- cbind(1, as.matrix(x.test) ) %*% coef(lm.all)
glmnet.preds <- predict(glmn, newx = datxte, s = "lambda.min")


rf.preds              <- predict(rf.all, as.matrix(x.test))
rf.preds.phd.hier     <- predict(rf.phd.hier, directions.phd.test)
rf.preds.sir.hier     <- predict(rf.sir.hier, directions.sir.test)
rf.preds.sir.non.hier <- predict(rf.sir.non.hier, dir.sir.non.hier.test)


##oracle
1 - mean((y.test - y.true.test) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#lm
1 - mean((y.test - lm.all.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - rf.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - glmnet.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

# sir
#1 - mean((y.test - preds.sir1) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.sir) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.sir3) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - rf.preds.sir.hier) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
# phd
#1 - mean((y.test - preds.phd1) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.phd) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.phd3) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - rf.preds.phd.hier) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
# non hierarchical versions
1 - mean((y.test - preds.sir.nh) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.sir.nh3) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - rf.preds.sir.non.hier) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
#
1 - mean((y.test - preds.phd.nh) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - preds.phd.nh3) ^ 2) / mean((y.test - mean(y.test)) ^ 2)



h <- seq(1, 15, length.out = 15)
system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir, y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.hier.sir <- smooth.lf(x = directions.sir, y = y, alpha = best.h,
                                     xev = directions.sir.test))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir[1:nobs, 1:1], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.hier.sir.a <- smooth.lf(x = directions.sir[1:nobs, 1:1], y = y, alpha = best.h,
                                         xev = directions.sir.test[1:nobs.test, 1:1]))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir[(nobs + 1):(2*nobs), 2:2], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.hier.sir.b <- smooth.lf(x = directions.sir[(nobs + 1):(2*nobs), 2:2], y = y, alpha = best.h,
                                           xev = directions.sir.test[(nobs.test + 1):(2*nobs.test), 2:2]))

system.time(gcv.vals   <- sapply(h, function(hv) gcv(x = directions.sir[(2 * nobs + 1):(3*nobs), ], y = y, alpha = hv)[4]))
best.h     <- h[which.min(gcv.vals)]
system.time(locfit.hier.sir.ab <- smooth.lf(x = directions.sir[(2 * nobs + 1):(3*nobs), ], y = y, alpha = best.h,
                                           xev = directions.sir.test[(2 * nobs.test + 1):(3*nobs.test), ]))

hier.sir.preds <- c(locfit.hier.sir.a$y, locfit.hier.sir.b$y, locfit.hier.sir.ab$y)

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

1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)









nsims <- 100
set.seed(123)
nobs.vec <- c(250, 500, 1000, 2000)
nobs.test <- 1e4
nvars <- 40

sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.list <- rep(list(sim.res), length(nobs.vec))

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

        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = 2)


        y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list.test[[1]] %*% beta.a, 1, sum))
        y.true.b <- +0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list.test[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))^2

        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = 2)

        var(y.true) / 2 ^ 2

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

round(colMeans(sim.res.list[[2]]), 4)
round(apply(sim.res.list[[2]], 2, sd), 4)



library(reshape2)
library(ggplot2)
df.m <- data.frame(res, nobs = rep(nobs.vec, each = nsims * ncol(sim.res.list[[1]])))

colnames(df.m)[2:3] <- c("Method", "R2")
df.m <- df.m[which(df.m$R2 > -0.25),]

df.m$Method <- factor(df.m$Method, levels = levels(df.m$Method)[c(1:2, 7:8, 3:6)])

pdf(paste0(fig.path, "sim_1_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_boxplot(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

ggplot(data = df.m, aes(x=Method, y=R2)) + geom_violin(aes(fill=Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()


apply(sim.res, 2, median)

#############################

nsims <- 100
set.seed(123)
nobs.vec <- c(250, 500, 1000, 2000)[2]
nobs.test <- 1e4
nvars <- 40

sim.res <- array(NA, dim = c(nsims, 8))
colnames(sim.res) <- c("hier sir", "hier phd",
                       "hier sir (big)", "hier phd (big)",
                       "sir (big)", "phd (big)", "sir separate", "phd separate")

sim.res.dir <- array(NA, dim = c(nsims, 4))
colnames(sim.res.dir) <- c("hier sir", "hier phd", "sir separate", "phd separate")

sim2.res.list <- rep(list(sim.res), length(nobs.vec))
sim2.direction.res.list <- rep(list(sim.res.dir), length(nobs.vec))

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


        mult.a <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)
        mult.b <- matrix(runif(nvars ^ 2, max = 1/sqrt(nvars)), ncol = nvars)

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

        y.true <- c(y.true.a, y.true.b, y.true.ab)
        y <- y.true + rnorm(nobs, sd = 2)


        y.true.a <- sin(apply(x.list.test[[1]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[1]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.b <- cos((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2 / ((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] - 1.5) ^ 2))
        y.true.b <- sin(apply(x.list.test[[2]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[2]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- (apply(x.list.test[[3]] %*% beta.ab[,2:3,drop=FALSE], 1, sum))^2 / ((0.5 + ( apply(x.list.test[[3]] %*% beta.ab[,,drop=FALSE], 1, sum) - 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.b, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.b[,1,drop=FALSE] + 1.5) ^ 2))
        y.true.ab <- y.true.ab + sin(apply(x.list.test[[3]] %*% beta.a, 1, sum)) ^ 2 / sqrt((0.5 + (x.list.test[[3]] %*% beta.a[,1,drop=FALSE] + 1.5) ^ 2))

        y.true.a <- 2 * cos(apply(x.list.test[[1]] %*% beta.a, 1, sum))
        y.true.b <- +0.1 * ((x.list.test[[2]] %*% beta.b[,1]) ^ 2) ^ 2
        y.true.ab <- 2 * cos(apply(x.list.test[[3]] %*% beta.a, 1, sum)) +
            0.1 * ((x.list.test[[3]] %*% beta.ab[,2]) ^ 2) ^ 2 +
            0.1 * (apply(sin(x.list.test[[3]] %*% beta.ab[,3]), 1, sum))^2

        y.true.test <- c(y.true.a, y.true.b, y.true.ab)
        y.test <- y.true.test + rnorm(nobs.test, sd = 2)

        var(y.true) / 2 ^ 2

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

        hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])

        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1],    beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])

        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     mult.a %*% beta.a, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), mult.a %*% beta.a, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,1], mult.a %*% beta.a, x.list.test[[3]])

        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     mult.b %*% beta.b, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), mult.b %*% beta.b, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,2], mult.b %*% beta.b, x.list.test[[3]])

        hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     eta.ab, x.list.test[[3]])
        hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), eta.ab, x.list.test[[3]])
        phd.cor      <- phd.cor      + cor.directions(phd.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])
        sir.cor      <- sir.cor      + cor.directions(sir.3$beta.hat[1:nvars,3], eta.ab, x.list.test[[3]])

        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5

        sim2.direction.res.list[[n]][s,1] <- hier.sir.cor
        sim2.direction.res.list[[n]][s,2] <- hier.phd.cor

        sim2.direction.res.list[[n]][s,3] <- sir.cor
        sim2.direction.res.list[[n]][s,4] <- phd.cor

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

        sim2.res.list[[n]][s,1] <- 1 - mean((y.test - hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim2.res.list[[n]][s,2] <- 1 - mean((y.test - hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim2.res.list[[n]][s,3] <- 1 - mean((y.test - locfit.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim2.res.list[[n]][s,4] <- 1 - mean((y.test - locfit.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim2.res.list[[n]][s,5] <- 1 - mean((y.test - locfit.non.hier.sir$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim2.res.list[[n]][s,6] <- 1 - mean((y.test - locfit.non.hier.phd$y) ^ 2) / mean((y.test - mean(y.test)) ^ 2)

        sim2.res.list[[n]][s,7] <- 1 - mean((y.test - non.hier.sir.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)
        sim2.res.list[[n]][s,8] <- 1 - mean((y.test - non.hier.phd.preds) ^ 2) / mean((y.test - mean(y.test)) ^ 2)


        if (s %% 5 == 0) cat("sim:", s, "complete \n")
    }
}

res2 <- do.call(rbind, lapply(sim2.res.list, melt))
res.dir2 <- do.call(rbind, lapply(sim2.direction.res.list, melt))


df.m2 <- data.frame(res2, nobs = rep(nobs.vec, each = nsims * ncol(sim2.res.list[[1]])))

colnames(df.m2)[2:3] <- c("Method", "R2")
df.m2 <- df.m2[which(df.m2$R2 > -0.25),]

df.m2$Method <- factor(df.m2$Method, levels = levels(df.m2$Method)[c(1:2, 7:8, 3:6)])

ggplot(data = df.m2, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw()

pdf(paste0(fig.path, "sim_2_boxplots.pdf"), height = 8, width = 10)

ggplot(data = df.m2, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2 Test Set") + facet_wrap(~ nobs) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))

dev.off()

## directions R^2

df.m.dir2 <- data.frame(res.dir2, nobs = rep(nobs.vec, each = nsims * ncol(sim2.direction.res.list[[1]])))

colnames(df.m.dir2)[2:3] <- c("Method", "R2")
#df.m.dir2 <- df.m.dir2[which(df.m.dir2$R2 > -0.25),]

df.m.dir2$Method <- factor(df.m.dir2$Method, levels = levels(df.m.dir2$Method)[c(1:2, 7:8, 3:6)])

ggplot(data = df.m.dir2, aes(x=Method, y=R2)) + geom_boxplot(aes(fill = Method)) +
    ylab("R^2") + facet_wrap(~ nobs) + theme_bw()

#############################

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
