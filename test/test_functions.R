
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
nobs <- 100
x <- matrix(rnorm(nobs * 40), ncol = 40)
y <- (apply(x, 1, function(xr) sum( (xr[1] + xr[2])^2 ))) +
    sum(apply(x, 1, function(xr) ( (xr[3] + xr[4])^2 ))) +
    (apply(x, 1, function(xr) ( sum(exp(xr[5:min(length(xr), 10)])) ))) +
    rnorm(nobs, sd = 0.25)
ordy <- order(y)


bandwidths <- exp(seq(log(0.25), log(30), length.out = 25))[10:25]
nw.cov.fits <- vector(mode = "list", length = length(bandwidths))
for (b in 1:length(bandwidths)) nw.cov.fits[[b]] <- nwsmoothcov(x = x, y = y, h = bandwidths[b])
gcvs.cov <- sapply(nw.cov.fits, function(f) f$gcv)
bandwidths[which.min(gcvs.cov)]

nw.cov.fit <- nw.cov.fits[[which.min(gcvs.cov)]]

#    nw.cov.fit$fitted[ordy]

cov2plot <- sapply(nw.cov.fit$fitted[ordy], function(x) x[9,3])

plot(x = y[ordy], y = cov2plot, type = "b")



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
    list(beta.hat = beta.hat, eta.hat = eta.hat)
}


phd <- function(x, y, d = 5L)
{
    cov <- cov(x)
    eig.cov <- eigen(cov)
    sqrt.inv.cov <- eig.cov$vectors %*% diag(1 / sqrt(eig.cov$values)) %*% t(eig.cov$vectors)
    x.tilde <- scale(x, scale = FALSE) %*% sqrt.inv.cov

    V.hat <- crossprod(x.tilde, drop(scale(y, scale = FALSE)) * x.tilde) / nrow(x)
    eig.V <- eigen(V.hat)
    eta.hat <- eig.V$vectors[,1:d]
    beta.hat <- t(t(eta.hat) %*% sqrt.inv.cov)
    list(beta.hat = beta.hat, eta.hat = eta.hat, M = V.hat, cov = cov, sqrt.inv.cov = sqrt.inv.cov)
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

beta.hat     <- sir(x, y, h = 20, d = 40)
beta.hat.phd <- phd(x, y, d = 40)

plot(x = (x %*% beta.hat$beta.hat)[,1], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,2], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,3], y = y)

plot(x = (x %*% beta.hat.phd$beta.hat)[,1], y = y)
plot(x = (x %*% beta.hat.phd$beta.hat)[,2], y = y)
plot(x = (x %*% beta.hat.phd$beta.hat)[,3], y = y)

library(dr)


df <- data.frame(y = y, x)

dr.mod <- dr(y ~ ., data = df, numdir = 4, nslices = 20, method = "sir")
dr.phd <- dr(y ~ ., data = df, numdir = 5, method = "phdy")

apply(x, 2, function(xx) cor(xx, y))
apply(x %*% dr.mod$evectors, 2, function(xx) cor(xx, y))
apply(x %*% beta.hat$beta.hat, 2, function(xx) cor(xx, y))


apply(x %*% dr.phd$evectors,       2, function(xx) cor(xx, y))[1:ncol(beta.hat.phd$beta.hat)]
apply(x %*% beta.hat.phd$beta.hat, 2, function(xx) cor(xx, y))


round(beta.hat.phd$beta.hat, 4)[1:7,1:5]
round(dr.phd$evectors, 4)[1:7,1:5]


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

summary(lm1ap <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:5] + I((x %*% beta.hat.phd$beta.hat[,1:5])^2 )
                   + I((x %*% beta.hat.phd$beta.hat[,1:5])^3 )) )

summary(lm2ap <- lm(y ~ x %*% beta.hat.phd$beta.hat[,1:5] + I((x %*% beta.hat.phd$beta.hat[,1:5])^2 )
                   + I((x %*% beta.hat.phd$beta.hat[,1:3])^3)
                   + I((x %*% beta.hat.phd$beta.hat[,1:3])^4)) )

anova(lm1a, lm2a)


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




set.seed(123)
nobs <- 100
x.list <- replicate(3, list(matrix(rnorm(nobs * 70), ncol = 70)))
x <- bdiag(x.list)
y <- (apply(x, 1, function(xr) sum( (xr[1] + xr[2])^2 ))) +
    sum(apply(x, 1, function(xr) ( (xr[3] + xr[4])^2 ))) +
    (apply(x, 1, function(xr) ( sum(cos(xr[5:min(length(xr), 40)]) ^ 2) ))) +
    rnorm(nobs * length(x.list), sd = 3.25)

y <- as.vector(sapply(x.list, function(xx) (apply(xx, 1, function(xr) sum( (xr[1] + xr[2])^2 ))) +
    sum(apply(xx, 1, function(xr) ( (xr[3] + xr[4])^2 ))) +
    (apply(xx, 1, function(xr) ( sum(cos(xr[5:min(length(xr), 40)]) ^ 2) ))) +
    rnorm(nobs, sd = 3.25)))

ordy <- order(y)

hier.sdr     <- hier.sir(x.list, y,  d = 2)
hier.sdr.phd <- hier.phd(x.list, y,  d = 2)
sdr.phd      <- phd(as.matrix(x), y, d = 2 * 3)

beta.c     <- t(hier.sdr$beta.hat) %*% hier.sdr$cov %*% hier.sdr$beta.hat
beta.c.phd <- t(hier.sdr.phd$beta.hat) %*% hier.sdr.phd$cov %*% hier.sdr.phd$beta.hat

round(beta.c, 4)
round(beta.c.phd, 4)


directions.sir   <- as.matrix(x %*% Re(hier.sdr$beta.hat))
directions.phd   <- as.matrix(x %*% Re(hier.sdr.phd$beta.hat))
dir.phd.non.hier <- as.matrix(x %*% Re(sdr.phd$beta.hat))

summary(lmsir   <- lm(y ~ directions.sir + I(directions.sir ^ 2)))
summary(lmphd   <- lm(y ~ directions.phd + I(directions.phd ^ 2)))
summary(lmphdnh <- lm(y ~ dir.phd.non.hier + I(dir.phd.non.hier ^ 2)))

summary(lmsir3   <- lm(y ~ directions.sir + I(directions.sir ^ 2) + I(directions.sir ^ 3)))
summary(lmphd3   <- lm(y ~ directions.phd + I(directions.phd ^ 2) + I(directions.phd ^ 3)))
summary(lmphdnh3 <- lm(y ~ dir.phd.non.hier + I(dir.phd.non.hier ^ 2) + I(dir.phd.non.hier ^ 3)))

anova(lmsir3, lmsir)
anova(lmphd3, lmphd)
anova(lmphdnh3, lmphdnh)

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
