
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
nobs <- 50
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

beta.hat <- sir(x, y, h = 20, d = 10)

plot(x = (x %*% beta.hat$beta.hat)[,1], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,2], y = y)
plot(x = (x %*% beta.hat$beta.hat)[,3], y = y)

library(dr)


df <- data.frame(y = y, x)

dr.mod <- dr(y ~ ., data = df, numdir = 4, nslices = 20, method = "sir")

apply(x, 2, function(xx) cor(xx, y))
apply(x %*% dr.mod$evectors, 2, function(xx) cor(xx, y))
apply(x %*% beta.hat$beta.hat, 2, function(xx) cor(xx, y))



summary(lm1 <- lm(y ~ x %*% dr.mod$evectors[,1:5] + I((x %*% dr.mod$evectors[,1:5])^2 )
           + I((x %*% dr.mod$evectors[,1:2])^3 )) )

summary(lm2 <- lm(y ~ x %*% dr.mod$evectors[,1:5] + I((x %*% dr.mod$evectors[,1:5])^2 )
           + I((x %*% dr.mod$evectors[,1:3])^3)
               + I((x %*% dr.mod$evectors[,1:3])^4)) )

anova(lm1, lm2)


summary(lm1a <- lm(y ~ x %*% beta.hat$beta.hat[,1:5] + I((x %*% beta.hat$beta.hat[,1:5])^2 )
                  + I((x %*% beta.hat$beta.hat[,1:5])^3 )) )

summary(lm2a <- lm(y ~ x %*% beta.hat$beta.hat[,1:5] + I((x %*% beta.hat$beta.hat[,1:5])^2 )
                  + I((x %*% beta.hat$beta.hat[,1:3])^3)
                  + I((x %*% beta.hat$beta.hat[,1:3])^4)) )

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


