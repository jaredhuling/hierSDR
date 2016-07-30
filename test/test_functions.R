
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
x <- matrix(rnorm(200 * 10), ncol = 10)
y <- sin(apply(x, 1, function(xr) exp( (xr[1] + xr[2])^2 ))) +
    exp(apply(x, 1, function(xr) cos( (xr[3] + xr[4])^2 ))) +
    cos(apply(x, 1, function(xr) ( sum(xr[5:length(xr)])^2 ))) +
    rnorm(200, sd = 0.25)
ordy <- order(y)


bandwidths <- exp(seq(log(0.25), log(30), length.out = 25))
nw.cov.fits <- vector(mode = "list", length = length(bandwidths))
for (b in 1:length(bandwidths)) nw.cov.fits[[b]] <- nwsmoothcov(x = x, y = y, h = bandwidths[b])
gcvs.cov <- sapply(nw.cov.fits, function(f) f$gcv)
bandwidths[which.min(gcvs.cov)]

nw.cov.fit <- nw.cov.fits[[which.min(gcvs.cov)]]

#    nw.cov.fit$fitted[ordy]

cov2plot <- sapply(nw.cov.fit$fitted[ordy], function(x) x[9,3])

plot(x = y[ordy], y = cov2plot, type = "b")
