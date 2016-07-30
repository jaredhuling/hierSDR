
y <- runif(200, -pi, pi)
x <- sin(y) + rnorm(200, sd = 0.25)

nw.fit <- nwsmooth(x = y, y = x, h = 0.75)
plot(x = y, y = x)
ordy <- order(y)
lines(x = y[ordy], y = nw.fit$fitted[ordy])


x <- matrix(rnorm(200 * 10), ncol = 10)
y <- sin(apply(x, 1, function(xr) exp( (xr[1] + xr[2])^2 ))) +
    exp(apply(x, 1, function(xr) cos( (xr[3] + xr[4])^2 ))) +
    cos(apply(x, 1, function(xr) ( sum(xr[5:length(xr)])^2 ))) +
    rnorm(200, sd = 0.25)
ordy <- order(y)
nw.cov.fit <- nwsmoothcov(x = x, y = y, h = 1)

nw.cov.fit$fitted[ordy]
