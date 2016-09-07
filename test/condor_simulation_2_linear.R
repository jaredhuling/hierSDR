
arg <- commandArgs()      # Gets the process number
print(arg)

print(eval(parse(text=arg[[8]])))
repnum <- eval(parse(text=arg[[8]]))  # Gets a workable character string from the argument
seed.num <- repnum
print(seed.num)


## load required packages
library(hierSDR)
library(Matrix)


## set up simulation grid
nobs.vec      <- c(250, 500, 1000)
nvars         <- 50
ncats         <- 2
x.type        <- c("simple", "complex")
beta.type.vec <- c("some.zero", "some.small")[1]
model.num.vec <- c(3, 4, 5)
nobs.test     <- 10000
simtype       <- "regr"
nsims         <- 100
sd.sim        <- 1
num.models    <- 7

#    set up grid of parameters over which to simulate #
grid           <- expand.grid(nobs.vec, nvars, x.type, beta.type.vec, model.num.vec)
colnames(grid) <- c("nobs", "nvars", "x.type", "beta.type", "model.num")

n.jobs <- nrow(grid)

## get index for the correct simulation parameters
## for this specific job (job id = seed.num) # job ids start at 0
sim.idx <- (seed.num + 1 - 1) %% nrow(grid) + 1

set.seed(seed.num)

nvars     <- grid$nvars[sim.idx]
x.type    <- grid$x.type[sim.idx]
nobs      <- grid$nobs[sim.idx]
beta.type <- grid$beta.type[sim.idx]
model.num <- grid$model.num[sim.idx]

## generate data
source("condor_generate_data_linear.R")


## fit SDR models

if (model.num == 1 | model.num == 3)
{
    hier.sdr      <- hier.sir(x.list, y,  d = c(1,1,1), h = 30L)
    hier.sdr.phd  <- hier.phd(x.list, y,  d = c(1,1,1))
    sdr.sir       <- sir(as.matrix(x), y, d = 1 * 3,    h = 30L)
    sdr.phd       <- phd(as.matrix(x), y, d = 1 * 3)

    semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = c(1,1,1), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 500, maxk = 1200)
    semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = c(1,1,1), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 500, maxk = 1200)


    sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
    sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
    sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = 30L)

    phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
    phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
    phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3)

    s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)
    s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)
    s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 3, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)

} else if (model.num == 2 | model.num == 4)
{
    hier.sdr      <- hier.sir(x.list, y,  d = c(1,1,0), h = 30L)
    hier.sdr.phd  <- hier.phd(x.list, y,  d = c(1,1,0))
    sdr.sir       <- sir(as.matrix(x), y, d = 1 * 2, h = 30L)
    sdr.phd       <- phd(as.matrix(x), y, d = 1 * 2)

    semi.hier.phd  <- hier.semi.phd(x.list, drop(y), d = c(1,1,0), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 500, maxk = 1200)
    semi.hier.phd2 <- semi.phd.hier(x.list, drop(y), d = c(1,1,0), h = exp(seq(log(0.5), log(25), length.out = 25)), maxit = 500, maxk = 1200)


    sir.1 <- sir(x.list[[1]], y[1:nobs],                d = 1, h = 30L)
    sir.2 <- sir(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = 30L)
    sir.3 <- sir(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2, h = 30L)

    phd.1 <- phd(x.list[[1]], y[1:nobs],                d = 1)
    phd.2 <- phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1)
    phd.3 <- phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2)

    s.phd.1 <- semi.phd(x.list[[1]], y[1:nobs],                d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)
    s.phd.2 <- semi.phd(x.list[[2]], y[(1 + nobs):(2*nobs)],   d = 1, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)
    s.phd.3 <- semi.phd(x.list[[3]], y[(1 + 2*nobs):(3*nobs)], d = 2, h = exp(seq(log(0.25), log(25), length.out = 25)), maxit = 500, maxk = 450)

}




if ( !(model.num %in% eta.zero.models) )
{
    ## beta A
    hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor   <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
    phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
    sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
    s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


    ## beta B
    hier.sir.cor <- hier.sir.cor + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),     beta.b, x.list.test[[2]])
    hier.phd.cor <- hier.phd.cor + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    semi.hier.phd.cor2  <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
    sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
    s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

    ## beta AB -> A
    hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
    hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
    sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
    s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

    ## beta AB -> B
    hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
    hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
    sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
    s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

    ## beta AB -> eta AB
    hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),3]),     beta.ab[,idx], x.list.test[[3]])))
    hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
    phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
    sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
    s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

    hier.sir.cor <- hier.sir.cor / 5
    hier.phd.cor <- hier.phd.cor / 5
    semi.hier.phd.cor <- semi.hier.phd.cor / 5
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 5
    semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
    phd.cor      <- phd.cor / 5
    sir.cor      <- sir.cor / 5
    s.phd.cor    <- s.phd.cor / 5

    ######## angles

    hier.sir.angle <- (1/5) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 3 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    hier.phd.angle <- (1/5) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 3 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    semi.hier.phd.angle <- (1/5) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 3 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    #sim.subsp.angle.res.list4[[n]][s,8] <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]),beta.b) )
    #                                                + 3 * angles(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),]), beta.ab))

    semi.hier.phd.angle2 <- (1/5) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 3 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

    phd.angle <- (1/5) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                    + 3 * angles(phd.3$beta.hat, beta.ab))

    sir.angle <- (1/5) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                    + 3 * angles(sir.3$beta.hat, beta.ab))

    semi.phd.angle <- (1/5) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                    + 3 * angles(s.phd.3$beta, beta.ab))


} else
{

    ## beta A
    hier.sir.cor <- cor.directions(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    hier.phd.cor <- cor.directions(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor   <- cor.directions(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a, x.list.test[[1]])
    #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
    semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
    phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
    sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
    s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


    ## beta B
    hier.sir.cor        <- hier.sir.cor        + cor.directions(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),                beta.b, x.list.test[[2]])
    hier.phd.cor        <- hier.phd.cor        + cor.directions(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),            beta.b, x.list.test[[2]])
    semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),             beta.b, x.list.test[[2]])
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]),        beta.b, x.list.test[[2]])
    semi.hier.phd.cor2  <- semi.hier.phd.cor2  + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),              beta.b, x.list.test[[2]])
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
    phd.cor             <- phd.cor             + cor.directions(phd.2$beta.hat[1:nvars,1],                                     beta.b, x.list.test[[2]])
    sir.cor             <- sir.cor             + cor.directions(sir.2$beta.hat[1:nvars,1],                                     beta.b, x.list.test[[2]])
    s.phd.cor           <- s.phd.cor           + cor.directions(s.phd.2$beta[1:nvars,1],                                       beta.b, x.list.test[[2]])

    ## beta AB -> A
    hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),1]),     beta.ab[,idx], x.list.test[[3]])))
    hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
    phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
    sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
    s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

    ## beta AB -> B
    hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),2]),     beta.ab[,idx], x.list.test[[3]])))
    hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
    phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
    sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
    s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

    hier.sir.cor <- hier.sir.cor / 4
    hier.phd.cor <- hier.phd.cor / 4
    semi.hier.phd.cor <- semi.hier.phd.cor / 4
    #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 4
    semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 4
    semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 4
    phd.cor      <- phd.cor / 4
    sir.cor      <- sir.cor / 4
    s.phd.cor    <- s.phd.cor / 4


    ##### angles
    hier.sir.angle <- (1/4) * ((angles(Re(hier.sdr$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 2 * angles(Re(hier.sdr$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    hier.phd.angle <- (1/4) * ((angles(Re(hier.sdr.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(hier.sdr.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 2 * angles(Re(hier.sdr.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    semi.hier.phd.angle <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 2 * angles(Re(semi.hier.phd$beta.hat[(2*nvars+1):(3 * nvars),]), beta.ab))

    #sim.subsp.angle.res.list3[[n]][s,8] <- (1/4) * ((angles(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]),beta.b) )
    #                                                + 2 * angles(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),]), beta.ab))

    semi.hier.phd.angle2 <- (1/4) * ((angles(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a) + angles(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]),beta.b) )
                                                    + 2 * angles(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),]), beta.ab))

    phd.angle <- (1/4) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                    + 2 * angles(phd.3$beta.hat, beta.ab))

    sir.angle <- (1/4) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                    + 2 * angles(sir.3$beta.hat, beta.ab))

    semi.phd.angle <- (1/4) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                    + 2 * angles(s.phd.3$beta, beta.ab))
}


result.mat.cor <- result.mat.angle <- data.frame(grid[sim.idx,], array(NA, dim = c(1, num.models + 1)))

colnames(result.mat.cor) <-
    colnames(result.mat.angle) <-
    c(colnames(grid), "SNR",
      "SIR", "Hier SIR", "PHD", "Hier PHD",
      "Semi PHD", "Semi Hier PHD", "Semi Hier PHD 2")


result.mat.cor[1,(ncol(grid) + 1):ncol(result.mat.cor)] <- c(snr,
                                                             sir.cor, hier.sir.cor,
                                                             phd.cor, hier.phd.cor,
                                                             s.phd.cor,
                                                             semi.hier.phd.cor,
                                                             semi.hier.phd.cor2)


result.mat.angle[1,(ncol(grid) + 1):ncol(result.mat.angle)] <- c(snr,
                                                                 sir.angle, hier.sir.angle,
                                                                 phd.angle, hier.phd.angle,
                                                                 semi.phd.angle,
                                                                 semi.hier.phd.angle,
                                                                 semi.hier.phd.angle2)


print("semi phd angle:")
print(semi.phd.angle)
print("semi hier phd angle:")
print(semi.hier.phd.angle)

print(result.mat.angle)

## save results
write.csv(result.mat.cor,   paste(seed.num, "_results_correlation.csv", sep = ""))
write.csv(result.mat.angle, paste(seed.num, "_results_angle.csv",       sep = ""))


