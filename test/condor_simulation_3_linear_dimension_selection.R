
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
nobs.vec      <- c(250, 500)
nvars         <- 10
ncats         <- 2
x.type        <- c("simple", "complex")[2]
beta.type.vec <- c("some.zero", "some.small")[1]
model.num.vec <- c(2, 5)
nobs.test     <- 10000
simtype       <- "regr"
nsims         <- 100
sd.sim        <- 1
num.models    <- 6
max.d         <- c(3, 3, 2)

d.mat <- data.matrix(expand.grid(1:max.d[1], 1:max.d[2], 0:max.d[3]))
dimnames(d.mat) <- NULL

d.nums <- 1:nrow(d.mat)

#    set up grid of parameters over which to simulate #
grid           <- expand.grid(nobs.vec, nvars, x.type, beta.type.vec, model.num.vec, d.nums)
colnames(grid) <- c("nobs", "nvars", "x.type", "beta.type", "model.num", "d.num")

n.jobs <- nrow(grid)

n.d <- nrow(d.mat)
n.unique.grids <- nrow(grid) / n.d

## get index for the correct simulation parameters
## for this specific job (job id = seed.num) # job ids start at 0
sim.idx <- (seed.num + 1 - 1) %% nrow(grid) + 1


rep.num <- (seed.num + 1 - sim.idx) / nrow(grid) + 1

grid.idx  <- c(sim.idx - 1) %% n.unique.grids + 1

set.seed(grid.idx * 1000 + rep.num)

nvars     <- grid$nvars[sim.idx]
x.type    <- grid$x.type[sim.idx]
nobs      <- grid$nobs[sim.idx]
beta.type <- grid$beta.type[sim.idx]
model.num <- grid$model.num[sim.idx]
d.cur     <- as.vector(d.mat[grid$d.num[sim.idx], ])

print(grid[sim.idx,])

Proj <- function(b) b %*% solve(crossprod(b), t(b))

proj.norm <- function(b1, b2)
{
    norm(Proj(b1) - Proj(b2), type = "F")
}

## generate data
#  source("C:/Users/Jared/Documents/GitHub/hierSDR/test/condor_generate_data_linear.R")
source("condor_generate_data_linear_2.R")

## fit SDR models

semi.hier.phd  <- semi.phd.hier.separate.dim.cv.k(x.list, drop(y),
                                                  d = d.cur, k = 5,
                                                  h = exp(seq(log(0.5), log(25), length.out = 25)),
                                                  maxit = 50, maxk = 2500)
#

gcvs  <- semi.hier.phd$gcvs


if ( !(model.num %in% eta.zero.models) )
{
    d.true <- c(1, 1, 1)

} else
{
    d.true <- c(1, 1, 0)
}



if (FALSE)
{
    if ( !(model.num %in% eta.zero.models) )
    {
        ## beta A
        hier.sir.cor <- cor.directions(hier.sir.betaAB[,1], beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(hier.phd.betaAB[,1], beta.a, x.list.test[[1]])
        semi.hier.phd.cor   <- cor.directions(semi.hier.phd$beta[[1]], beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])

        cor.directions(semi.hier.phd$beta[[1]], beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(hier.sir.betaAB[,2],     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(hier.phd.betaAB[,2], beta.b, x.list.test[[2]])
        semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(semi.hier.phd$beta[[2]], beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor2  <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        cor.directions(semi.hier.phd$beta[[2]], beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(hier.sir.betaAB[,1],     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(hier.phd.betaAB[,1], beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,1], beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(hier.sir.betaAB[,2],     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(hier.phd.betaAB[,2], beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,2], beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))

        max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,2], beta.ab[,idx], x.list.test[[3]])))


        ## beta AB -> eta AB
        hier.sir.cor <- hier.sir.cor + max(sapply(1:3, function(idx) cor.directions(hier.sir.betaAB[,3],     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(hier.phd.betaAB[,3], beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,3], beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:3, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),3]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:3, function(idx) cor.directions(phd.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:3, function(idx) cor.directions(sir.3$beta.hat[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:3, function(idx) cor.directions(s.phd.3$beta[1:nvars,3], beta.ab[,idx], x.list.test[[3]])))

        max(sapply(1:3, function(idx) cor.directions(semi.hier.phd$beta[[3]][,3], beta.ab[,idx], x.list.test[[3]])))


        hier.sir.cor <- hier.sir.cor / 5
        hier.phd.cor <- hier.phd.cor / 5
        semi.hier.phd.cor <- semi.hier.phd.cor / 5
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 5
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
        phd.cor      <- phd.cor / 5
        sir.cor      <- sir.cor / 5
        s.phd.cor    <- s.phd.cor / 5

        ######## angles

        hier.sir.angle <- (1/5) * ((angles(hier.sir.betaAB[,1], beta.a) + angles(hier.sir.betaAB[,2], beta.b) )
                                                        + 3 * angles(hier.sir.betaAB, beta.ab))

        hier.phd.angle <- (1/5) * ((angles(hier.phd.betaAB[,1], beta.a) + angles(hier.phd.betaAB[,2], beta.b) )
                                                        + 3 * angles(hier.phd.betaAB, beta.ab))

        semi.hier.phd.angle <- (1/5) * ((angles(semi.hier.phd$beta[[1]], beta.a) + angles(semi.hier.phd$beta[[2]],beta.b) )
                                                        + 3 * angles(semi.hier.phd$beta[[3]], beta.ab))

        phd.angle <- (1/5) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 3 * angles(phd.3$beta.hat, beta.ab))

        sir.angle <- (1/5) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                                                        + 3 * angles(sir.3$beta.hat, beta.ab))

        semi.phd.angle <- (1/5) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                                        + 3 * angles(s.phd.3$beta, beta.ab))


        (1/5) * (( angles(semi.hier.phd$beta[[1]], beta.a) + angles(semi.hier.phd$beta[[2]],beta.b) )
                 + 3 * angles(semi.hier.phd$beta[[3]], beta.ab))

        ######## projection difference norm

        hier.sir.proj <- (1/5) * ((proj.norm(hier.sir.betaAB[,1], beta.a) + proj.norm(hier.sir.betaAB[,2], beta.b) )
                                   + 3 * proj.norm(hier.sir.betaAB, beta.ab))

        hier.phd.proj <- (1/5) * ((proj.norm(hier.phd.betaAB[,1], beta.a) + proj.norm(hier.phd.betaAB[,2], beta.b) )
                                   + 3 * proj.norm(hier.phd.betaAB, beta.ab))

        semi.hier.phd.proj <- (1/5) * ((proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
                                        + 3 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

        phd.proj <- (1/5) * ((proj.norm(phd.1$beta.hat[1:nvars,1], beta.a) + proj.norm(phd.2$beta.hat[1:nvars,1],beta.b) )
                              + 3 * proj.norm(phd.3$beta.hat, beta.ab))

        sir.proj <- (1/5) * (( proj.norm(sir.1$beta.hat[1:nvars,1], beta.a) + proj.norm(sir.2$beta.hat[1:nvars,1],beta.b) )
                              + 3 * proj.norm(sir.3$beta.hat, beta.ab))

        semi.phd.proj <- (1/5) * ((proj.norm(s.phd.1$beta[1:nvars,1], beta.a)+ proj.norm(s.phd.2$beta[1:nvars,1],beta.b) )
                                   + 3 * proj.norm(s.phd.3$beta, beta.ab))

        (1/5) * (( proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
                 + 3 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

    } else
    {

        ## beta A
        hier.sir.cor <- cor.directions(hier.sir.betaAB[,1], beta.a, x.list.test[[1]])
        hier.phd.cor <- cor.directions(hier.phd.betaAB[,1], beta.a, x.list.test[[1]])
        semi.hier.phd.cor   <- cor.directions(semi.hier.phd$beta[[1]], beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.W <- cor.directions(Re(semi.hier.phd$beta.hat.W[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor2  <- cor.directions(Re(semi.hier.phd2$beta[1:nvars,1]), beta.a, x.list.test[[1]])
        #semi.hier.phd.cor.u <- cor.directions(Re(semi.hier.phd$beta.unconstrained[1:nvars,1]), beta.a, x.list.test[[1]])
        phd.cor      <- cor.directions(phd.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        sir.cor      <- cor.directions(sir.1$beta.hat[1:nvars,1], beta.a, x.list.test[[1]])
        s.phd.cor    <- cor.directions(s.phd.1$beta[1:nvars,1],   beta.a, x.list.test[[1]])


        ## beta B
        hier.sir.cor <- hier.sir.cor + cor.directions(hier.sir.betaAB[,2],     beta.b, x.list.test[[2]])
        hier.phd.cor <- hier.phd.cor + cor.directions(hier.phd.betaAB[,2], beta.b, x.list.test[[2]])
        semi.hier.phd.cor   <- semi.hier.phd.cor + cor.directions(semi.hier.phd$beta[[2]], beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + cor.directions(Re(semi.hier.phd$beta.hat.W[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor2  <- semi.hier.phd.cor2 + cor.directions(Re(semi.hier.phd2$beta[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + cor.directions(Re(semi.hier.phd$beta.unconstrained[(nvars+1):(2 * nvars),2]), beta.b, x.list.test[[2]])
        phd.cor      <- phd.cor      + cor.directions(phd.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        sir.cor      <- sir.cor      + cor.directions(sir.2$beta.hat[1:nvars,1], beta.b, x.list.test[[2]])
        s.phd.cor    <- s.phd.cor    + cor.directions(s.phd.2$beta[1:nvars,1],   beta.b, x.list.test[[2]])

        ## beta AB -> A
        hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(hier.sir.betaAB[,1],     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(hier.phd.betaAB[,1], beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(semi.hier.phd$beta[[3]][,1], beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),1]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,1], beta.ab[,idx], x.list.test[[3]])))

        ## beta AB -> B
        hier.sir.cor <- hier.sir.cor + max(sapply(1:2, function(idx) cor.directions(hier.sir.betaAB[,2],     beta.ab[,idx], x.list.test[[3]])))
        hier.phd.cor <- hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(hier.phd.betaAB[,2], beta.ab[,idx], x.list.test[[3]])))
        semi.hier.phd.cor   <- semi.hier.phd.cor + max(sapply(1:2, function(idx) cor.directions(semi.hier.phd$beta[[3]][,2], beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.hat.W[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor2  <- semi.hier.phd.cor2 + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd2$beta[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u + max(sapply(1:2, function(idx) cor.directions(Re(semi.hier.phd$beta.unconstrained[(2*nvars+1):(3 * nvars),2]), beta.ab[,idx], x.list.test[[3]])))
        phd.cor      <- phd.cor      + max(sapply(1:2, function(idx) cor.directions(phd.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        sir.cor      <- sir.cor      + max(sapply(1:2, function(idx) cor.directions(sir.3$beta.hat[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))
        s.phd.cor    <- s.phd.cor    + max(sapply(1:2, function(idx) cor.directions(s.phd.3$beta[1:nvars,2], beta.ab[,idx], x.list.test[[3]])))


        hier.sir.cor <- hier.sir.cor / 4
        hier.phd.cor <- hier.phd.cor / 4
        semi.hier.phd.cor <- semi.hier.phd.cor / 4
        #semi.hier.phd.cor.W <- semi.hier.phd.cor.W / 5
        #semi.hier.phd.cor2 <- semi.hier.phd.cor2 / 5
        #semi.hier.phd.cor.u <- semi.hier.phd.cor.u / 5
        phd.cor      <- phd.cor / 4
        sir.cor      <- sir.cor / 4
        s.phd.cor    <- s.phd.cor / 4

        ######## angles

        hier.sir.angle <- (1/4) * ((angles(hier.sir.betaAB[,1], beta.a) + angles(hier.sir.betaAB[,2], beta.b) )
                                   + 2 * angles(hier.sir.betaAB, beta.ab))

        hier.phd.angle <- (1/4) * ((angles(hier.phd.betaAB[,1], beta.a) + angles(hier.phd.betaAB[,2], beta.b) )
                                   + 2 * angles(hier.phd.betaAB, beta.ab))

        semi.hier.phd.angle <- (1/4) * ((angles(semi.hier.phd$beta[[1]], beta.a) + angles(semi.hier.phd$beta[[2]],beta.b) )
                                        + 2 * angles(semi.hier.phd$beta[[3]], beta.ab))

        phd.angle <- (1/4) * ((angles(phd.1$beta.hat[1:nvars,1], beta.a) + angles(phd.2$beta.hat[1:nvars,1],beta.b) )
                              + 2 * angles(phd.3$beta.hat, beta.ab))

        sir.angle <- (1/4) * (( angles(sir.1$beta.hat[1:nvars,1], beta.a) + angles(sir.2$beta.hat[1:nvars,1],beta.b) )
                              + 2 * angles(sir.3$beta.hat, beta.ab))

        semi.phd.angle <- (1/4) * ((angles(s.phd.1$beta[1:nvars,1], beta.a)+ angles(s.phd.2$beta[1:nvars,1],beta.b) )
                                   + 2 * angles(s.phd.3$beta, beta.ab))


        ######## projection difference norm

        hier.sir.proj <- (1/4) * ((proj.norm(hier.sir.betaAB[,1], beta.a) + proj.norm(hier.sir.betaAB[,2], beta.b) )
                                  + 2 * proj.norm(hier.sir.betaAB, beta.ab))

        hier.phd.proj <- (1/4) * ((proj.norm(hier.phd.betaAB[,1], beta.a) + proj.norm(hier.phd.betaAB[,2], beta.b) )
                                  + 2 * proj.norm(hier.phd.betaAB, beta.ab))

        semi.hier.phd.proj <- (1/4) * ((proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
                                       + 2 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

        phd.proj <- (1/4) * ((proj.norm(phd.1$beta.hat[1:nvars,1], beta.a) + proj.norm(phd.2$beta.hat[1:nvars,1],beta.b) )
                             + 2 * proj.norm(phd.3$beta.hat, beta.ab))

        sir.proj <- (1/4) * (( proj.norm(sir.1$beta.hat[1:nvars,1], beta.a) + proj.norm(sir.2$beta.hat[1:nvars,1],beta.b) )
                             + 2 * proj.norm(sir.3$beta.hat, beta.ab))

        semi.phd.proj <- (1/4) * ((proj.norm(s.phd.1$beta[1:nvars,1], beta.a)+ proj.norm(s.phd.2$beta[1:nvars,1],beta.b) )
                                  + 2 * proj.norm(s.phd.3$beta, beta.ab))

    }

    result.mat.cor <- result.mat.angle <- result.mat.proj <- data.frame(grid[sim.idx,], array(NA, dim = c(1, num.models + 1)))

    colnames(result.mat.cor) <-
        colnames(result.mat.angle) <-
        colnames(result.mat.proj) <-
        c(colnames(grid), "SNR",
          "SIR", "Hier SIR", "PHD", "Hier PHD",
          "Semi PHD", "Semi Hier PHD")


    result.mat.cor[1,(ncol(grid) + 1):ncol(result.mat.cor)] <- c(snr,
                                                                 sir.cor, hier.sir.cor,
                                                                 phd.cor, hier.phd.cor,
                                                                 s.phd.cor,
                                                                 semi.hier.phd.cor)


    result.mat.angle[1,(ncol(grid) + 1):ncol(result.mat.angle)] <- c(snr,
                                                                     sir.angle, hier.sir.angle,
                                                                     phd.angle, hier.phd.angle,
                                                                     semi.phd.angle,
                                                                     semi.hier.phd.angle)

    result.mat.proj[1,(ncol(grid) + 1):ncol(result.mat.proj)] <- c(snr,
                                                                   sir.proj, hier.sir.proj,
                                                                   phd.proj, hier.phd.proj,
                                                                   semi.phd.proj,
                                                                   semi.hier.phd.proj)


    print("semi phd angle:")
    print(semi.phd.angle)
    print("semi hier phd angle:")
    print(semi.hier.phd.angle)

    print(result.mat.angle)

    ## save results
    write.csv(result.mat.cor,   paste(seed.num, "_results_correlation.csv", sep = ""))
    write.csv(result.mat.angle, paste(seed.num, "_results_angle.csv",       sep = ""))
    write.csv(result.mat.proj,  paste(seed.num, "_results_proj_norm.csv",   sep = ""))
}



result.mat <- data.frame(grid[sim.idx,], array(NA, dim = c(1, 2)))

colnames(result.mat) <- c(colnames(grid), "SNR", "cv")

result.mat[1,(ncol(grid) + 1):ncol(result.mat)] <- c(snr, gcvs)

write.csv(result.mat, paste(seed.num, "_results_dim_selection.csv", sep = ""))
