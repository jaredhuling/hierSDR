
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
d.mat <- d.mat[which(rowSums(d.mat) < 5),]

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

if (model.num == 5)
{
    sd.sim <- 0.5
}

print(grid[sim.idx,])

Proj <- function(b) b %*% solve(crossprod(b), t(b))

proj.norm <- function(b1, b2)
{
    norm(Proj(b1) - Proj(b2), type = "F")
}

## generate data
#  source("C:/Users/Jared/Documents/GitHub/hierSDR/test/condor_generate_data_linear.R")
source("condor_generate_data_linear_3.R")

## fit SDR models

semi.hier.phd  <- semi.phd.hier.separate(x.list, drop(y),
                                         d = d.cur, k = 5,
                                         h = exp(seq(log(0.5), log(25), length.out = 25)),
                                         maxit = 50, maxk = 2500)
#

vic  <- semi.hier.phd$vic
mse  <- semi.hier.phd$mse
sse  <- semi.hier.phd$sse
dd1   <- d.cur[1]
dd2   <- d.cur[2]
dd3   <- d.cur[3]
print(d.cur)

if ( !(model.num %in% eta.zero.models) )
{
    d.true <- c(1, 1, 1)

} else
{
    d.true <- c(1, 1, 0)
}




if ( !(model.num %in% eta.zero.models) )
{


    ######## projection difference norm

    #hier.sir.proj <- (1/5) * ((proj.norm(hier.sir.betaAB[,1], beta.a) + proj.norm(hier.sir.betaAB[,2], beta.b) )
    #                           + 3 * proj.norm(hier.sir.betaAB, beta.ab))

    #hier.phd.proj <- (1/5) * ((proj.norm(hier.phd.betaAB[,1], beta.a) + proj.norm(hier.phd.betaAB[,2], beta.b) )
    #                           + 3 * proj.norm(hier.phd.betaAB, beta.ab))

    semi.hier.phd.proj <- (1/5) * ((proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
                                    + 3 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

    #phd.proj <- (1/5) * ((proj.norm(phd.1$beta.hat[1:nvars,1], beta.a) + proj.norm(phd.2$beta.hat[1:nvars,1],beta.b) )
    #                      + 3 * proj.norm(phd.3$beta.hat, beta.ab))

    #sir.proj <- (1/5) * (( proj.norm(sir.1$beta.hat[1:nvars,1], beta.a) + proj.norm(sir.2$beta.hat[1:nvars,1],beta.b) )
    #                      + 3 * proj.norm(sir.3$beta.hat, beta.ab))

    #semi.phd.proj <- (1/5) * ((proj.norm(s.phd.1$beta[1:nvars,1], beta.a)+ proj.norm(s.phd.2$beta[1:nvars,1],beta.b) )
    #                           + 3 * proj.norm(s.phd.3$beta, beta.ab))

    (1/5) * (( proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
             + 3 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

} else
{


    ######## projection difference norm

    #hier.sir.proj <- (1/4) * ((proj.norm(hier.sir.betaAB[,1], beta.a) + proj.norm(hier.sir.betaAB[,2], beta.b) )
    #                          + 2 * proj.norm(hier.sir.betaAB, beta.ab))

    #hier.phd.proj <- (1/4) * ((proj.norm(hier.phd.betaAB[,1], beta.a) + proj.norm(hier.phd.betaAB[,2], beta.b) )
    #                          + 2 * proj.norm(hier.phd.betaAB, beta.ab))

    semi.hier.phd.proj <- (1/4) * ((proj.norm(semi.hier.phd$beta[[1]], beta.a) + proj.norm(semi.hier.phd$beta[[2]],beta.b) )
                                   + 2 * proj.norm(semi.hier.phd$beta[[3]], beta.ab))

    #phd.proj <- (1/4) * ((proj.norm(phd.1$beta.hat[1:nvars,1], beta.a) + proj.norm(phd.2$beta.hat[1:nvars,1],beta.b) )
    #                     + 2 * proj.norm(phd.3$beta.hat, beta.ab))

    #sir.proj <- (1/4) * (( proj.norm(sir.1$beta.hat[1:nvars,1], beta.a) + proj.norm(sir.2$beta.hat[1:nvars,1],beta.b) )
    #                     + 2 * proj.norm(sir.3$beta.hat, beta.ab))

    #semi.phd.proj <- (1/4) * ((proj.norm(s.phd.1$beta[1:nvars,1], beta.a)+ proj.norm(s.phd.2$beta[1:nvars,1],beta.b) )
    #                          + 2 * proj.norm(s.phd.3$beta, beta.ab))

}



result.mat <- data.frame(grid[sim.idx,], array(NA, dim = c(1, 11)))

colnames(result.mat) <- c(colnames(grid), "SNR", "vic", "proj.norm", "mse1", "mse2", "mse3", "sse1", "sse2", "sse3", "d1", "d2", "d3")

print(c(snr, vic, semi.hier.phd.proj, mse, sse, dd1, dd2, dd3))

result.mat[1,(ncol(grid) + 1):ncol(result.mat)] <- c(snr, vic, semi.hier.phd.proj, mse, sse, dd1, dd2, dd3)

write.csv(result.mat, paste(seed.num, "_results_dim_selection.csv", sep = ""))
