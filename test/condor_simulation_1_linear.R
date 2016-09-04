
arg <- commandArgs()      # Gets the process number
repnum <- eval(parse(text=arg[[3]]))  # Gets a workable character string from the argument
seed.num <- repnum
print(seed.num)


## load required packages
library(hierSDR)



## set up simulation grid
nobs.vec <- c(250, 500, 1000, 2000)
nvars <- 50
ncats <- 2
x.type <- c("simple", "complicated")
nobs.test <- 10000
simtype <- "regr"
nsims <- 500

#    set up grid of parameters over which to simulate #
grid <- expand.grid(nobs.vec, nvars, x.type)
colnames(grid) <- c("nobs", "nvars", "x.type")

n.jobs <- nrow(grid)

## get index for the correct simulation parameters
## for this specific job (job id = seed.num) # job ids start at 0
sim.idx <- (seed.num + 1 - 1) %% nrow(grid) + 1

set.seed(seed.num)

nvars  <- grid$nvars[sim.idx]
x.type <- grid$x.type[sim.idx]
nobs   <- grid$nobs[sim.idx]




