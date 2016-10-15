

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


## these models have eta as zero
eta.zero.models <- c(2, 4)


d.mat <- data.matrix(expand.grid(1:max.d[1], 1:max.d[2], 0:max.d[3]))
dimnames(d.mat) <- NULL
d.mat <- d.mat[which(rowSums(d.mat) < 5),]

d.nums <- 1:nrow(d.mat)

d.nums <- 1:nrow(d.mat)

#    set up grid of parameters over which to simulate #
grid           <- expand.grid(nobs.vec, nvars, x.type, beta.type.vec, model.num.vec, d.nums)
colnames(grid) <- c("nobs", "nvars", "x.type", "beta.type", "model.num", "d.num")

n.jobs <- nrow(grid)

n.d <- nrow(d.mat)


n.sims <- 100
total.sims <- nrow(grid) * n.sims



n.unique.grids <- nrow(grid) / n.d

path <- "C:/Users/Jared/Dropbox/ACO/rehospitalization/sdr_hierarchical/simulation/results/"

res.mat <- read.csv(paste0(path, "results_dim_sel_4_10_15_2016.csv"), stringsAsFactors = FALSE)


df.2.use <- res.mat[1:n.d,]

for (i in 1:ncol(df.2.use))
{
    if (typeof(df.2.use[,i]) == "integer")
    {
        df.2.use[,i] <- NA_integer_
    } else if (is.character(df.2.use[,i]))
    {
        df.2.use[,i] <- NA_character_
    } else if (is.numeric(df.2.use[,i]))
    {
        df.2.use[,i] <- NA_real_
    } else
    {
        df.2.use[,i] <- NA
    }
}



res.list <- rep(list(rep(list(df.2.use), n.sims)), n.unique.grids)

for (i in 1:nrow(res.mat))
{

    res.cur <- res.mat[i,,drop=FALSE]

    seed.num <- res.cur$sim_num

    ## get index for the correct simulation parameters
    ## for this specific job (job id = seed.num) # job ids start at 0
    sim.idx <- (seed.num + 1 - 1) %% nrow(grid) + 1

    rep.num <- (seed.num + 1 - sim.idx) / nrow(grid) + 1
    set.seed(seed.num)

    grid.idx  <- c(sim.idx - 1) %% n.unique.grids + 1

    nvars     <- grid$nvars[sim.idx]
    x.type    <- grid$x.type[sim.idx]
    nobs      <- grid$nobs[sim.idx]
    beta.type <- grid$beta.type[sim.idx]
    model.num <- grid$model.num[sim.idx]
    d.cur     <- as.vector(d.mat[grid$d.num[sim.idx], ])

    res.list[[grid.idx]][[rep.num]][res.cur$d.num,] <- res.cur
}

res.frame <- array(NA, dim = c(n.sims, 8L))
colnames(res.frame) <- c("rank", "num.correct", "top1", "top3", "top5", "top10", "big.enough.dim", "dims.completed")
rank.res <- rep(list(res.frame), n.unique.grids)

for (g in 1:n.unique.grids)
{
    for (s in 1:n.sims)
    {
        res.cur <- na.omit(res.list[[g]][[s]])
        if (nrow(res.cur) > 1)
        {
            d.mat.cur <- d.mat[res.cur$d.num,,drop=FALSE]
            min.idx <- which.min(res.cur$vic)
            best.d  <- d.mat.cur[min.idx,]

            d.vec <- apply(d.mat.cur, 1, function(r) paste(r, collapse = ","))


            ranks <- rank(res.cur$vic)

            model.num <- res.cur$model.num[1]

            if ( !(model.num %in% eta.zero.models) )
            {
                d.true <- c(1, 1, 1)

            } else
            {
                d.true <- c(1, 1, 0)
            }

            d.true.v  <- paste(d.true, collapse = ",")

            rank.true <- ranks[which(d.vec == d.true.v)]
            #rank.true <- sample(ranks, 1)

            if (length(rank.true))
            {


                top.10 <- 1 * (rank.true <= 10)
                top.5  <- 1 * (rank.true <= 5)
                top.3  <- 1 * (rank.true <= 3)
                top.1  <- 1 * (rank.true == 1)


                #best.d <- d.mat.cur[sample.int(nrow(d.mat.cur), 1),]
                num.correct <- sum(best.d == d.true)
                rank.res[[g]][s,] <- c(rank.true, num.correct, top.1, top.3, top.5, top.10,
                                       1 * all(best.d >= d.true),
                                       nrow(res.cur))
            }
        }
    }
}

round(res.avg <- sapply(rank.res, function(mat) colMeans(mat[mat[,match("dims.completed", colnames(mat))] > 7,], na.rm = TRUE)), 3)
grid[1:6,]

res.2.table <- cbind(grid[1:6,1:2], Model=c(1, 1, 1, 2, 2, 2), round(t(res.avg)[,1:5], 3) )

library(stargazer)

stargazer(res.2.table, type = "text", summary = FALSE)
stargazer(res.2.table, type = "latex", summary = FALSE)

