


opt.est.eqn <- function(init, est.eqn, est.eqn.grad,
                        opt.method = c("lbfgs2",
                                       "ucminf",
                                       "lbfgs.x",
                                       "bfgs.x",
                                       "bfgs",
                                       "lbfgs",
                                       "spg"),
                        nn = 0.75,
                        optimize.nn = FALSE, maxit = 100, verbose = FALSE)
{
    opt.method <- match.arg(opt.method)

    if (optimize.nn)
    {
        init.par <- c(log(nn / (1 - nn)), init)
    } else
    {
        init.par <- init
    }

    if (opt.method == "bfgs")
    {
        slver <-   optim(par     = init.par, # beta.init[(d+1):nrow(beta.init),],
                         fn      = est.eqn,
                         gr      = est.eqn.grad,
                         method  = "BFGS",
                         nn.val  = nn,
                         optimize.nn = optimize.nn,
                         control = list(maxit = maxit, abstol = 1e-10))
        if (optimize.nn) slver$par <- slver$par[-1]
    } else if (opt.method == "lbfgs")
    {
        slver <-   optim(par     = init.par,
                         fn      = est.eqn,
                         gr      = est.eqn.grad,
                         method  = "L-BFGS",
                         nn.val  = nn,
                         optimize.nn = optimize.nn,
                         control = list(maxit = maxit, factr = 1e-10))
        if (optimize.nn) slver$par <- slver$par[-1]
    } else if (opt.method == "lbfgs2")
    {
        slver <-   lbfgs(vars      = init.par,
                         call_eval = est.eqn,
                         call_gra  = est.eqn.grad,
                         nn.val    = nn,
                         optimize.nn = optimize.nn,
                         invisible = 1 * (!verbose),
                         max_iterations = maxit)
        if (optimize.nn) slver$par <- slver$par[-1]
    } else if (opt.method == "bfgs.x")
    {
        slver <-   optimx(par     = init.par,
                          fn      = est.eqn,
                          #gr      = est.eqn.grad,
                          control = list(trace = 1 * verbose,
                                         maxit = maxit,
                                         kkt = FALSE),
                          method  = "BFGS",
                          nn.val  = nn,
                          optimize.nn = optimize.nn)
        if (is.null(slver$par) & names(slver)[1] == "p1")
        {
            if (optimize.nn)
            {
                if (verbose) print(paste("optimized nn:", exp(slver$p1) / (1 + exp(slver$p1))))
                slver <- list(par = unlist(slver[2:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            } else
            {
                slver <- list(par = unlist(slver[1:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            }
        }
    } else if (opt.method == "lbfgs.x")
    {

        slver <-   optimx(par     = init.par, # c(log(nn / (1 - nn)), init.par),
                          fn      = est.eqn,
                          #gr      = est.eqn.grad,
                          control = list(trace = 1 * verbose,
                                         maxit = maxit,
                                         kkt = FALSE),
                          method  = "L-BFGS-B",
                          nn.val  = nn,
                          optimize.nn = optimize.nn)
        if (is.null(slver$par) & names(slver)[1] == "p1")
        {
            if (optimize.nn)
            {
                if (verbose) print(paste("optimized nn:", exp(slver$p1) / (1 + exp(slver$p1))))
                slver <- list(par = unlist(slver[2:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            } else
            {
                slver <- list(par = unlist(slver[1:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            }
        }
    } else
    {
        slver <-   optimx(par     = init.par,
                          fn      = est.eqn,
                          #gr      = est.eqn.grad,
                          control = list(trace = 1 * verbose,
                                         maxit = maxit,
                                         kkt = FALSE),
                          method  = opt.method,
                          nn.val  = nn,
                          optimize.nn = optimize.nn)
        if (is.null(slver$par) & names(slver)[1] == "p1")
        {
            if (optimize.nn)
            {
                if (verbose) print(paste("optimized nn:", exp(slver$p1) / (1 + exp(slver$p1))))
                slver <- list(par = unlist(slver[2:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            } else
            {
                slver <- list(par = unlist(slver[1:length(init.par)]),
                              value = slver$value,
                              convcode = slver$convcode)
            }
        }
    }
    slver
}

try.nn <- function(nn.vals = c(0.1, 0.25, 0.5, 0.75, 0.9, 0.95),
                   init, est.eqn, est.eqn.grad,
                   opt.method = c("lbfgs.x", "bfgs", "lbfgs2",
                                  "bfgs.x",
                                  "lbfgs",
                                  "spg", "ucminf"),
                   optimize.nn = FALSE, maxit = 10, verbose = FALSE)
{
    n.try  <- length(nn.vals)
    values <- numeric(n.try)
    for (i in 1:n.try)
    {
        slver.cur <- opt.est.eqn(init         = init,
                                 est.eqn      = est.eqn,
                                 est.eqn.grad = est.eqn.grad,
                                 opt.method   = opt.method,
                                 nn           = nn.vals[i],
                                 optimize.nn  = optimize.nn,
                                 maxit        = maxit,
                                 verbose      = verbose > 1)

        if (verbose) cat("nn:", nn.vals[i], "; f(x) =", slver.cur$value, "\n")

        values[i] <- slver.cur$value
    }
    nn.vals[which.min(values)]
}
