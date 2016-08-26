
set.seed(123)
x <- matrix(rnorm(2000 * 5), ncol = 5)
y <- x[,1] - x[,2] + rnorm(nrow(x))

grad.fun <- function(beta, x, y) {
    gradient <- (1 / nrow(x)) * (crossprod(x, (x %*% beta) - y )  )
    drop(gradient)
}

grad.fun.constr <- function(beta, x, y, constr) {
    Proj <- constr %*% solve(crossprod(constr), t(constr))
    gradient <- drop((1 / nrow(x)) * (crossprod(x, (x %*% beta) - y )  ))
    adj.fact <- drop(Proj %*% gradient)
    gradient - adj.fact
}




# define gradient descent update algorithm
grad.descent <- function(x, y, alpha = 0.05, maxit = 50L, tol = 1e-8){
    beta <- numeric(ncol(x)) # Initialize the parameters

    for (i in 1:maxit) {
        prev <- beta
        beta <- beta - alpha  * grad.fun(beta, x, y)
        if (all(abs( (beta - prev) / abs(prev) ) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    beta
}

backtracking.linesearch <- function(beta, direction, fn, gr, c = 1e-4, rho = 0.9, maxit = 50L, ...)
{
    alpha = 1
    itr <- 0

    f_x         <- fn(beta , ...)
    f_x_scratch <- fn(beta + direction, ...)
    gr_scratch  <- gr(beta + direction, ...)
    gxp <- drop(crossprod(drop(gr_scratch), drop(direction)))

    while(f_x_scratch > f_x + c * alpha * gxp)
    {
        itr <- itr + 1
        if (itr > maxit)
        {
            break
        }

        alpha <- alpha * rho

        f_x_scratch <- fn(beta + alpha * direction, ...)
    }
    alpha
}

backtracking.linesearch.constr <- function(beta, direction, fn, gr, constraints, c = 1e-4, rho = 0.9, maxit = 50L, ...)
{
    alpha = 1
    itr <- 0

    f_x         <- fn(beta , ...)
    f_x_scratch <- fn(beta + direction, ...)
    gr_scratch  <- gr(beta + direction, constr = constraints, ...)
    gxp <- drop(crossprod(drop(gr_scratch), drop(direction)))

    while(f_x_scratch > f_x + c * alpha * gxp)
    {
        itr <- itr + 1
        if (itr > maxit)
        {
            break
        }

        alpha <- alpha * rho

        f_x_scratch <- fn(beta + alpha * direction, ...)
    }
    alpha
}


bfgs <- function(x, y, fn = NULL, gr = grad.fun, maxit = 50L, tol = 1e-8, alpha = 1, c = 1e-4, rho = 0.9)
{
    p <- ncol(x)
    B.inv.cur <- diag(p)
    beta <- numeric(p)
    grad.cur <- gr(beta, x, y)
    direc <- -grad.cur
    if (!is.null(fn))
    {
        alpha <- backtracking.linesearch(beta, direction = direc, fn = fn, gr = gr, x = x, y = y, c = c, rho = rho)
        #cat("alpha", alpha, "\n")
    }
    beta <- beta + alpha * direc

    for (i in 1:maxit)
    {
        s <- alpha * direc
        gradient <- grad.cur
        B.inv <- B.inv.cur
        grad.cur <- gr(beta, x, y)
        gy <- grad.cur - gradient
        norm.s.gy <- drop(crossprod(s, gy))
        #B.inv.cur <- (diag(p) - tcrossprod(s, gy) / norm.s.gy) %*% B.inv %*% (diag(p) - tcrossprod(gy, s) / norm.s.gy)
        #B.inv.cur <- B.inv.cur + tcrossprod(s) / norm.s.gy
        B.inv.y <- drop(B.inv %*% gy)
        B.inv.cur <- B.inv + ((norm.s.gy + drop(crossprod(gy, B.inv.y)) ) / (norm.s.gy ^ 2) ) * tcrossprod(s)
        B.inv.cur <- B.inv.cur - (tcrossprod(B.inv.y, s) + tcrossprod(s, B.inv.y)) / norm.s.gy

        direc <- -B.inv.cur %*% grad.cur
        if (!is.null(fn))
        {
            alpha <- backtracking.linesearch(beta, direction = direc, fn = fn, gr = gr, x = x, y = y, c = c, rho = rho)
            #cat("alpha", alpha, "\n")
        }
        beta <- beta + alpha * direc
        if (all(abs(grad.cur) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    drop(beta)
}


lbfgs <- function(x, y, fn = NULL, gr = grad.fun, maxit = 50L, tol = 1e-8, alpha = 1, c = 1e-4, rho = 0.9)
{
    p <- ncol(x)
    B.inv.cur <- diag(p)
    beta <- numeric(p)
    grad.cur <- gr(beta, x, y)
    direc <- -grad.cur
    if (!is.null(fn))
    {
        alpha <- backtracking.linesearch(beta, direction = direc, fn = fn, gr = gr, x = x, y = y, c = c, rho = rho)
        #cat("alpha", alpha, "\n")
    }
    beta <- beta + alpha * direc

    for (i in 1:maxit)
    {
        s <- alpha * direc
        gradient <- grad.cur
        B.inv <- B.inv.cur
        grad.cur <- gr(beta, x, y)
        gy <- grad.cur - gradient
        norm.s.gy <- drop(crossprod(s, gy))
        #B.inv.cur <- (diag(p) - tcrossprod(s, gy) / norm.s.gy) %*% B.inv %*% (diag(p) - tcrossprod(gy, s) / norm.s.gy)
        #B.inv.cur <- B.inv.cur + tcrossprod(s) / norm.s.gy
        B.inv.y <- drop(B.inv %*% gy)
        B.inv.cur <- B.inv + ((norm.s.gy + drop(crossprod(gy, B.inv.y)) ) / (norm.s.gy ^ 2) ) * tcrossprod(s)
        B.inv.cur <- B.inv.cur - (tcrossprod(B.inv.y, s) + tcrossprod(s, B.inv.y)) / norm.s.gy

        direc <- -B.inv.cur %*% grad.cur
        if (!is.null(fn))
        {
            alpha <- backtracking.linesearch(beta, direction = direc, fn = fn, gr = gr, x = x, y = y, c = c, rho = rho)
            #cat("alpha", alpha, "\n")
        }
        beta <- beta + alpha * direc
        if (all(abs(grad.cur) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    drop(beta)
}




# define gradient descent update algorithm
grad.descent.constr <- function(x, y, constraints = matrix(0, nrow = ncol(x)), alpha = 0.05, maxit = 50L, tol = 1e-8){
    beta <- numeric(ncol(x)) # Initialize the parameters

    for (i in 1:maxit) {
        prev <- beta
        gradient <- grad.fun(beta, x, y)
        Proj <- constraints %*% solve(crossprod(constraints), t(constraints))
        adj.fact <- Proj %*% gradient
        beta <- beta - alpha  *  (gradient - adj.fact)
        if (all(abs( (beta - prev) / abs(prev) ) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    drop(beta)
}

bfgs.constr <- function(x, y, fn = NULL, gr = grad.fun.constr,
                        constraints,
                        maxit = 50L, tol = 1e-8,
                        alpha = 1, c = 1e-4, rho = 0.9)
{
    p <- ncol(x)
    B.inv.cur <- diag(p)
    beta <- numeric(p)
    grad.cur <- gr(beta, x, y, constraints)
    direc <- -grad.cur
    if (!is.null(fn))
    {
        alpha <- backtracking.linesearch.constr(beta, direction = direc,
                                                fn = fn, gr = gr,
                                                constraints = constraints,
                                                c = c, rho = rho,
                                                x = x, y = y)
        #cat("alpha", alpha, "\n")
    }
    beta <- beta + alpha * direc

    for (i in 1:maxit)
    {
        s <- direc
        gradient <- grad.cur
        B.inv <- B.inv.cur
        grad.cur <- gr(beta, x, y, constraints)
        gy <- grad.cur - gradient
        norm.s.gy <- drop(crossprod(s, gy))
        B.inv.cur <- (diag(p) - tcrossprod(s, gy) / norm.s.gy) %*% B.inv %*% (diag(p) - tcrossprod(gy, s) / norm.s.gy)
        B.inv.cur <- B.inv.cur + tcrossprod(s) / norm.s.gy

        direc <- -B.inv.cur %*% grad.cur
        if (!is.null(fn))
        {
            alpha <- backtracking.linesearch.constr(beta, direction = direc,
                                                    constraints = constraints,
                                                    fn = fn, gr = gr,
                                                    c = c, rho = rho,
                                                    x = x, y = y)
            #cat("alpha", alpha, "\n")
        }
        beta <- beta + alpha * direc
        if (all(abs(grad.cur) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    drop(beta)
}


func.2.min <- function(beta)
{
    sum((x %*% beta - y)^2 )
}


constr.mat <- cbind(c(0,0,1,-1,0),
                    c(0,0,1,0,-1))

eq.fun <- function(beta, x, y)
{
    constr.mat <- cbind(c(0,0,1,-1,0),
                        c(0,0,1,0,-1))
    print(str(beta))
    if (is.null(dim(beta))) beta <- matrix(beta, ncol = 5)
    print(dim(t(constr.mat)))
    drop(t(constr.mat) %*% t(beta))
}

(beta.ls <- drop(solve(crossprod(x), crossprod(x, y))))

bfgs(x, y)
grad.descent(x, y, maxit = 500, alpha = 0.21)
grad.descent.constr(x, y, maxit = 500, alpha = 0.21, constraints = constr.mat)
bfgs.constr(x, y, constraints = constr.mat)

library(Rsolnp)




fn1=function(x)
{
    exp(x[1]*x[2]*x[3]*x[4]*x[5])
}
eqn1=function(x){
    z1=x[1]*x[1]+x[2]*x[2]+x[3]*x[3]+x[4]*x[4]+x[5]*x[5]
    z2=x[2]*x[3]-5*x[4]*x[5]
    z3=x[1]*x[1]*x[1]+x[2]*x[2]*x[2]
    return(c(z1,z2,z3))
}
x0 = c(-2, 2, 2, -1, -1)
powell=solnp(x0, fun = fn1, eqfun = eqn1, eqB = c(10, 0, -1))


eq.fun <- function(beta)
{
    z1 <- beta[3] - beta[4]
    z2 <- beta[3] - beta[5]
    c(z1, z2)
}

x0 <- rnorm(5)
lsminc <- solnp(x0, fun = func.2.min, eqfun = eq.fun, eqB = c(0, 0))


LHS <- rbind(cbind(crossprod(x), constr.mat),
             cbind(t(constr.mat), matrix(0, ncol = 2, nrow = 2)) )
RHS <- c(crossprod(x, y), rep(0, 2))

drop(solve(LHS, RHS))[1:ncol(x)]
lsminc$pars
grad.descent.constr(x, y, maxit = 500, alpha = 0.21, constraints = constr.mat)
bfgs.constr(x, y, constraints = constr.mat)



logistic.loss <- function(beta, x, y, ...)
{
    y <- 2 * y - 1
    sum(log(1 + exp(-y * (x %*% beta))))
}

logistic.grad <- function(beta, x, y)
{
    y <- 2 * y - 1
    -colSums( ((  1/( 1 + exp(y * drop(x %*% beta)) )  ) * y) * x)
}

logistic.grad.nd <- function(beta, xx, yy)
{
    y <- 2 * yy - 1
    -colSums( ((  1/( 1 + exp(y * drop(xx %*% beta)) )  ) * y) * xx)
}


logistic.grad.descent <- function(x, y, alpha = 0.05, maxit = 50L, tol = 1e-8)
{
    beta <- numeric(ncol(x)) # Initialize the parameters

    for (i in 1:maxit) {
        prev <- beta
        beta <- beta - alpha  * logistic.grad(beta, x, y)
        if (all(abs( (beta - prev) / abs(prev) ) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    beta
}

logistic.grad.descent <- function(x, y, alpha = 0.05, maxit = 50L, tol = 1e-8)
{
    beta <- numeric(ncol(x)) # Initialize the parameters

    for (i in 1:maxit) {
        prev <- beta
        beta <- beta - alpha  * logistic.grad(beta, x, y)
        if (all(abs( (beta - prev) / abs(prev) ) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    beta
}

logistic.grad.descent.constr <- function(x, y, alpha = 0.05, maxit = 50L, tol = 1e-8, constraints)
{
    beta <- numeric(ncol(x)) # Initialize the parameters

    for (i in 1:maxit) {
        prev <- beta
        beta <- beta - alpha  * logistic.grad.constr(beta, x, y, constraints)
        if (all(abs( (beta - prev) / abs(prev) ) <= tol))
        {
            cat("Converged in ", i, "iterations \n")
            break
        }
    }
    beta
}

set.seed(123)
prob <- 1 / (1 + exp(-x[,1] * 0.5 + x[,2] * 0.5))
yl <- rbinom(nrow(x), 1, prob)

logistic.loss.np <- function(beta)
{
    y2 <- 2 * yl - 1
    sum(log(1 + exp(-y2 * (x %*% beta))))
}



logistic.grad.constr <- function(beta, x, y, constr)
{
    Proj <- constr %*% solve(crossprod(constr), t(constr))
    gradient <- logistic.grad(beta, x, y)
    adj.fact <- drop(Proj %*% gradient)
    gradient - adj.fact
}

x0 <- rnorm(5)
lsminc.logistic <- solnp(x0, fun = logistic.loss.np, eqfun = eq.fun, eqB = c(0, 0))

## unconstrained
(coef(glm(yl ~ x - 1, family = binomial())))
logistic.grad.descent(x, yl, maxit = 5000, alpha = 0.001)
optim(rep(0, ncol(x)), fn = logistic.loss, gr = logistic.grad, x = x, y = yl, method = "BFGS")$par
system.time(replicate(100, optim(rep(0, ncol(x)), fn = logistic.loss, gr = logistic.grad, x = x, y = yl, method = "BFGS")$par))
system.time(replicate(100, optim(rep(0, ncol(x)), fn = logistic.loss, gr = logistic.grad, x = x, y = yl, method = "L-BFGS", control = list(factr = 1e-12))$par))
system.time((dm <- bfgs(x = x, y = yl, gr = logistic.grad, maxit = 500, alpha = 0.5)))
system.time((dm <- bfgs(x = x, y = yl, fn = logistic.loss, gr = logistic.grad, maxit = 500, alpha = 0.5)))
system.time((dm <- bfgs(x = x, y = yl, fn = logistic.loss, gr = logistic.grad, maxit = 500, alpha = 0.5, c = 0.01, rho = 0.2)))
system.time(replicate(100, bfgs(x = x, y = yl, fn = logistic.loss, gr = logistic.grad, maxit = 500, alpha = 0.5, c = 0.01, rho = 0.2)))

(dm <- bfgs(x = x, y = yl, fn = logistic.loss, gr = logistic.grad, maxit = 500, alpha = 0.5, c = 0.01, rho = 0.2))

## constrained
logistic.grad.descent.constr(x, yl, maxit = 5000, alpha = 0.001, constraints = constr.mat)
lsminc.logistic <- solnp(x0, fun = logistic.loss.np, eqfun = eq.fun, eqB = c(0, 0));lsminc.logistic$par

(dmc <- bfgs.constr(x = x, y = yl, gr = logistic.grad.constr, constraints = constr.mat, maxit = 500, alpha = 0.01))
(dmc <- bfgs.constr(x = x, y = yl, fn = logistic.loss, gr = logistic.grad.constr, constraints = constr.mat, maxit = 500, alpha = 0.01, c = 0.01, rho = 0.2))
optim(rep(0, ncol(x)), fn = logistic.loss, gr = logistic.grad.constr, constr = constr.mat, x = x, y = yl, method = "L-BFGS", control = list(factr = 1e-12))$par
optim(runif(ncol(x)), fn = logistic.loss, gr = logistic.grad.constr, constr = constr.mat, x = x, y = yl, method = "L-BFGS", control = list(factr = 1e-12))$par


logistic.grad(rep(0, ncol(x)), x = x, y = y)
library(numDeriv)
numDeriv::grad(logistic.grad.nd, rep(0, ncol(x)), xx = x, yy = y)
