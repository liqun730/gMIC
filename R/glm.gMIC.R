# transform a covariance matrix into a block diagonal matrix, each block is the
# covariance matrix (or its inverse if inv=TRUE) of a group of variables
cov.bdiag <- function(M, group, inv=TRUE){
  uniq <- unique(group)
  L <- list()
  i <- 1
  for(g in uniq){
    ix <- group == g
    L[[i]] <- if(inv) solve(M[ix,ix]) else M[ix,ix]
    i <- i+1
  }
  bdiag(L)
}

# get pvalue and coverage probability for grouped variables:
group.pval <- function(formula, data, family, intercept, gamma_hat, group){
  if(intercept) group <- c(min(group)-1,group)
  uniq.grp <- unique(group); n <- nrow(data)
  grp.size <- sapply(uniq.grp, function(x) sum(group==x))
  fit0 <- glm(formula, data = data, family = family,
              start=gamma_hat, control=glm.control(epsilon=1e20, maxit=1, trace=F))
  VCOV.gamma <- summary(fit0)$"cov.scaled"
  VCOV.inv <- cov.bdiag(VCOV.gamma,group)
  lst_hat <- split(gamma_hat,group); lst_hat <- bdiag(lst_hat)
  stat_pval <- diag(t(lst_hat)%*%VCOV.inv%*%(lst_hat))

  pvalues <- pchisq(stat_pval, df=grp.size, lower.tail=F)
  return(pvalues)
}

#' The gMIC Function for (Group) Variable Selection in Generalized Linear Model
#'
#' @param formula An object of class \code{\link[stats]{formula}}, with the response on the left of a \code{~} operator, and the terms on the right.
#' @param family A description of the error distribution and link function to be used in the model. Preferably for computational speed,
#' this is a character string naming a family function among the following three choices: \code{"gaussian"},
#' \code{"binomial"}, or \code{"poisson"}.  Otherwise, it has to be a family function or the result of a call to a family function that
#' can be called for by \code{\link[stats]{glm.fit}}. See \code{\link[stats]{family}} for details of family functions.
#' @param data A data.frame in which to interpret the variables named in the \code{formula} argument.
#' @param group A vector indicating the group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
#' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
#' @param beta0 A vector eqaul to the initial value for the model parameter, default is NULL.
#' @param criterion A string indicating the type of information criterion ("AIC" or "BIC") to approximate. Default is "BIC".
#' @param lambda0 A number, the user-specified penalty parameter for model complexity. If \code{criterion="AIC"} or \code{"BIC"}, the value
#' of \code{lambda0} will be ignored.
#' @param a0 The approximation parameter of the gMIC method.
#' @param scale.x A boolean indicating whether or not to studentize the features. Default is \code{TRUE}.
#' @param orthogonal.x A boolean indicating whether or not to orthogonalize the features within each group. Default is true. See \code{link{orthogonalize}} for details.
#' @param rounding.digits Number of digits after the decimal point for rounding-up estiamtes. Default value is 4.
#' @param optim.method Optimization method for gMIC, one of c("GenSA", "BFGS", "ADAM"), indicating we use GenSA, BFGS, or ADAM for gMIC optimmization.
#' Default is BFGS. For unknown methods specified by user, the default with be used.
#' @param lower The lower bounds for the search space in \code{GenSA}. The default is -10 (\eqn{p} by \eqn{1} vector).
#' @param upper The upper bounds for the search space in \code{GenSA}. The default is +10 (\eqn{p} by \eqn{1} vector).
#' @param maxit.global  Maximum number of iterations allowed for the global optimization algorithm \code{SANN}. Default value is 100.
#' @param maxit.local Maximum number of iterations allowed for the local optimizaiton algorithm \code{BFGS}. Default value is 100.
#' @param epsilon The convergence tolerance.
#' @param stepsize The stepsize (or learning rate) for optim.method = "GD" and "ADAM".
#' @param details Logical value: if \code{TRUE}, detailed results will be printed out when running \code{glm.gMIC}.
#' @return A list of objects as follows,
#' \describe{
#' \item{coefficients}{The estimates for gamma, the standard error of gamma, pvalues of gamma, and estimates of beta.}
#' \item{group.pvalues}{The group-level p-values for each group of variables.}
#' }
#' @import GenSA
#' @import Matrix
#' @examples
#' library(MASS)
#' library(Matrix)
#' n=500;a=100
#' sig <- function(k, rho){
#'   m = matrix(rho,nrow=k,ncol=k)
#'   diag(m) <- 1
#'   return(m)
#' }
#' bt = c(1,1,1, .5,.5,.5, rep(0,24)); p = length(bt)
#' group = c(1,1,1, 2,2,2,rep(3:10,each=3))
#' rho1 = 0.1; rho2 = 0.6
#' COV = sig(p,rho1) + bdiag(rep(list(sig(3,rho2)),p/3))
#' set.seed(1234)
#' X = mvrnorm(n,rep(0,p),COV)
#' z = X%*%bt; pr = 1/(1+exp(-z))
#' y = rbinom(n,1,pr)
#' Xy = as.data.frame(cbind(X,y))
#' dim(Xy)
#' names(Xy) <- c(paste("X",1:30,sep=""),"y")
#' fit <- glm.gMIC(y~.-1,group=group,family="binomial",a0=a,data=Xy,orthogonal.x=T)
#'
#' @export
#'
glm.gMIC <- function(formula, family = c("gaussian", "binomial", "poisson"), data, group=NULL,
                     beta0=NULL, criterion ="BIC", lambda0=0, a0=NULL,
                     scale.x=FALSE, orthogonal.x=FALSE, # Do we orthogonalize X?
                     rounding.digits = 4,				  # Rounding digits for final result
                     optim.method = "BFGS", lower=NULL, upper=NULL,
                     maxit.global=100, maxit.local=100, epsilon=1e-6, 
                     stepsize = 0.01, details=FALSE)
{
  call <- match.call()
  # CHECK THE family= ARGUMENT
  family0 <- family
  if (is.character(family)) family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family)) family <- family()
  if (is.null(family$family)) {print(family); stop("'family0' not recognized")}

  # CHECK THE data= ARGUMENT
  if (missing(data)) data <- environment(formula)
  # OBTAIN THE DESIGN MATRIX X AND RESPONSE y
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1L) {
    nm <- rownames(Y); dim(Y) <- NULL
    if (!is.null(nm)) names(Y) <- nm
  }
  yname <- as.character(formula[[2]]); # EXTRACT RESPONSE NAME FROM THE FORMULA
  X <- if (!is.empty.model(mt)) model.matrix(mt, mf)
  else matrix(, NROW(Y), 0L)
  
  Xnames <- colnames(X)
  intercept <- Xnames[1]=="(Intercept)"
  if(is.null(group)) {
    if(intercept) group <- (1:(ncol(X)-1))
    else group <- (1:ncol(X))
  }

  # Standardize X
  if (scale.x || orthogonal.x) {
    if (intercept){
      X1 <- X[, -1]
      if(scale.x) { X1 <- scale(X1, center = TRUE, scale = TRUE); sc = attr(X1,"scaled:scale"); adj <- attr(X1,"scaled:center") }
      if(orthogonal.x) X1 <- orthogonalize(X1, group)
      X[, -1] <- X1
    } else {
      if(scale.x) { X <- scale(X, center = TRUE, scale = TRUE); sc = attr(X,"scaled:scale")}
      if(orthogonal.x) X <- orthogonalize(X, group)
    }
  }

  # DETERMINE lambda
  n <- NROW(Y);
  if (is.null(a0)) a0 <- min(n, 100)
  lambda <- if (criterion =="BIC") log(n) else if (criterion =="AIC") 2 else lambda0
  
  # DETERMINE OPTIMIZATION METHODS
  if(!optim.method %in% c("BFGS", "GenSA", "ADAM")) {
    warning(paste0("Unknow optimization method: \"", optim.method, "\", \"BFGS\" will be used!"))
    optim.method = "BFGS"
  }

  # OBTAIN MLE AS STARTING VALUE
  if (is.null(beta0)) {
    fit <- eval(call("glm.fit", x =X, y = Y, family = family))
    beta0 <- fit$coef
    if (details) {
      print("Initial beta:")
      print(beta0)
    }
  }
  p <- length(beta0)

  # THE OBJECTIVE FUNCTION
  fun <- ifelse(is.element(family0, c("gaussian", "binomial", "poisson")), LoglikPen, LoglikPenGLM)

  # GRADIENT
  grad <- NULL

  # USING GENERALIZED SIMULATED ANNEALING {GenSA} PLUS BFGS
  if (optim.method == "BFGS") {
    # OPTIMIZATION USING SIMULATED ANNEALING, FOLLOWED BY BFGS
    opt.fit1 <- optim(par=beta0, fn=fun, gr = grad,
                      method = "SANN", control = list(maxit=maxit.global, trace=F, reltol=epsilon),
                      group=group,X=X, y=Y, lambda=lambda, a=a0, family = family0)
    opt.fit2 <- optim(par=opt.fit1$par, fn=fun, gr = grad,
                      method = "BFGS", control = list(maxit=maxit.local, trace=F, reltol=epsilon),
                      group=group, X=X, y=Y, lambda=lambda, a=a0, family = family0)
    gamma <- opt.fit2$par
    # min.Q <- opt.fit2$value
  } else if(optim.method == "GenSA") {
    # THE LOWER AND UPPER BOUNDS FOR SEARCH SPACE IN OPTIMIZATION
    if (is.null(lower)) {lower=rep(-10, p); upper <- rep(10, p)}
    else if (length(lower)==1 && p > 1) {lower <- rep(lower, p); upper <- rep(upper, p)}
    else if (length(lower)!= p) stop("Wrong specification for lower= and upper=. Be aware of its appropriate dimension!")
    opt.fit <- GenSA::GenSA(par=beta0, fn=fun, lower = lower, upper = upper,
                      control=list(maxit=maxit.global, nb.stop.improvement=5),
                      group=group, X=X, y=Y, lambda=lambda, a=a0, family = family0)
    # min.Q <- opt.fit1$value
    gamma <-  opt.fit$par
  } else{
    family.str <- ifelse(is.function(family), family()$family, family$family)
    gamma <- as.vector(adam_gmic(X, Y, a0, lambda, beta0, group, family.str, stepsize, epsilon, maxit.global))
  }
  if (details) {
    print("Fitted Gamma:")
    print(gamma)
  }
  # Reverse the scaling and/or orthogonalization
  if(scale.x || orthogonal.x){
    if (intercept) {
      if(orthogonal.x) {
        gamma[-1] <- unorthogonalize(gamma[-1],X1)
      }
      if(scale.x) {
        gamma[-1] <- gamma[-1]/sc
        gamma[1] <- gamma[1] - sum(gamma[-1]*adj)
      }
    }
    else{
      if(orthogonal.x)  {
        gamma <- unorthogonalize(gamma, X)
      }
      if(scale.x) {
        gamma <- gamma/sc
      }
    }
  }

  # Prepare the variable selection result
  beta.hat <- rep(0, length(gamma))
  grps <- unique(group)
  for(grp in grps) {
    ix <- if(intercept) which(group == grp) + 1 else which(group == grp)
    m = length(ix)
    beta.hat[ix] <- gamma[ix] * tanh(a0 / m *sum(gamma[ix]^2))
  }

  # Obtain SE for gamma
  pvals <- group.pval(formula, data, family, intercept, gamma, group)
  fit0 <- glm(formula, data = data, family = family,
              start=gamma, control=glm.control(epsilon=1e20, maxit=1, trace=F))
  VCOV.gamma <- summary(fit0)$"cov.scaled"
  se.gamma <- sqrt(diag(VCOV.gamma))
  z.gamma <- gamma/se.gamma
  pvalue.gamma <- 2*pnorm(abs(z.gamma), lower.tail=F)

  beta.hat <- round(beta.hat, rounding.digits)  # Rounding

  ##
  result <- cbind(gamma=gamma,se.gamma=se.gamma,pvalue.gamma=pvalue.gamma,beta=beta.hat)
  list(coefficients=data.frame(result), group.pvalues=pvals)
}
