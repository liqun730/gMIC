norm2 <-function(x) sqrt(sum(x^2))

setupLambda <- function(X, y, group, family, lambda.min, nlambda) {
  
  ## Fit to unpenalized covariates
  n <- length(y); intercept <- ncol(X) != length(group)
  fit <- glm(y~1, family=family)

  ## Determine lambda.max
  if (family=="gaussian") {
    r <- fit$residuals
  } else {
    w <- fit$weights
    if (max(w) < 1e-4) stop("Unpenalized portion of model is already saturated; exiting...")
    r <- residuals(fit, "working")*w
  }
  lambda.max <- lambda.min
  for(g in unique(group)){
    ix <- which(group == g)
    z <- norm2(t(X[,ix])%*%r) / (n * sqrt(length(ix)))
    lambda.max <- max(lambda.max, z)
  }
  
  if (lambda.min==0) {
    lambda <- c(exp(seq(log(lambda.max), log(.001*lambda.max), length=nlambda-1)), 0)
  } else {
    lambda <- exp(seq(log(lambda.max), log(lambda.min*lambda.max), length=nlambda))
  }
  rev(lambda)
}


#' The gLasso, gMCP, gSCAD Variable Selection in Generalized Linear Model
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
#' @param penalty What penalty to use, should be one of 'gLASSO', 'gSCAD' and 'gMCP'.
#' @param optim.method Optimization method for non-converx group variable selection, one of c("GD", "GCD", "ADAM").
#' indicating we use, gradient descent, ADAM, or group coordinate descent for grpreg optimmization. Default is "GD".
#' @param scale.x A boolean indicating whether or not to studentize the features. Default is \code{TRUE}.
#' @param orthogonal.x A boolean indicating whether or not to orthogonalize the features within each group. Default is true. See \code{link{orthogonalize}} for details.
#' @param rounding.digits Number of digits after the decimal point for rounding-up estiamtes. Default value is 4.
#' @param maxit.global  Maximum number of iterations allowed for the global optimization algorithm \code{SANN}. Default value is 100.
#' @param maxit.local Maximum number of iterations allowed for the local optimizaiton algorithm \code{BFGS}. Default value is 100.
#' @param epsilon The convergence tolerance.
#' @param stepsize The stepsize (or learning rate) for optim.method = "GD".
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
#' fit <- glm.grpreg(y~.-1,group=group,family="binomial",data=Xy,penalty='gMCP',orthogonal.x=T)
#'
#' @export
#'
glm.grpreg <- function(formula, family = c("gaussian", "binomial", "poisson"), data, group=NULL,
                     beta0=NULL, penalty = c("gLasso", "gMCP", "gSCAD"),
                     optim.method = "GD",
                     scale.x=FALSE, orthogonal.x=FALSE, # Do we orthogonalize X?
                     rounding.digits = 4,				  # Rounding digits for final result
                     maxit.global=100, maxit.local=100, epsilon=1e-6, 
                     stepsize = 0.01, details=FALSE
)
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
  if (details) print(X[1:10, ])  ################
  
  #
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
  
  # Standardize Y
  family.str <- ifelse(is.function(family), family()$family, family$family)
  YY <- if(family.str == "gaussian" & intercept) Y - mean(Y) else Y
  
  # DETERMINE lambdas
  nlambda=100; lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}
  log.lambda = TRUE
  alpha=1
  lambdas <- setupLambda(X, YY, group, family.str, lambda.min, nlambda)
  
  # OBTAIN MLE AS STARTING VALUE
  if (is.null(beta0)) {
    fit <- eval(call("glm.fit", x =X, y = Y, family = family))
    beta0 <- fit$coef
    if (details) print(beta0)
  }
  p <- length(beta0)
  
  # Solving the optimization
  a <- ifelse(penalty=="gSCAD", 4, 3)
  beta <- c(); loss <- df <- aic <- bic <- c()
  if(!optim.method %in% c('GD', 'ADAM', 'GCD')){
    warning(paste0("Unknow optimization method: \"", optim.method, "\", \"GD\" will be used!"))
    optim.method = "GD"
  }
  for(lambda in lambdas){
    bt_lambda <- switch(optim.method, 
                        "GD" = as.vector(gd_grpreg(X, y, beta0, lambda, a, group, family.str, penalty, stepsize, epsilon, maxit.global)),
                        "ADAM" = as.vector(adam_grpreg(X, y, beta0, lambda, a, group, family.str, penalty, stepsize, epsilon, maxit.global)),
                        "GCD" = as.vector(gcd_grpreg(X, y, beta0, lambda, a, group, family.str, penalty, epsilon, maxit.global))
                        )
    bt_lambda <- round(bt_lambda, rounding.digits)
    
    # Reverse the scaling and/or orthogonalization
    if(scale.x || orthogonal.x){
      if (intercept) {
        if(orthogonal.x) bt_lambda[-1] <- unorthogonalize(bt_lambda[-1],X1)
        if(scale.x) {
          bt_lambda[-1] <- bt_lambda[-1]/sc
          bt_lambda[1] <- bt_lambda[1]-sum(bt_lambda[-1]*adj)
        }
      }
      else{
        if(orthogonal.x)  bt_lambda <- unorthogonalize(bt_lambda,X)
        if(scale.x) bt_lambda <- bt_lambda/sc
      }
    }
    
    dfl = sum(bt_lambda != 0)
    ls_lambda = loglik_grpreg(X, y, bt_lambda, family.str)
    beta <- cbind(beta, bt_lambda)
    loss <- c(loss, ls_lambda)
    aic <- c(aic, ls_lambda + 2 * dfl)
    bic <- c(bic, ls_lambda + log(length(y)) * dfl)
    df <- c(df, dfl)
    beta0 <- bt_lambda
  }
  
  val = list()
  val$beta <- beta
  val$lambdas <- lambdas
  val$aic <- aic
  val$bic <- bic
  val$df <- df
  val$loss <- loss
  
  val
}