#' Compute the penalized log likelihood for GLM with MIC penalty (Self-Written)
#'
#' @param beta A p-dimensional vector containing the regression ceofficients.
#' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
#'   It means the first 2 features form a group of variables and the last 2 features form another group of variables.
#' @param X An \eqn{n} by \eqn{p} design matrix.
#' @param y The \eqn{n} by 1 response vector
#' @param lambda The penalty parameter euqals either 2 in AIC or ln(n) in BIC (by default).
#' It can be specified as any value of the user's own choice.
#' @param a The scale parameter in the hyperbolic tangent function of the MIC penalty. By default, \eqn{a = 50}.
#' @param family a description of the error distribution. To use this function, \code{family} needs to be one of \code{"gaussian"},
#'
#' @return The value of the penalized log likelihood function evaluated at beta.
#' @details
#' This function is much faster than \code{\link{LoglikPenGLM}}, but it is only applicable for Gaussian linear regression, logistic regression,
#' and loglinear or Poisson regression models. To take advantage, sepcify  \code{family} as \code{family="gaussian"}, \code{family="binomial"}, or \code{family="poisson"} only.
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#'
#' @export


LoglikPen <- function(beta, group, X, y, lambda, a, family = "gaussian")
{

  # THE PANALTY PART
  intercept <- length(beta) != length(group)
  bt <- if(intercept) beta[-1] else beta
  tmp <- .Call('_gMIC_grp_mic', PACKAGE = 'gMIC', bt, group, a)
  w <- tmp$w; pen <- lambda*tmp$pen
  beta.prime <- bt*(w)
  if(intercept) beta.prime <- c(beta[1],beta.prime)

  # THE (-2)*LOGLIKELIHOOD PART
  eta <- X%*%beta.prime
  if (family=="gaussian")
    loglik <- NROW(X)* log(sum((y - eta)^2))
  else if (family=="binomial")
    loglik <- -2*sum(y*eta - log(1+exp(eta)))
  else if (family=="poisson")
    loglik <- -2*sum(y*eta - exp(eta))
  else stop("Hmmm. How did you get here? Wrong specification of family.")

  # OUTPUT THE PENALZIED LOGLIKELIHOOD OR SMOOTHED AIC/BIC
  loglik.pen <- loglik + pen;
  return(loglik.pen)
}
