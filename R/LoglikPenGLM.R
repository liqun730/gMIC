#' Compute the penalized log likelihood for GLM with MIC penalty via R function \code{glm}
#'
#' @param beta A p-dimensional vector containing the regression ceofficients.
#' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
#'   It means the first 2 features form a group of variables and the last 2 features form another group of variables.
#' @param X An \eqn{n} by \eqn{p} design matrix.
#' @param y The \eqn{n} by 1 response vector
#' @param lambda The penalty parameter euqals either 2 in AIC or ln(n) in BIC (by default).
#' It can be specified as any value of the user's own choice.
#' @param a The scale parameter in the hyperbolic tangent function of the MIC penalty. By default, \eqn{a = 50}.
#' @param family a description of the error distribution and link function to be used in the model. It needs to be the result
#' of a call to a family function since \code{glm.fit} is used here. In other words, it can NOT be specified as a character
#' string naming a family function nor a family function. See \code{\link[stats]{glm}} and \code{\link[stats]{family}} for
#' details of family functions.
#' @return The value of the penalized log likelihood function evaluated at beta.
#' @seealso \code{\link[stats]{glm}}, \code{\link[stats]{family}}
#'
#' @export
#'

LoglikPenGLM <- function(beta, group, X, y, lambda, a, family = gaussian(link = "identity"))
{
  # THE PANALTY PART
  intercept <- length(beta) != length(group)
  bt <- if(intercept) beta[-1] else beta
  tmp <- .Call('_gMIC_pen_gmic', PACKAGE = 'gMIC', bt, group, a)
  w <- tmp$w; pen <- lambda*tmp$pen
  beta.prime <- bt*(w)
  if(intercept) beta.prime <- c(beta[1],beta.prime)

  # THE LOGLIKELIHOOD PART OBTAINED FROM R FUNCTION glm()
  eta <- X%*%beta.prime; n <- NROW(X);
  fit <- eval(call("glm.fit", x = matrix(, n, 0L), y = y, offset = eta,
                   family = family,
                   control = glm.control(epsilon = 1e+8, maxit = 0.5, trace = T),
                   intercept = F))
  dev <- fit$deviance

  # OUTPUT THE PENALZIED LOGLIKELIHOOD
  loglik.pen <- dev + pen
  return(loglik.pen)
}
