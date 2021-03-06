% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/LoglikPen.R
\name{LoglikPen}
\alias{LoglikPen}
\title{Compute the penalized log likelihood for GLM with MIC penalty (Self-Written)}
\usage{
LoglikPen(gamma, group, X, y, lambda, a, family = "gaussian")
}
\arguments{
\item{group}{The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
It means the first 2 features form a group of variables and the last 2 features form another group of variables.}

\item{X}{An \eqn{n} by \eqn{p} design matrix.}

\item{y}{The \eqn{n} by 1 response vector}

\item{lambda}{The penalty parameter euqals either 2 in AIC or ln(n) in BIC (by default).
It can be specified as any value of the user's own choice.}

\item{a}{The scale parameter in the hyperbolic tangent function of the MIC penalty. By default, \eqn{a = 50}.}

\item{family}{a description of the error distribution. To use this function, \code{family} needs to be one of \code{"gaussian"},}

\item{beta}{A p-dimensional vector containing the regression ceofficients.}
}
\value{
The value of the penalized log likelihood function evaluated at beta.
}
\description{
Compute the penalized log likelihood for GLM with MIC penalty (Self-Written)
}
\details{
This function is much faster than \code{\link{LoglikPenGLM}}, but it is only applicable for Gaussian linear regression, logistic regression,
and loglinear or Poisson regression models. To take advantage, sepcify  \code{family} as \code{family="gaussian"}, \code{family="binomial"}, or \code{family="poisson"} only.
}
\seealso{
\code{\link[stats]{glm}}, \code{\link[stats]{family}}
}
