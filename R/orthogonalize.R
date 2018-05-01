#' Orthogonalize the Design Matrix with Specified Group Structure
#'
#' @param X The design matrix to be orthogonalized.
#' @param group The group structure based on which X is orthogonalized. For example, assume that X has 4 columns and group=c(1,1,2,2).
#'     It means the first 2 features form a group of variables and the last 2 features form another group of variables.
#' @return The orthogonalized design. Assume that X_g denotes columns corresponding to the g-th group, then t(X_g)%*%(X_g) = I.
#' @examples
#' X <- matrix(runif(1200),100,12)
#' group <- rep(1:3,each=4)
#' X <- orthogonalize(X,group)
#'
#' @export
#'
orthogonalize <- function(X, group) {
  n <- nrow(X)
  J <- max(group)
  T1 <- vector("list", J)
  X.tmp <- X
  for (j in seq_along(numeric(J))) {
    ind <- which(group==j)
    SVD <- svd(X[, ind, drop=FALSE], nu=0)
    T1[[j]] <- sweep(SVD$v, 2, sqrt(n)/SVD$d, "*")
    X[,ind] <- X.tmp[,ind]%*%T1[[j]]
  }
  attr(X, "T1") <- T1
  X
}


# Unorthogonalize fitted beta within each group
unorthogonalize <- function(b, XX) {
  T1 <- Matrix::bdiag(attr(XX, "T1"))
  val <- as.vector(T1 %*% b)
  return(val)
}
