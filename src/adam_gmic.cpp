#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

//' Gradient Function for gMIC
//' 
//' @param X The design matrix.
//' @param y The response vector
//' @param a The approximation parameter for gMIC.
//' @param lambda The gMIC penalization parameter.
//' @param gamma The optimization parameter.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @export
//'
// [[Rcpp::export]]
arma::vec grad_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family) {
  
  bool intercept = gamma.n_elem != group.n_elem;
  arma::vec grps = unique(group), w, gr, pen;
  int n_group = grps.n_elem, start = 0;
  arma::sp_mat M(gamma.n_elem, gamma.n_elem);
  
  // Construct w and M
  if(intercept) { w.insert_rows(0, arma::ones(1)); M(0, 0) = 1; start += 1; }
  for (int g = 0; g < n_group; g++) {
    arma::uvec ix = (intercept)? arma::uvec(find(group == grps(g)) + 1) : find(group == grps(g));
    arma::vec gamma_g = gamma(ix);
    int m = gamma_g.n_elem;
    double th = std::tanh(a / m * accu(square(gamma_g)));
    M.submat(start, start, start + m - 1, start + m - 1) = (1 - th * th) * 2 * a / m * gamma_g * gamma_g.t() + 
      arma::diagmat(arma::ones(m) * th);
    //double th = std::tanh(a * accu(square(gamma_g)));
    //M.submat(start, start, start + m - 1, start + m - 1) = (1 - th * th) * 2 * a * gamma_g * gamma_g.t();
    w.insert_rows(start, th * arma::ones(m));
    start += m;
  }
  
  arma::vec eta = X * (w % gamma);
  pen = 2 * a * lambda * (1 - square(w)) % gamma;
  
  if(family == "gaussian") {
    gr = 2 * X.n_rows *M * X.t() * (y - eta) / accu(square(y - eta));
  }
  else if(family == "binomial") {
    gr = 2 * M * X.t() * (arma::exp(eta) / (1 + arma::exp(eta)) - y);
  }
  else {
    gr = 2 * M * X.t() * (arma::exp(eta) - y);
  }
  
  return gr + pen;
}

//' The ADAM Algorithm for gMIC Optimization (experimental)
//' 
//' @param X Design matrix.
//' @param y The response vector.
//' @param a The approximation parameter for gMIC.
//' @param lambda The penalization parameter for gMIC, e.g., 2 for AIC and long(n) for BIC.
//' @param gamma The optimization parameter gamma.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @param stepsize Stepsize for group coordinate descent.
//' @param tol Convergence tolerance.
//' @param maxit Maximum number of iterations.
//' @param b1 ADAM hyperparameter, default is 0.9.
//' @param b2 ADAM hyperparameter, default is 0.999.
//' @param e ADAM hyperparameter, default is to 1e-8.
//' @export
//'
// [[Rcpp::export]]
arma::vec adam_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family,
                    double stepsize, double tol, int maxit, double b1 = 0.7, double b2 = 0.9, double e = 1.0e-8) {
  
  // initialization for ADAM:
  arma::vec mt = arma::zeros(gamma.n_elem), vt = arma::zeros(gamma.n_elem);
  arma::uvec ix;
  
  int t = 0;
  while(t < maxit) {
    arma::vec gamma0 = gamma;
    arma::vec grad = grad_gmic(X, y, a, lambda, gamma, group, family);
    mt = (b1 * mt + (1 - b1) * grad) / (1 - std::pow(b1, t+1));
    vt = (b2 * vt + (1 - b2) * grad % grad) / (1 - std::pow(b2, t+1));
    gamma -= stepsize * mt / (arma::sqrt(vt) + e);
    if(sqrt(accu(square(gamma - gamma0))) < tol) t = maxit;
    else t += 1;
  }

  return gamma;
}
