#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' The Group MIC Penalty Function
//'
//' @param beta A p-dimensional vector containing the regression ceofficients.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param a The approximation parameter of the gMIC.
//' @return The gMIC penalty evaluated at beta.
//' @export
//'
// [[Rcpp::export]]
List grp_mic(arma::vec beta, arma::vec group, double a){
  arma::vec grps = unique(group);
  int n_group = grps.n_elem, count = 0;
  arma::vec w;
  double pen;
  for(int i = 0; i < n_group; i++){
    arma::uvec ix = find(group == grps(i));
    int m = ix.n_elem;
    double th = std::tanh(a / m * accu(square(beta(ix))));
    w.insert_rows(count, th * arma::ones(m));
    count += m;
  }
  pen = sum(w);
  return Rcpp::List::create(Rcpp::Named("w") = w,Rcpp::Named("pen") = pen);
}


//' The gradient function for gMIC
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
  if(intercept) { w.insert_rows(0, 1); M(0, 0) = 1; start += 1; }
  for (int g = 0; g < n_group; g++) {
    arma::uvec ix = (intercept)? arma::uvec(find(group == grps(g)) + 1) : find(group == grps(g));
    arma::vec gamma_g = gamma(ix);
    int m = gamma_g.n_elem;
    double th = std::tanh(a / m * accu(square(gamma_g)));
    M.submat(start, start, start + m - 1, start + m - 1) = (1 - th * th) * 2 * a / m * gamma_g * gamma_g.t();
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



//' The gradient descent algorithm for gMIC optimization
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
//' @export
//'
// [[Rcpp::export]]
arma::vec gd_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family, 
                  double stepsize, double tol, int maxit) {

  int k = 0;
  // double thresh = 1.0e-5;
  while(k < maxit) {
    arma::vec gamma0 = gamma;
    arma::vec grad = grad_gmic(X, y, a, lambda, gamma, group, family);
    gamma -= stepsize * grad;
    // arma::uvec ix = find(gamma < thresh && gamma > -thresh);
    // gamma.elem(ix).fill(0.0);
    if(sqrt(accu(square(gamma - gamma0))) < tol) k = maxit;
    else k += 1;
  }
  
  return gamma;
}



//' The ADAM algorithm for gMIC optimization (experimental)
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
//' @param b1 ADAM hyperparameter, default set to 0.9.
//' @param b2 ADAM hyperparameter, default set to 0.999.
//' @param e ADAM hyperparameter, default set to 1e-8.
//' @export
//'
// [[Rcpp::export]]
arma::vec adam_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family,
                    double stepsize, double tol, int maxit, double b1 = 0.7, double b2 = 0.9, double e = 1.0e-8) {
  
  // initialization for ADAM:
  arma::vec mt = arma::zeros(gamma.n_elem), vt = arma::zeros(gamma.n_elem);
  // double thresh = 1.0e-5;
  
  int t = 0;
  while(t < maxit) {
    arma::vec gamma0 = gamma;
    arma::vec grad = grad_gmic(X, y, a, lambda, gamma, group, family);
    mt = (b1 * mt + (1 - b1) * grad) / (1 - std::pow(b1, t+1));
    vt = (b2 * vt + (1 - b2) * grad % grad) / (1 - std::pow(b2, t+1));
    gamma -= stepsize * mt / (arma::sqrt(vt) + e);
    // arma::uvec ix = find(gamma < thresh && gamma > -thresh);
    // gamma.elem(ix).fill(0.0);
    if(sqrt(accu(square(gamma - gamma0))) < tol) t = maxit;
    else t += 1;
  }
  
  return gamma;
}
