#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


//' The -2*Log Likelihood Function for GLM
//' 
//' @param X The design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @export
//'
// [[Rcpp::export]]
double loglik_grpreg(arma::mat X, arma::vec y, arma::vec beta, String family) {
  
  arma::vec eta = X * beta;
  double l;
  
  if(family == "gaussian") {
    l = accu(square(y - eta));
  }
  else if(family == "binomial") {
    l = -2*accu(y % eta - arma::log(1+arma::exp(eta)));
  }
  else {
    l = -2*accu(y % eta - arma::exp(eta));
  }
  return l;

}


//' The Gradient Function for GLM
//' 
//' @param X The design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
arma::vec grad_loglik(arma::mat X, arma::vec y, arma::vec beta, String family) {
  
  arma::vec gr;
  arma::vec eta = X * beta;

  if(family == "gaussian") {
    gr = -X.t() * (y - eta) / 2;
  }
  else if(family == "binomial") {
    gr = X.t() * (y - arma::exp(eta) / (1 + arma::exp(eta)));
  }
  else {
    gr = X.t() * (y - arma::exp(eta));
  }
  
  return gr / X.n_rows;
}


//' The Gradient Function for gLASSO, gSCAD, gMCP penalties
//' 
//' @param beta The model coefficients.
//' @param lambda The gMIC penalization parameter.
//' @param a The hyper parameter of MCP or SCAD.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param penalty One of "gLASSO", "gSCAD" and "gMCP"
arma::vec grad_pen(arma::vec beta, double lambda, double a, arma::vec group, String penalty){
  
  arma::vec grps = unique(group), gr;
  int n_group = grps.n_elem, start = 0;
  bool intercept = beta.n_elem != group.n_elem;
  
  if(intercept) { gr.insert_rows(0, zeros(1)); start += 1; }
  for (int g = 0; g < n_group; g++) {
    arma::uvec ix = (intercept)? arma::uvec(find(group == grps(g)) + 1) : find(group == grps(g));
    arma::vec beta_g = beta(ix);
    int m = beta_g.n_elem;
    double gnorm = sqrt(accu(square(beta_g)));
    if (penalty == "gLASSO"){
      gr.insert_rows(start, lambda * sqrt(m) * beta_g);
    } else if (penalty == "gMCP") {
      if (gnorm == 0) gr.insert_rows(start, lambda * sqrt(m) * ones(m));
      else if (gnorm <= a * lambda * sqrt(m)) gr.insert_rows(start, (lambda / gnorm - 1.0 / a) * beta_g);
      else gr.insert_rows(start, zeros(m));
    } else {
      if (gnorm == 0) gr.insert_rows(start, lambda * sqrt(m) * ones(m));
      else if(gnorm <= lambda * sqrt(m)) gr.insert_rows(start, lambda * beta_g / gnorm);
      else if(gnorm <= a * lambda * sqrt(m)) gr.insert_rows(start, (a * lambda * sqrt(m) - gnorm) / (gnorm * (a - 1.0)) * beta_g);
      else gr.insert_rows(start, zeros(m));
    }
    start += m;
  }

  return gr;
}


//' The Gradient Function for GLM with gLASSO, gSCAD, gMCP penalties
//' 
//' @param X The design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param lambda The gMIC penalization parameter.
//' @param a The hyper parameter of MCP or SCAD.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @param penalty One of "gLASSO", "gSCAD" and "gMCP".
arma::vec grad_grpreg(arma::mat X, arma::vec y, arma::vec beta, double lambda, double a, arma::vec group, String family, String penalty) {
  arma::vec gr = grad_loglik(X, y, beta, family);
  arma::vec pen = grad_pen(beta, lambda, a, group, penalty);
  return -2 * gr + pen;
}

//' The Gradient Descent Algorithm for GLM with gLASSO, gSCAD, gMCP penalties
//' 
//' @param X Design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param lambda The penalization parameter for gMIC, e.g., 2 for AIC and long(n) for BIC.
//' @param a The approximation parameter for MCP or SCAD.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @param penalty One of "gLASSO", "gSCAD" and "gMCP".
//' @param stepsize Stepsize for group coordinate descent.
//' @param tol Convergence tolerance.
//' @param maxit Maximum number of iterations.
//' @export
//'
// [[Rcpp::export]]
arma::vec gd_grpreg(arma::mat X, arma::vec y, arma::vec beta, double lambda, double a, arma::vec group, String family, 
                  String penalty, double stepsize, double tol, int maxit) {
  
  int k = 0;
  // double thresh = 1.0e-5;
  while(k < maxit) {
    arma::vec beta0 = beta;
    arma::vec grad = grad_grpreg(X, y, beta, lambda, a, group, family, penalty);
    beta -= stepsize * beta;
    // arma::uvec ix = find(gamma < thresh && gamma > -thresh);
    // gamma.elem(ix).fill(0.0);
    if(sqrt(accu(square(beta - beta0))) < tol) k = maxit;
    else k += 1;
  }
  
  return beta;
}



//' The ADAM Algorithm for GLM with gLASSO, gSCAD, gMCP penalties (experimental)
//' 
//' @param X Design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param lambda The penalization parameter for gMIC, e.g., 2 for AIC and long(n) for BIC.
//' @param a The approximation parameter for MCP or SCAD.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @param penalty One of "gLASSO", "gSCAD" and "gMCP".
//' @param stepsize Stepsize for group coordinate descent.
//' @param tol Convergence tolerance.
//' @param maxit Maximum number of iterations.
//' @param b1 ADAM hyperparameter, default set to 0.9.
//' @param b2 ADAM hyperparameter, default set to 0.999.
//' @param e ADAM hyperparameter, default set to 1e-8.
//' @export
//'
// [[Rcpp::export]]
arma::vec adam_grpreg(arma::mat X, arma::vec y, arma::vec beta, double lambda, double a, arma::vec group, String family,
                      String penalty, double stepsize, double tol, int maxit, double b1 = 0.7, double b2 = 0.9, double e = 1.0e-8) {
  
  // initialization for ADAM:
  arma::vec mt = arma::zeros(beta.n_elem), vt = arma::zeros(beta.n_elem);
  // double thresh = 1.0e-5;
  
  int t = 0;
  while(t < maxit) {
    arma::vec beta0 = beta;
    arma::vec grad = grad_grpreg(X, y, beta, lambda, a, group, family, penalty);
    mt = (b1 * mt + (1 - b1) * grad) / (1 - std::pow(b1, t+1));
    vt = (b2 * vt + (1 - b2) * grad % grad) / (1 - std::pow(b2, t+1));
    beta -= stepsize * mt / (arma::sqrt(vt) + e);
    // arma::uvec ix = find(gamma < thresh && gamma > -thresh);
    // gamma.elem(ix).fill(0.0);
    if(sqrt(accu(square(beta - beta0))) < tol) t = maxit;
    else t += 1;
  }
  
  return beta;
}


// Soft thresholding
arma::vec soft(arma::vec beta_g, double t) {
  double gnorm = sqrt(accu(square(beta_g)));
  if (gnorm <= t) return arma::zeros(beta_g.n_elem);
  return (1 - t / gnorm) * beta_g;
}


// Proximal
arma::vec proximal(arma::vec beta_g, double lambda, double a, String penalty) {
  
  arma::vec res;
  double gnorm = sqrt(accu(square(beta_g)));
  
  if (penalty == "gLASSO") {
    res = soft(beta_g, lambda);
  } else if (penalty == "gMCP") {
    if (gnorm <= a * lambda) res = (a / (a - 1)) * soft(beta_g, lambda);
    else res = beta_g;
  } else {
    if (gnorm <= 2 * lambda) res = soft(beta_g, lambda);
    else if (gnorm <= a * lambda) res = ((a - 1) / (a - 2)) * soft(beta_g, a * lambda / (a - 1));
    else res = beta_g;
  }
  
  return res;
}


//' The Group Coordinate Descent Algorithm for GLM with gLASSO, gSCAD, gMCP penalties
//' 
//' @param X Orthonormalized design matrix.
//' @param y The response vector.
//' @param beta The model coefficients.
//' @param lambda The penalization parameter for gMIC, e.g., 2 for AIC and long(n) for BIC.
//' @param a The approximation parameter for gMIC.
//' @param group The group structure of the model. For example, assume that X has 4 columns and group=c(1,1,2,2).
//' It means the first 2 features form a group of variables and the last 2 features form another group of variables.
//' @param family The type of glm model, should be one of "gaussian", "binomial" or "poisson".
//' @param penalty One of "gLASSO", "gSCAD" and "gMCP".
//' @param maxit Maximum number of iterations.
//' @export
//'
// [[Rcpp::export]]
arma::vec gcd_grpreg(arma::mat X, arma::vec y, arma::vec beta, double lambda, double a, arma::vec group,
                     String family, String penalty, double tol, int maxit) {
  
  
  arma::vec grps = unique(group);
  int n_group = grps.n_elem, t = 0;
  bool intercept = X.n_cols != group.n_elem;
  
  while(t < maxit) {
    
    arma::vec beta0 = beta;
    
    for (int g = 0; g < n_group; g++) {
      arma::uvec ixg = (intercept)? arma::uvec(find(group == grps(g)) + 1) : find(group == grps(g));
      arma::vec betag = beta(ixg);
      arma::mat Xg = X.cols(ixg);
      int mg = betag.n_elem;
      if (family == "gaussian") {
        arma::vec vg = y - X * beta + Xg * betag;
        betag = proximal(Xg.t() * vg / X.n_rows, lambda * sqrt(mg), a, penalty);
      } else {
        arma::vec Xbeta = X * beta;
        arma::vec hgs = (family == "binomial")? arma::vec(arma::exp(Xbeta) / square(1 + arma::exp(Xbeta))) : square(arma::exp(Xbeta));
        double hg = 1e-9;
        for (int i = 0; i < mg; i++) {
          double h = accu(square(Xg.col(i)) % hgs);
          if (h > hg) hg = h;
        }
        arma::vec grad = (family == "binomial")? arma::vec(Xg.t() * (y - arma::exp(Xbeta) / (1 + arma::exp(Xbeta)))) :
          Xg.t() * (y - arma::exp(Xbeta));
        betag = proximal(grad / (X.n_rows * hg), lambda * sqrt(mg) / hg, a, penalty);
      }
      beta(ixg) = betag;
    }
    
    if(sqrt(accu(square(beta - beta0))) < tol) t = maxit;
    else t += 1;
    
  }
  
  return beta;
}
