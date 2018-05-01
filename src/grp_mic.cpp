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
  int n_group = grps.n_elem, count=0;
  arma::vec mm = arma::zeros(n_group), ww = arma::zeros(n_group), w;
  double pen;
  for(int i=0; i<n_group; i++){
    arma::uvec ix = find(group == grps(i));
    mm(i) = ix.n_elem;
    ww(i) = std::tanh(a/mm(i)*accu(square(beta(ix))));
    w.insert_rows(count,ww(i)*arma::ones(mm(i)));
    count += mm(i);
  }
  pen = sum(w);
  return Rcpp::List::create(Rcpp::Named("w") = w,Rcpp::Named("pen") = pen);
}
