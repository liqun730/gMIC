// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// grad_gmic
arma::vec grad_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family);
RcppExport SEXP _gMIC_grad_gmic(SEXP XSEXP, SEXP ySEXP, SEXP aSEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP groupSEXP, SEXP familySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    rcpp_result_gen = Rcpp::wrap(grad_gmic(X, y, a, lambda, gamma, group, family));
    return rcpp_result_gen;
END_RCPP
}
// adam_gmic
arma::vec adam_gmic(arma::mat X, arma::vec y, double a, double lambda, arma::vec gamma, arma::vec group, String family, double stepsize, double tol, int maxit, double b1, double b2, double e);
RcppExport SEXP _gMIC_adam_gmic(SEXP XSEXP, SEXP ySEXP, SEXP aSEXP, SEXP lambdaSEXP, SEXP gammaSEXP, SEXP groupSEXP, SEXP familySEXP, SEXP stepsizeSEXP, SEXP tolSEXP, SEXP maxitSEXP, SEXP b1SEXP, SEXP b2SEXP, SEXP eSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< String >::type family(familySEXP);
    Rcpp::traits::input_parameter< double >::type stepsize(stepsizeSEXP);
    Rcpp::traits::input_parameter< double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< int >::type maxit(maxitSEXP);
    Rcpp::traits::input_parameter< double >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< double >::type b2(b2SEXP);
    Rcpp::traits::input_parameter< double >::type e(eSEXP);
    rcpp_result_gen = Rcpp::wrap(adam_gmic(X, y, a, lambda, gamma, group, family, stepsize, tol, maxit, b1, b2, e));
    return rcpp_result_gen;
END_RCPP
}
// pen_gmic
List pen_gmic(arma::vec beta, arma::vec group, double a);
RcppExport SEXP _gMIC_pen_gmic(SEXP betaSEXP, SEXP groupSEXP, SEXP aSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type group(groupSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    rcpp_result_gen = Rcpp::wrap(pen_gmic(beta, group, a));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_gMIC_grad_gmic", (DL_FUNC) &_gMIC_grad_gmic, 7},
    {"_gMIC_adam_gmic", (DL_FUNC) &_gMIC_adam_gmic, 13},
    {"_gMIC_pen_gmic", (DL_FUNC) &_gMIC_pen_gmic, 3},
    {NULL, NULL, 0}
};

RcppExport void R_init_gMIC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
