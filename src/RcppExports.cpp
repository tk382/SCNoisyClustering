// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// ic_c
arma::vec ic_c(int k, arma::mat Y, arma::vec phi, arma::rowvec xi, arma::mat theta, int p, int n);
RcppExport SEXP _SCNoisyClustering_ic_c(SEXP kSEXP, SEXP YSEXP, SEXP phiSEXP, SEXP xiSEXP, SEXP thetaSEXP, SEXP pSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type k(kSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type phi(phiSEXP);
    Rcpp::traits::input_parameter< arma::rowvec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type theta(thetaSEXP);
    Rcpp::traits::input_parameter< int >::type p(pSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(ic_c(k, Y, phi, xi, theta, p, n));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_SCNoisyClustering_ic_c", (DL_FUNC) &_SCNoisyClustering_ic_c, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_SCNoisyClustering(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
