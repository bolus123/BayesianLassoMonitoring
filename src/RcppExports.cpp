// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// rinvgaussiancpp
arma::colvec rinvgaussiancpp(int n, double mu, double lambda);
RcppExport SEXP _BayesianLassoMonitoring_rinvgaussiancpp(SEXP nSEXP, SEXP muSEXP, SEXP lambdaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< double >::type lambda(lambdaSEXP);
    rcpp_result_gen = Rcpp::wrap(rinvgaussiancpp(n, mu, lambda));
    return rcpp_result_gen;
END_RCPP
}
// getV
arma::mat getV(arma::colvec Y, int q);
RcppExport SEXP _BayesianLassoMonitoring_getV(SEXP YSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(getV(Y, q));
    return rcpp_result_gen;
END_RCPP
}
// rmvnorm
arma::colvec rmvnorm(arma::colvec Mean, arma::mat Sigma);
RcppExport SEXP _BayesianLassoMonitoring_rmvnorm(SEXP MeanSEXP, SEXP SigmaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec >::type Mean(MeanSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Sigma(SigmaSEXP);
    rcpp_result_gen = Rcpp::wrap(rmvnorm(Mean, Sigma));
    return rcpp_result_gen;
END_RCPP
}
// rtwosegnorm
arma::colvec rtwosegnorm(int n, double a, double b, double mean, double sd);
RcppExport SEXP _BayesianLassoMonitoring_rtwosegnorm(SEXP nSEXP, SEXP aSEXP, SEXP bSEXP, SEXP meanSEXP, SEXP sdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< double >::type a(aSEXP);
    Rcpp::traits::input_parameter< double >::type b(bSEXP);
    Rcpp::traits::input_parameter< double >::type mean(meanSEXP);
    Rcpp::traits::input_parameter< double >::type sd(sdSEXP);
    rcpp_result_gen = Rcpp::wrap(rtwosegnorm(n, a, b, mean, sd));
    return rcpp_result_gen;
END_RCPP
}
// getHMatMT
arma::mat getHMatMT(int T, int q);
RcppExport SEXP _BayesianLassoMonitoring_getHMatMT(SEXP TSEXP, SEXP qSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    rcpp_result_gen = Rcpp::wrap(getHMatMT(T, q));
    return rcpp_result_gen;
END_RCPP
}
// getHMatSustained
arma::mat getHMatSustained(int T, int q, int w);
RcppExport SEXP _BayesianLassoMonitoring_getHMatSustained(SEXP TSEXP, SEXP qSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(getHMatSustained(T, q, w));
    return rcpp_result_gen;
END_RCPP
}
// getHMatIsolated
arma::mat getHMatIsolated(int T, int q, int w);
RcppExport SEXP _BayesianLassoMonitoring_getHMatIsolated(SEXP TSEXP, SEXP qSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(getHMatIsolated(T, q, w));
    return rcpp_result_gen;
END_RCPP
}
// getHMatGradual
arma::mat getHMatGradual(int T, int q, int w);
RcppExport SEXP _BayesianLassoMonitoring_getHMatGradual(SEXP TSEXP, SEXP qSEXP, SEXP wSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type T(TSEXP);
    Rcpp::traits::input_parameter< int >::type q(qSEXP);
    Rcpp::traits::input_parameter< int >::type w(wSEXP);
    rcpp_result_gen = Rcpp::wrap(getHMatGradual(T, q, w));
    return rcpp_result_gen;
END_RCPP
}
// GibbsRFLSM
Rcpp::List GibbsRFLSM(arma::colvec& Y, int& q, arma::mat& A, double& a, double& b, double& alpha, double& beta, double& theta1, double& theta2, double& xi2, Rcpp::String& method, double& bound0, double& boundqplus1, int& nsim, int& by, int& burnin, double& tol, Rcpp::Nullable<Rcpp::NumericMatrix> H);
RcppExport SEXP _BayesianLassoMonitoring_GibbsRFLSM(SEXP YSEXP, SEXP qSEXP, SEXP ASEXP, SEXP aSEXP, SEXP bSEXP, SEXP alphaSEXP, SEXP betaSEXP, SEXP theta1SEXP, SEXP theta2SEXP, SEXP xi2SEXP, SEXP methodSEXP, SEXP bound0SEXP, SEXP boundqplus1SEXP, SEXP nsimSEXP, SEXP bySEXP, SEXP burninSEXP, SEXP tolSEXP, SEXP HSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::colvec& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< int& >::type q(qSEXP);
    Rcpp::traits::input_parameter< arma::mat& >::type A(ASEXP);
    Rcpp::traits::input_parameter< double& >::type a(aSEXP);
    Rcpp::traits::input_parameter< double& >::type b(bSEXP);
    Rcpp::traits::input_parameter< double& >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< double& >::type theta1(theta1SEXP);
    Rcpp::traits::input_parameter< double& >::type theta2(theta2SEXP);
    Rcpp::traits::input_parameter< double& >::type xi2(xi2SEXP);
    Rcpp::traits::input_parameter< Rcpp::String& >::type method(methodSEXP);
    Rcpp::traits::input_parameter< double& >::type bound0(bound0SEXP);
    Rcpp::traits::input_parameter< double& >::type boundqplus1(boundqplus1SEXP);
    Rcpp::traits::input_parameter< int& >::type nsim(nsimSEXP);
    Rcpp::traits::input_parameter< int& >::type by(bySEXP);
    Rcpp::traits::input_parameter< int& >::type burnin(burninSEXP);
    Rcpp::traits::input_parameter< double& >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<Rcpp::NumericMatrix> >::type H(HSEXP);
    rcpp_result_gen = Rcpp::wrap(GibbsRFLSM(Y, q, A, a, b, alpha, beta, theta1, theta2, xi2, method, bound0, boundqplus1, nsim, by, burnin, tol, H));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_BayesianLassoMonitoring_rinvgaussiancpp", (DL_FUNC) &_BayesianLassoMonitoring_rinvgaussiancpp, 3},
    {"_BayesianLassoMonitoring_getV", (DL_FUNC) &_BayesianLassoMonitoring_getV, 2},
    {"_BayesianLassoMonitoring_rmvnorm", (DL_FUNC) &_BayesianLassoMonitoring_rmvnorm, 2},
    {"_BayesianLassoMonitoring_rtwosegnorm", (DL_FUNC) &_BayesianLassoMonitoring_rtwosegnorm, 5},
    {"_BayesianLassoMonitoring_getHMatMT", (DL_FUNC) &_BayesianLassoMonitoring_getHMatMT, 2},
    {"_BayesianLassoMonitoring_getHMatSustained", (DL_FUNC) &_BayesianLassoMonitoring_getHMatSustained, 3},
    {"_BayesianLassoMonitoring_getHMatIsolated", (DL_FUNC) &_BayesianLassoMonitoring_getHMatIsolated, 3},
    {"_BayesianLassoMonitoring_getHMatGradual", (DL_FUNC) &_BayesianLassoMonitoring_getHMatGradual, 3},
    {"_BayesianLassoMonitoring_GibbsRFLSM", (DL_FUNC) &_BayesianLassoMonitoring_GibbsRFLSM, 18},
    {NULL, NULL, 0}
};

RcppExport void R_init_BayesianLassoMonitoring(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
