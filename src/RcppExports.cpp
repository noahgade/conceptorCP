// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// NRMSE
double NRMSE(arma::mat output, arma::mat target);
RcppExport SEXP _conceptorCP_NRMSE(SEXP outputSEXP, SEXP targetSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type output(outputSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type target(targetSEXP);
    rcpp_result_gen = Rcpp::wrap(NRMSE(output, target));
    return rcpp_result_gen;
END_RCPP
}
// runRNN
arma::cube runRNN(arma::mat input, arma::cube W, float RSInitial, float bscale, float iscale);
RcppExport SEXP _conceptorCP_runRNN(SEXP inputSEXP, SEXP WSEXP, SEXP RSInitialSEXP, SEXP bscaleSEXP, SEXP iscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type RSInitial(RSInitialSEXP);
    Rcpp::traits::input_parameter< float >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< float >::type iscale(iscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(runRNN(input, W, RSInitial, bscale, iscale));
    return rcpp_result_gen;
END_RCPP
}
// conceptorCalc
arma::cube conceptorCalc(arma::cube RS, float aperture);
RcppExport SEXP _conceptorCP_conceptorCalc(SEXP RSSEXP, SEXP apertureSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type RS(RSSEXP);
    Rcpp::traits::input_parameter< float >::type aperture(apertureSEXP);
    rcpp_result_gen = Rcpp::wrap(conceptorCalc(RS, aperture));
    return rcpp_result_gen;
END_RCPP
}
// runCRNN
arma::cube runCRNN(arma::mat input, arma::cube W, arma::cube C, float RSInitial, int washoutL, float bscale, float iscale);
RcppExport SEXP _conceptorCP_runCRNN(SEXP inputSEXP, SEXP WSEXP, SEXP CSEXP, SEXP RSInitialSEXP, SEXP washoutLSEXP, SEXP bscaleSEXP, SEXP iscaleSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type C(CSEXP);
    Rcpp::traits::input_parameter< float >::type RSInitial(RSInitialSEXP);
    Rcpp::traits::input_parameter< int >::type washoutL(washoutLSEXP);
    Rcpp::traits::input_parameter< float >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< float >::type iscale(iscaleSEXP);
    rcpp_result_gen = Rcpp::wrap(runCRNN(input, W, C, RSInitial, washoutL, bscale, iscale));
    return rcpp_result_gen;
END_RCPP
}
// WoutCalc
arma::cube WoutCalc(arma::mat input, arma::cube ResStates, float regular);
RcppExport SEXP _conceptorCP_WoutCalc(SEXP inputSEXP, SEXP ResStatesSEXP, SEXP regularSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type ResStates(ResStatesSEXP);
    Rcpp::traits::input_parameter< float >::type regular(regularSEXP);
    rcpp_result_gen = Rcpp::wrap(WoutCalc(input, ResStates, regular));
    return rcpp_result_gen;
END_RCPP
}
// outputCalc
arma::cube outputCalc(arma::cube Wout, arma::cube CResStates);
RcppExport SEXP _conceptorCP_outputCalc(SEXP WoutSEXP, SEXP CResStatesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type Wout(WoutSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type CResStates(CResStatesSEXP);
    rcpp_result_gen = Rcpp::wrap(outputCalc(Wout, CResStates));
    return rcpp_result_gen;
END_RCPP
}
// RNNParamFit
arma::mat RNNParamFit(arma::mat input, arma::cube W, float RSInitial, int washoutL, int trainL, arma::vec bscales, arma::vec iscales);
RcppExport SEXP _conceptorCP_RNNParamFit(SEXP inputSEXP, SEXP WSEXP, SEXP RSInitialSEXP, SEXP washoutLSEXP, SEXP trainLSEXP, SEXP bscalesSEXP, SEXP iscalesSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type input(inputSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W(WSEXP);
    Rcpp::traits::input_parameter< float >::type RSInitial(RSInitialSEXP);
    Rcpp::traits::input_parameter< int >::type washoutL(washoutLSEXP);
    Rcpp::traits::input_parameter< int >::type trainL(trainLSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type bscales(bscalesSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type iscales(iscalesSEXP);
    rcpp_result_gen = Rcpp::wrap(RNNParamFit(input, W, RSInitial, washoutL, trainL, bscales, iscales));
    return rcpp_result_gen;
END_RCPP
}
// angleCalc
arma::mat angleCalc(arma::cube RS, arma::cube CRS);
RcppExport SEXP _conceptorCP_angleCalc(SEXP RSSEXP, SEXP CRSSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type RS(RSSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type CRS(CRSSEXP);
    rcpp_result_gen = Rcpp::wrap(angleCalc(RS, CRS));
    return rcpp_result_gen;
END_RCPP
}
// KSstatCalc
arma::vec KSstatCalc(arma::vec angles, float kappa);
RcppExport SEXP _conceptorCP_KSstatCalc(SEXP anglesSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type angles(anglesSEXP);
    Rcpp::traits::input_parameter< float >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(KSstatCalc(angles, kappa));
    return rcpp_result_gen;
END_RCPP
}
// CRNNBootstrap
arma::vec CRNNBootstrap(arma::cube bootinput, arma::cube W, arma::cube C, float RSInitial, int washoutL, int trainL, float bscale, float iscale, float kappa);
RcppExport SEXP _conceptorCP_CRNNBootstrap(SEXP bootinputSEXP, SEXP WSEXP, SEXP CSEXP, SEXP RSInitialSEXP, SEXP washoutLSEXP, SEXP trainLSEXP, SEXP bscaleSEXP, SEXP iscaleSEXP, SEXP kappaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::cube >::type bootinput(bootinputSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type W(WSEXP);
    Rcpp::traits::input_parameter< arma::cube >::type C(CSEXP);
    Rcpp::traits::input_parameter< float >::type RSInitial(RSInitialSEXP);
    Rcpp::traits::input_parameter< int >::type washoutL(washoutLSEXP);
    Rcpp::traits::input_parameter< int >::type trainL(trainLSEXP);
    Rcpp::traits::input_parameter< float >::type bscale(bscaleSEXP);
    Rcpp::traits::input_parameter< float >::type iscale(iscaleSEXP);
    Rcpp::traits::input_parameter< float >::type kappa(kappaSEXP);
    rcpp_result_gen = Rcpp::wrap(CRNNBootstrap(bootinput, W, C, RSInitial, washoutL, trainL, bscale, iscale, kappa));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_hello_world
arma::mat rcpparma_hello_world();
RcppExport SEXP _conceptorCP_rcpparma_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpparma_hello_world());
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_outerproduct
arma::mat rcpparma_outerproduct(const arma::colvec& x);
RcppExport SEXP _conceptorCP_rcpparma_outerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_outerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_innerproduct
double rcpparma_innerproduct(const arma::colvec& x);
RcppExport SEXP _conceptorCP_rcpparma_innerproduct(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_innerproduct(x));
    return rcpp_result_gen;
END_RCPP
}
// rcpparma_bothproducts
Rcpp::List rcpparma_bothproducts(const arma::colvec& x);
RcppExport SEXP _conceptorCP_rcpparma_bothproducts(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rcpparma_bothproducts(x));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_conceptorCP_NRMSE", (DL_FUNC) &_conceptorCP_NRMSE, 2},
    {"_conceptorCP_runRNN", (DL_FUNC) &_conceptorCP_runRNN, 5},
    {"_conceptorCP_conceptorCalc", (DL_FUNC) &_conceptorCP_conceptorCalc, 2},
    {"_conceptorCP_runCRNN", (DL_FUNC) &_conceptorCP_runCRNN, 7},
    {"_conceptorCP_WoutCalc", (DL_FUNC) &_conceptorCP_WoutCalc, 3},
    {"_conceptorCP_outputCalc", (DL_FUNC) &_conceptorCP_outputCalc, 2},
    {"_conceptorCP_RNNParamFit", (DL_FUNC) &_conceptorCP_RNNParamFit, 7},
    {"_conceptorCP_angleCalc", (DL_FUNC) &_conceptorCP_angleCalc, 2},
    {"_conceptorCP_KSstatCalc", (DL_FUNC) &_conceptorCP_KSstatCalc, 2},
    {"_conceptorCP_CRNNBootstrap", (DL_FUNC) &_conceptorCP_CRNNBootstrap, 9},
    {"_conceptorCP_rcpparma_hello_world", (DL_FUNC) &_conceptorCP_rcpparma_hello_world, 0},
    {"_conceptorCP_rcpparma_outerproduct", (DL_FUNC) &_conceptorCP_rcpparma_outerproduct, 1},
    {"_conceptorCP_rcpparma_innerproduct", (DL_FUNC) &_conceptorCP_rcpparma_innerproduct, 1},
    {"_conceptorCP_rcpparma_bothproducts", (DL_FUNC) &_conceptorCP_rcpparma_bothproducts, 1},
    {NULL, NULL, 0}
};

RcppExport void R_init_conceptorCP(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
