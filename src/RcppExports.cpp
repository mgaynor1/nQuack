// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include "../inst/include/nQuack.h"
#include <RcppArmadillo.h>
#include <Rcpp.h>
#include <string>
#include <set>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// estepB3
Rcpp::List estepB3(const Rcpp::List parmlist, const arma::vec xi, const arma::vec trunc);
RcppExport SEXP _nQuack_estepB3(SEXP parmlistSEXP, SEXP xiSEXP, SEXP truncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type trunc(truncSEXP);
    rcpp_result_gen = Rcpp::wrap(estepB3(parmlist, xi, trunc));
    return rcpp_result_gen;
END_RCPP
}
// emstepB3
Rcpp::List emstepB3(Rcpp::List parmlist, arma::vec xi, int niter, double epsilon, arma::vec trunc);
RcppExport SEXP _nQuack_emstepB3(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepB3(parmlist, xi, niter, epsilon, trunc));
    return rcpp_result_gen;
END_RCPP
}
// emstepN
Rcpp::List emstepN(Rcpp::List parmlist, const arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepN(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepN(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepNA
Rcpp::List emstepNA(Rcpp::List parmlist, const arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepNA(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepNA(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepNU
Rcpp::List emstepNU(Rcpp::List parmlist, const arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepNU(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepNU(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepNUA
Rcpp::List emstepNUA(Rcpp::List parmlist, const arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepNUA(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepNUA(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepB
Rcpp::List emstepB(Rcpp::List parmlist, arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepB(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepB(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepBB
Rcpp::List emstepBB(Rcpp::List parmlist, arma::mat xm, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepBB(SEXP parmlistSEXP, SEXP xmSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xm(xmSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepBB(parmlist, xm, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepBBU
Rcpp::List emstepBBU(Rcpp::List parmlist, arma::mat xm, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepBBU(SEXP parmlistSEXP, SEXP xmSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type xm(xmSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepBBU(parmlist, xm, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// emstepBU
Rcpp::List emstepBU(Rcpp::List parmlist, arma::vec xi, int niter, double epsilon, arma::vec trunc, std::string type);
RcppExport SEXP _nQuack_emstepBU(SEXP parmlistSEXP, SEXP xiSEXP, SEXP niterSEXP, SEXP epsilonSEXP, SEXP truncSEXP, SEXP typeSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::List >::type parmlist(parmlistSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type xi(xiSEXP);
    Rcpp::traits::input_parameter< int >::type niter(niterSEXP);
    Rcpp::traits::input_parameter< double >::type epsilon(epsilonSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< std::string >::type type(typeSEXP);
    rcpp_result_gen = Rcpp::wrap(emstepBU(parmlist, xi, niter, epsilon, trunc, type));
    return rcpp_result_gen;
END_RCPP
}
// alphabetacalc
arma::vec alphabetacalc(const double mu, const double var);
RcppExport SEXP _nQuack_alphabetacalc(SEXP muSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(alphabetacalc(mu, var));
    return rcpp_result_gen;
END_RCPP
}
// alphabetacalcvec
arma::mat alphabetacalcvec(arma::vec mu, arma::vec var);
RcppExport SEXP _nQuack_alphabetacalcvec(SEXP muSEXP, SEXP varSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type var(varSEXP);
    rcpp_result_gen = Rcpp::wrap(alphabetacalcvec(mu, var));
    return rcpp_result_gen;
END_RCPP
}
// alphabetacalctau
arma::vec alphabetacalctau(const double mu, const double tau, const double error);
RcppExport SEXP _nQuack_alphabetacalctau(SEXP muSEXP, SEXP tauSEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const double >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(alphabetacalctau(mu, tau, error));
    return rcpp_result_gen;
END_RCPP
}
// alphabetacalctauvec
arma::mat alphabetacalctauvec(arma::vec mu, const double tau, const double error);
RcppExport SEXP _nQuack_alphabetacalctauvec(SEXP muSEXP, SEXP tauSEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type mu(muSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    Rcpp::traits::input_parameter< const double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(alphabetacalctauvec(mu, tau, error));
    return rcpp_result_gen;
END_RCPP
}
// prepare_data
void prepare_data(std::string name, std::string inpath, std::string outpath, std::string threads);
RcppExport SEXP _nQuack_prepare_data(SEXP nameSEXP, SEXP inpathSEXP, SEXP outpathSEXP, SEXP threadsSEXP) {
BEGIN_RCPP
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< std::string >::type name(nameSEXP);
    Rcpp::traits::input_parameter< std::string >::type inpath(inpathSEXP);
    Rcpp::traits::input_parameter< std::string >::type outpath(outpathSEXP);
    Rcpp::traits::input_parameter< std::string >::type threads(threadsSEXP);
    prepare_data(name, inpath, outpath, threads);
    return R_NilValue;
END_RCPP
}
// resample_xm
arma::mat resample_xm(arma::mat xm, int n);
RcppExport SEXP _nQuack_resample_xm(SEXP xmSEXP, SEXP nSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type xm(xmSEXP);
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    rcpp_result_gen = Rcpp::wrap(resample_xm(xm, n));
    return rcpp_result_gen;
END_RCPP
}
// nQuire_reformat
NumericMatrix nQuire_reformat(NumericMatrix xm);
RcppExport SEXP _nQuack_nQuire_reformat(SEXP xmSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type xm(xmSEXP);
    rcpp_result_gen = Rcpp::wrap(nQuire_reformat(xm));
    return rcpp_result_gen;
END_RCPP
}
// process_rcpp
NumericMatrix process_rcpp(NumericMatrix x, int mindepth, double maxprob, NumericVector trunc, double error);
RcppExport SEXP _nQuack_process_rcpp(SEXP xSEXP, SEXP mindepthSEXP, SEXP maxprobSEXP, SEXP truncSEXP, SEXP errorSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    Rcpp::traits::input_parameter< int >::type mindepth(mindepthSEXP);
    Rcpp::traits::input_parameter< double >::type maxprob(maxprobSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type trunc(truncSEXP);
    Rcpp::traits::input_parameter< double >::type error(errorSEXP);
    rcpp_result_gen = Rcpp::wrap(process_rcpp(x, mindepth, maxprob, trunc, error));
    return rcpp_result_gen;
END_RCPP
}

// validate (ensure exported C++ functions exist before calling them)
static int _nQuack_RcppExport_validate(const char* sig) { 
    static std::set<std::string> signatures;
    if (signatures.empty()) {
    }
    return signatures.find(sig) != signatures.end();
}

// registerCCallable (register entry points for exported C++ functions)
RcppExport SEXP _nQuack_RcppExport_registerCCallable() { 
    R_RegisterCCallable("nQuack", "_nQuack_RcppExport_validate", (DL_FUNC)_nQuack_RcppExport_validate);
    return R_NilValue;
}

static const R_CallMethodDef CallEntries[] = {
    {"_nQuack_estepB3", (DL_FUNC) &_nQuack_estepB3, 3},
    {"_nQuack_emstepB3", (DL_FUNC) &_nQuack_emstepB3, 5},
    {"_nQuack_emstepN", (DL_FUNC) &_nQuack_emstepN, 6},
    {"_nQuack_emstepNA", (DL_FUNC) &_nQuack_emstepNA, 6},
    {"_nQuack_emstepNU", (DL_FUNC) &_nQuack_emstepNU, 6},
    {"_nQuack_emstepNUA", (DL_FUNC) &_nQuack_emstepNUA, 6},
    {"_nQuack_emstepB", (DL_FUNC) &_nQuack_emstepB, 6},
    {"_nQuack_emstepBB", (DL_FUNC) &_nQuack_emstepBB, 6},
    {"_nQuack_emstepBBU", (DL_FUNC) &_nQuack_emstepBBU, 6},
    {"_nQuack_emstepBU", (DL_FUNC) &_nQuack_emstepBU, 6},
    {"_nQuack_alphabetacalc", (DL_FUNC) &_nQuack_alphabetacalc, 2},
    {"_nQuack_alphabetacalcvec", (DL_FUNC) &_nQuack_alphabetacalcvec, 2},
    {"_nQuack_alphabetacalctau", (DL_FUNC) &_nQuack_alphabetacalctau, 3},
    {"_nQuack_alphabetacalctauvec", (DL_FUNC) &_nQuack_alphabetacalctauvec, 3},
    {"_nQuack_prepare_data", (DL_FUNC) &_nQuack_prepare_data, 4},
    {"_nQuack_resample_xm", (DL_FUNC) &_nQuack_resample_xm, 2},
    {"_nQuack_nQuire_reformat", (DL_FUNC) &_nQuack_nQuire_reformat, 1},
    {"_nQuack_process_rcpp", (DL_FUNC) &_nQuack_process_rcpp, 5},
    {"_nQuack_RcppExport_registerCCallable", (DL_FUNC) &_nQuack_RcppExport_registerCCallable, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_nQuack(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
