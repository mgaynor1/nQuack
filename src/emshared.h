#ifndef EMSHARED_H
#define EMSHARED_H


// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;


double exp_sum(const arma::vec v);
arma::vec mlogitcU(const arma::vec avec, const int nmixt);
double mlogitcD(const arma::vec avec, const int nmixt);
arma::vec invmlogitcA(const arma::vec betas);
arma::vec invmlogitcU(const arma::vec betas);
arma::vec invmlogitcD(const double betas);
double prob_unif_nq(double minfrac, double x);
arma::vec prob_unif_vec(arma::vec xi, int logit = 0);
arma::vec prob_unif_vec_bb(arma::vec xm, int logit = 0);
arma::vec alphabetacalcA(const double mu, const double var);
arma::vec mutransformA(const arma::vec mvec);
arma::vec mvecinA(const arma::vec mvec);
arma::vec stransformA(const arma::vec svec);
arma::vec svecinA(const arma::vec svec);
arma::vec dbetaV(arma::vec xi, double a, double b, int log);
arma::vec sampFT(std::string type, int m);
arma::vec sampFTBU(std::string type, int m);

#endif
