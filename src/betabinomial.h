#ifndef BETABINOMIAL_h
#define BETABINOMIAL_h

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;


arma::mat t1veccalc(arma::vec mu, arma::vec var);
arma::mat t2veccalc(arma::vec mu, arma::vec var);
arma::vec dbbinomV(arma::vec xm1, arma::vec xm2, const double t1, const double t2, bool log);
arma::vec pbbinomV(arma::vec ttxm, arma::vec xm1, const double t1, const double t2);

#endif
