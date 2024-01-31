#ifndef USHARED_H
#define USHARED_H

// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

arma::vec dnormV(arma::vec xi, double a, double b, int log);
double llcalcfinalNU(Rcpp::List eout);
double llcalcfinalN(Rcpp::List eout);


#endif
