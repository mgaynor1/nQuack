#ifndef EXTRADISTRBBD_H
#define EXTRADISTRBBD_H

#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]


using Rcpp::NumericVector;


NumericVector cpp_dbbinom( const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& log_prob = false);
NumericVector cpp_pbbinom( const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail = true, const bool& log_prob = false);

#endif
