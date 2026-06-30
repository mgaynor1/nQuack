#ifndef EXTRADISTRBBD_H
#define EXTRADISTRBBD_H

#include <Rcpp.h>

/*
 * This file is derived from the extraDistr R package.
 *
 * Original source:
 * extraDistr:
 * https://github.com/twolodzko/extraDistr/blob/master/src/beta-binomial-distribution.cpp
 * https://github.com/twolodzko/extraDistr/blob/master/src/shared.h
 *
 * License: GPL-2
 */

using Rcpp::NumericVector;


NumericVector cpp_dbbinom( const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& log_prob = false);
NumericVector cpp_pbbinom( const NumericVector& x, const NumericVector& size, const NumericVector& alpha, const NumericVector& beta, const bool& lower_tail = true, const bool& log_prob = false);

#endif
