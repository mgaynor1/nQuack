// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include "extradistrbbd.h"

// t1vec calc
arma::mat t1veccalc(arma::vec mu, arma::vec var){
  vec y(mu.size(), fill::ones);
  arma::vec t1vec = (( ( ( mu % (y - mu) ) / var ) - y) % mu);
  return (t1vec);
}

//t2vec
arma::mat t2veccalc(arma::vec mu, arma::vec var){
  vec y(mu.size(), fill::ones);
  arma::vec t2vec = (( ( ( mu % (y - mu) ) / var ) - y) % (y - mu));
  return (t2vec);
}

// dbbinomV
// param xm1 vector of x.
// param xm2 vector of size
// param t1 alpha
// param t2 beta
// param log logical; if TRUE, dbbinom is given as log.
arma::vec dbbinomV(arma::vec xm1, arma::vec xm2, const double t1, const double t2, bool log){


  Rcpp::NumericVector rxm1 = Rcpp::NumericVector(Rcpp::wrap(xm1));
  Rcpp::NumericVector rxm2 = Rcpp::NumericVector(Rcpp::wrap(xm2));
  Rcpp::NumericVector t1v =  Rcpp::NumericVector::create(t1);
  Rcpp::NumericVector t2v =  Rcpp::NumericVector::create(t2);


  Rcpp::NumericVector samp =  cpp_dbbinom(rxm1, rxm2, t1v, t2v, log);

  arma::vec out = Rcpp::as<arma::vec>(Rcpp::wrap(samp));
  return(out);
}

// pbbinomV
// param ttxm vector of x.
// param xm2 vector of size
// param t1 alpha
// param t2 beta
arma::vec pbbinomV(arma::vec ttxm, arma::vec xm1, const double t1, const double t2){

  Rcpp::NumericVector rttx = Rcpp::NumericVector(Rcpp::wrap(ttxm));
  Rcpp::NumericVector rxm1 = Rcpp::NumericVector(Rcpp::wrap(xm1));
  Rcpp::NumericVector t1v = Rcpp::NumericVector::create(t1);
  Rcpp::NumericVector t2v = Rcpp::NumericVector::create(t2);


  Rcpp::NumericVector samp = cpp_pbbinom(rttx, rxm1, t1v, t2v);

  arma::vec out = Rcpp::as<arma::vec>(Rcpp::wrap(samp));
  return(out);
}
