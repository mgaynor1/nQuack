// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include "emshared.h"

// dnormV
//
// param xi Vector of doubles.
// param a Mean parameter.
// param b Standard deviation parameter
// param log Interger indicating to log or not to log.
//
// return Numeric vector from dbeta.
arma::vec dnormV(arma::vec xi, double a, double b, int log){
  int xil = xi.size();
  arma::vec out(xil);
  for(int i = 0; i < xil; i++){
    out(i) = R::dnorm(xi(i), a, b, log);
    //out(i) = (1/(b * sqrt(2*M_PI))) * exp( -( pow(xi(i) - a, 2) / (2 * pow(b, 2)) ));
  }
  return out;
}

// Log-Likelihood for BIC - Normal
double llcalcfinalN(Rcpp::List eout){
  Rcpp::List parmlist = eout["parm.list"];
  arma::vec avec = parmlist["avec"];
  arma::vec mvec = parmlist["mvec"];
  arma::vec svec = parmlist["svec"];
  arma::vec xi = eout["xi"];
  arma::vec trunc = eout["trunc"];

  const int m = mvec.size();
  int xl = xi.size();

  arma::mat mixit(xl, m);
  // double ua =  1.0 - sum(avec);
  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
    for(int i = 0; i < m; i++){
      mixit.col(i) =  avec(i)*(dnormV(xi, mvec(i), svec(i), 0));
    }
  }else{
    for(int i  = 0; i < m; i++) {
      double F2F1 = R::pnorm(trunc(1), mvec(i), svec(i), 1, 0) - R::pnorm(trunc(0), mvec(i), svec(i), 1, 0);
      arma::vec pp = (avec(i)*dnormV(xi, mvec(i), svec(i), 0));
      mixit.col(i) =  (pp/F2F1);
    }
  }

  // Uniform distribution
  //mixit.col(m) = (avec(m)*prob_unif_vec(xi));

  arma::vec sumrows = arma::sum(mixit, 1);
  if(all(sumrows) == false){ //if all are non.zero = false
    for(int i = 0; i < xl; i++){
      if(sumrows(i) == 0){
        sumrows(i) = DBL_MIN;
      }
    }
  }
  double sumit = sum(log(sumrows));
  return sumit;
}

// Log-Likelihood for BIC - Normal + Uniform
double llcalcfinalNU(Rcpp::List eout){
  Rcpp::List parmlist = eout["parm.list"];
  arma::vec avec = parmlist["avec"];
  arma::vec mvec = parmlist["mvec"];
  arma::vec svec = parmlist["svec"];
  arma::vec xi = eout["xi"];
  arma::vec trunc = eout["trunc"];

  const int m = mvec.size();
  int xl = xi.size();

  arma::mat mixit(xl, (m+1));
  // double ua =  1.0 - sum(avec);
  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
    for(int i = 0; i < m; i++){
      mixit.col(i) =  avec(i)*(dnormV(xi, mvec(i), svec(i), 0));
    }
  }else{
    for(int i  = 0; i < m; i++) {
      double F2F1 = R::pnorm(trunc(1), mvec(i), svec(i), 1, 0) - R::pnorm(trunc(0), mvec(i), svec(i), 1, 0);
      arma::vec pp = (avec(i)*dnormV(xi, mvec(i), svec(i), 0));
      mixit.col(i) =  (pp/F2F1);
    }
  }

  // Uniform distribution
  mixit.col(m) = (avec(m)*prob_unif_vec(xi));

  arma::vec sumrows = arma::sum(mixit, 1);
  if(all(sumrows) == false){ //if all are non.zero = false
    for(int i = 0; i < xl; i++){
      if(sumrows(i) == 0){
        sumrows(i) = DBL_MIN;
      }
    }
  }
  double sumit = sum(log(sumrows));
  return sumit;
}
