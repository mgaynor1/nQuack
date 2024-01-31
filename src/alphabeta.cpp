// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;


//' Calculate Alpha and Beta from Mean and Variance
//'
//' @param mu Mean.
//' @param var Variance.
//'
//' @return Numeric vector of alpha and beta.
// [[Rcpp::export]]
arma::vec alphabetacalc(const double mu, const double var){
 arma::vec ab(2);
 const double comm = ( ( ( mu * (1.0 - mu) ) / var ) - 1.0);
 ab[0] = (mu * comm);
 ab[1] = ( ( 1.0 - mu ) * comm);
 return (ab);
}


//' Vector-based - Calculate Alpha and Beta from Mean and Variance
//'
//' @param mu Vector of mean.
//' @param var Vector of variance.
//'
//' @return Numeric matrix of alpha and beta.
// [[Rcpp::export]]
arma::mat alphabetacalcvec(arma::vec mu, arma::vec var){
  const int ml = mu.size();
  arma::mat ab(ml, 2);
  vec y(mu.size(), fill::ones);
  arma::vec comm = ( ( ( mu % (y - mu) ) / var ) - y);
  ab.col(0) = (mu % comm);
  ab.col(1) = ( ( y - mu ) % comm);
  return (ab);
}


//' Calculate Alpha and Beta from Mean, Tau, and Error rate.
//'
//' @param mu Mean.
//' @param tau Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less overdispersion and 1 indicates high overdispersion.  Here tau must be greater than 0.
//' @param error Sequencing error rate.
//'
//' @return Numeric vector of alpha and beta.
// [[Rcpp::export]]
 arma::vec alphabetacalctau(const double mu, const double tau, const double error){
   arma::vec ab(2);
   const double mm = ((mu*(1.0-error)) + ((1.0-mu)*error));
   ab[0] = (((1.0-tau)*mm)/tau);
   ab[1] = (((1.0-mm)*(1.0-tau))/tau);
   return (ab);
 }


//' Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error rate.
//'
//' @param mu Vector of mean.
//' @param tau Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less overdispersion and 1 indicates high overdispersion.  Here tau must be greater than 0.
//' @param error Sequencing error rate. Ranges from 0 to 1.
//'
//' @return Numeric matrix of alpha and beta.
// [[Rcpp::export]]
arma::mat alphabetacalctauvec(arma::vec mu, const double tau, const double error){
  const int ml = mu.size();
  arma::mat ab(ml, 2);
  arma::vec y(mu.size(), fill::ones);
  arma::vec t(mu.size(), fill::value(tau));
  arma::vec erv(mu.size(), fill::value(error));
  arma::vec mm = ((mu%(y-erv)) + ((y-mu)%erv));
  ab.col(0) = (((y-t)%mm)/t);
  ab.col(1) = (((y-mm)%(y-t))/t);
  return (ab);
}





