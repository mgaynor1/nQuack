// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>

using namespace arma;

//' Calculate Alpha and Beta from Mean and Variance
//'
//' @param xm Matrix with total coverage and coverage at a randomly sampled allele.
//' @param n Length of matrix.
//'
//' @returns Randomly sampled matrix.
//' @examples
//' outdf <- resample_xm(as.matrix(xm), n = 10)
// [[Rcpp::export]]
arma::mat resample_xm(arma::mat xm, int n){
  arma::uvec index = arma::linspace<arma::uvec>(0, n-1, n);
  arma::uvec locs = Rcpp::RcppArmadillo::sample(index,n,true);
  arma::mat  resample = xm.rows(locs);
  return resample;
}



