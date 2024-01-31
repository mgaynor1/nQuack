#include <Rcpp.h>
using namespace Rcpp;

//' @title Data Preparation - Use nQuire's Data
//'
//' @description This function reduce a three column data frame to
//' two columns by randomly sampling allele A or B for every site. This is used in our function `process_nquire()`
//'
//' @param xm A matrix with three columns: Total Coverage, Counts for Allele A, and Counts for Allele B.
//'
//' @returns Numeric Matrix with total coverage and coverage for a randomly sampled allele.
// [[Rcpp::export]]
NumericMatrix nQuire_reformat(NumericMatrix xm){
  int xl = xm.nrow();
  NumericMatrix xmo(xl, 2);
  xmo(_, 0) = xm(_, 0);
  NumericVector prob = {0.5, 0.5};
  for(int i = 0; i < xl; i++){
    NumericVector y(2);
    y[0] = xm(i, 1);
    y[1] = xm(i, 2);
    double samp = as<double>(sample(y, 1, false, prob));
    xmo(i, 1) = samp;
  }
  return(xmo);
}

