// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// Helper
// Numeric vector exponential sum
//
// param v Numeric vector.
//
// return Sum of exponential of each value in a numeric vector.
double exp_sum(const arma::vec v){
  arma::vec vv = exp(v);
  return(std::accumulate(vv.begin(), vv.end(), 0.0));
}


// LINK FUNCTIONS
// mlogitc
// description To force the final alpha estimates to be between 0 and 1, and add up to 1,
//  we use logit transformations.
// param avec List of alpha of nmixt length.
// param nmixt Length of alphas.
//
// return Logit-transform of the first (nmixt-1) alpha values
arma::vec mlogitcU(const arma::vec avec, const int nmixt){
  arma::vec pis = avec(span(0, (nmixt-2)));
  double pim = avec((nmixt-1));
  return log(pis/pim);
}

double mlogitcD(const arma::vec avec, const int nmixt){
  double pis = avec(0);
  double pim = avec(1);
  return log(pis/pim);
}


// invmlogitc
//
// param betas List of betas to transform.
//
// return Numeric vector of alpha and beta.
arma::vec invmlogitcA(const arma::vec betas){
 double denom = (1.0 + exp_sum(betas));
 arma::vec  nom = exp(betas);
 arma::vec  out = nom/denom;
 return(out);
}
arma::vec invmlogitcU(const arma::vec betas){
  int b = betas.size();
  arma::vec out(b+1);
  double denom = (1.0 + exp_sum(betas));
  arma::vec  nom = exp(betas);
  out(span(0, b-1)) = (nom/denom);
  out(b) = (1.0 - sum(out(span(0, b-1))));
  return(out);
}

arma::vec invmlogitcD(const double betas){
  arma::vec out(2);
  double denom = (1.0 + exp(betas));
  out(0) = (exp(betas)/denom);
  out(1) = (1.0 - out(0));
  return(out);
}

// prob_unif_nq - from nQuire
//
// param minfrac Minimum allele frequency.
// param x Allele frequency at a site.
// return Uniform distribution.
double prob_unif_nq(double minfrac, double x){
  double unif = 0;
  double upper = 1.0 - minfrac;
  if(x >= minfrac && x <= upper){
    unif = 1 / (upper - minfrac);
  }
  return(unif);
}

// prob_unif_vec - vectorized
//
// param xi Vector of allele frequency.
//
// return Sampled uniform distribution.
arma::vec prob_unif_vec(arma::vec xi, int logit = 0){
  double minfrac = min(xi);
  double upper = max(xi);
  double dif = (upper - minfrac);

  arma::vec unif(xi.size());

  for(int i = 0; i < xi.size(); i++){
    if(xi(i) >= minfrac && xi(i) <= upper){
      if(logit == 0){
      unif(i) = (1 / dif);
      } else if (logit == 1){
        unif(i) = log(1.0/dif);
      }
    } else{
      unif(i) = 0;
    }
  }
  return(unif);
}



// prob_unif_vec - for Beta-Binomial
//
// param xi Vector of allele frequency.
//
// return Sampled uniform distribution.
arma::vec prob_unif_vec_bb(arma::vec xm, int logit = 0){
  arma::vec unif(xm.size());

  for(int i = 0; i < xm.size(); i++){
     if(logit == 0){
        unif(i) = (1 / xm(i));
      } else if (logit == 1){
        unif(i) = log(1.0/xm(i));
      }
  }
  return(unif);
}

// Alpha and Beta calculation
//
// param mu Mean.
// param var Variance.
//
// return Numeric vector of alpha and beta.
arma::vec alphabetacalcA(const double mu, const double var){
  arma::vec ab(2);
  const double comm = ( ( ( mu * (1.0 - mu) ) / var ) - 1.0);
  ab[0] = (mu * comm);
  ab[1] = ( ( 1.0 - mu ) * comm);
  return (ab);
}

// mutransform
//
// param mvec List of means.
//
// return Numeric vector transformed
arma::vec mutransformA(const arma::vec mvec){
  arma::vec mveout = ((exp(mvec * -1.0)) + 1.0);
  return 1/mveout;
}

// mvecin
//
// param mvec List of means.
//
// return Numeric vector transformed
arma::vec mvecinA(const arma::vec mvec){
  arma::vec l = log(mvec);
  arma::vec r = log(1.0 - mvec);
  return (l-r);
}

// stransform
//
// param svec List of variance.
//
// return Numeric vector transformed.
arma::vec stransformA(const arma::vec svec){
  return(exp(svec));
}

// svecin
//
// param svec List of variance.
//
// return Numeric vector transformed
arma::vec svecinA(const arma::vec svec){
  return(log(svec));
}



// dbetaV
//
// param xi Vector of doubles.
// param a Alpha parameter.
// param b Beta parameter
// param log Interger indicating to log or not to log.
//
// return Numeric vector from dbeta.
arma::vec dbetaV(arma::vec xi, double a, double b, int log){
  int xil = xi.size();
  arma::vec out(xil);
  for(int i = 0; i < xil; i++){
    out(i) = R::dbeta(xi(i), a, b, log);
  }
  return out;
}


// Flag sample
arma::vec sampFT(std::string type, int m){
  arma::vec iglpl(3);
  if(type == "free"){
    iglpl(0) = 0;
    if(m == 1){
      iglpl(1) = 2;
      iglpl(2) = 1;
    } else{
      iglpl(1) = ((3*m)-1);
      iglpl(2) = 1;
    }
  } else if(type == "fixed"){ //alpha
    iglpl(0) = 1;
    iglpl(1) = (m-1);
    iglpl(2) = (2*m);
  } else if(type == "fixed_2"){ //both
    iglpl(0) = 2;
    if(m == 1){
      iglpl(1) = 1;
      iglpl(2) = 2;
    } else{
      iglpl(1) = ((2*m)-1);
      iglpl(2) = m;
    }
  } else if(type == "fixed_3"){
    iglpl(0) = 3;
    iglpl(1) = m;
    iglpl(2) = (2*m);
  } else{
    Rcpp::Rcerr << "Parameter is not allowed. Type may equal free, fixed, fixed_2, or fixed_3" << std::endl;
  }
  return iglpl;
}

// Flag type
// Distribution-Uniform-only
arma::vec sampFTBU(std::string type, int m){
  arma::vec iglpl(3);
  if(type == "free"){
    iglpl(0) = 0;
    iglpl(1) = (3*m);
    iglpl(2) = 1;
  } else if(type == "fixed"){ //alpha
    iglpl(0) = 1;
    iglpl(1) = m;
    iglpl(2) = (2*m);
  } else if(type == "fixed_2"){ //both
    iglpl(0) = 2;
    iglpl(1) = (2*m);
    iglpl(2) = m;
  } else if(type == "fixed_3"){
    iglpl(0) = 3;
    iglpl(1) = m;
    iglpl(2) = ((2*m)+1);
  } else{
    Rcpp::Rcerr << "Parameter is not allowed. Type may equal free, fixed, fixed_2, or fixed_3" << std::endl;
  }
  return iglpl;
}

