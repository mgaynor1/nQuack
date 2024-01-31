// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include "emshared.h"
#include "ushared.h"
#include <math.h>


// title E-Step for Expected Maximization - Normal Distribution
// description This function is used in expected maximization. Here we complete
//  the E-Step and calculate the log-likelihood. Modifications include a correction for
//  the truncated distribution.
Rcpp::List estepN(const Rcpp::List parmlist, const arma::vec xi, std::string type,
                  const arma::vec trunc){

  const arma::vec avec = parmlist["avec"];
  const arma::vec mvec = parmlist["mvec"];
  const arma::vec svec = parmlist["svec"];

  const int m = mvec.size();

  //int pl = (n == m ? 1 : 2);

  const int xl = xi.size();
  arma::mat zprob( xl, m );

  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0){
    for(int i  = 0; i < m; i++) {
      arma::vec ppv = dnormV(xi, mvec(i), svec(i), 0);
      zprob.col(i) = (avec(i) * ppv);
    }
  } else{
  for(int i  = 0; i < m; i++) {
    double F2F1 = R::pnorm(trunc(1), mvec(i), svec(i), 1, 0) - R::pnorm(trunc(0), mvec(i), svec(i), 1, 0);
    arma::vec pp = dnormV(xi,  mvec(i), svec(i), 0);
    pp = (pp/F2F1);
    zprob.col(i) = (avec(i) * pp);
  }
}


  arma::vec denom = arma::sum(zprob, 1);

  arma::mat zprobb( xl, m );
  arma::vec sjvec(m );
  for(int i = 0; i < m ; i++){
    zprobb.col(i) = ( zprob.col(i) /denom ) ;
    sjvec(i) = sum(zprobb.col(i));
  }
  double ll = sum(log(denom));

  Rcpp::List P = Rcpp::List::create( Rcpp::Named("avec") = avec, Rcpp::_["mvec"] = mvec, Rcpp::_["svec"] = svec);
  return( Rcpp::List::create( Rcpp::_["zprob"] = zprobb,  Rcpp::_["parm.list"] = P, Rcpp::_["xi"] = xi, Rcpp::_["denom"] = denom, Rcpp::_["trunc"] = trunc, Rcpp::_["type"]=type, Rcpp::_["sjvec"]=sjvec, Rcpp::_["LL"]=ll));

}





// M-Step with Numerical Optimization for Expected Maximization - Normal Distribution
//
// description This function is used in expected maximization to maximize the parameter values.
// param eout List with output from the estep
Rcpp::List mstepN(Rcpp::List eout){

  const Rcpp::List parmlist = eout["parm.list"];
  const arma::vec avec = parmlist["avec"];
  const arma::vec mvec = parmlist["mvec"];
  const arma::vec svec = parmlist["svec"];
  const arma::mat Zprobsmat = eout["zprob"];
  const arma::vec xi = eout["xi"];
  //const arma::vec denom = eout["denom"];
  std::string type = eout["type"];
  //const arma::vec trunc = eout["trunc"];
  const arma::vec sjvec = eout["sjvec"];

  //int nmixt = avec.size();
  int m = mvec.size();
  int n = xi.n_rows;

  arma::vec ahats(m);
  arma::vec mhats(m);
  arma::vec shats(m);

  if(type == "free"){
    for(int i = 0; i < m; i++){
      arma::vec gamma = Zprobsmat.col(i);
      double sj = sjvec(i);
      arma::vec sp = arma::square(xi - mvec(i));
      ahats(i) = (sj/n);
      mhats(i) = ((1/sj)*(sum(gamma % xi)));
      shats(i) = sqrt((1/sj)*(sum(gamma % sp)));

    }
  } else if(type == "fixed"){
      for(int i = 0; i < m; i++){
        arma::vec gamma = Zprobsmat.col(i);
        double sj = sjvec(i);
        arma::vec sp = square(xi - mvec(i));
        ahats(i) = (sj/n);
        mhats(i) = mvec(i);
        shats(i) = svec(i);
    }
  } else if(type == "fixed_2"){
      for(int i = 0; i < m; i++){
        arma::vec gamma = Zprobsmat.col(i);
        double sj = sjvec(i);
        arma::vec sp = square(xi - mvec(i));
        ahats(i) = (sj/n);
        mhats(i) = mvec(i);
        shats(i) = sqrt((1/sj)*(sum(gamma % sp)));
      }
  } else if(type == "fixed_3"){
      for(int i = 0; i < m; i++){
        arma::vec gamma = Zprobsmat.col(i);
        double sj = sjvec(i);
        arma::vec sp = square(xi - mvec(i));
        ahats(i) = avec(i);
        mhats(i) = mvec(i);
        shats(i) = sqrt((1/sj)*(sum(gamma % sp)));
      }
  }

    // Logic check
    if(ahats.has_nan() == true || all(ahats) == false){ //if all are non.zero = false
      //Rcpp::Rcout << ahats << std::endl;
      ahats = avec;
    }
     if(shats.has_nan() == true || all(shats) == false){
       //Rcpp::Rcout << shats << std::endl;
       shats = svec;
     }
     if(mhats.has_nan() == true || all(mhats) == false){
      //Rcpp::Rcout << mhats << std::endl;
       mhats = mvec;

     }
     // Add the Uniform
     //double ah = sum(ahats);
    // ahats.resize(nmixt);
     //ahats(nmixt-1) = (1 - ah);

  Rcpp::List res = Rcpp::List::create(Rcpp::_["avec"] = ahats, Rcpp::_["mvec"] = mhats, Rcpp::_["svec"] = shats);
  return(res);
}




//' @title Expected maximization - Normal Distribution
//'
//' @description This function calculates the log-likelihood using
//'  the expected maximization algorithm with the Normal Distribution.
//'  This code follows nQuire and does not use an augmented likelihood.
//'
//' @param parmlist A list containing initial alpha, mean, and variance values.
//' @param xi List of observations, in this case allele frequencies.
//' @param niter Max number of iterates.
//' @param epsilon Epsilon value for convergence tolerance. When the absolute delta log-likelihood is
//'    below this value, convergence is reached.
//' @param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
//' @param type String indicating model type. Options: "free" (estimated parameter(s): alpha, mean, and variance), "fixed" (estimated parameter(s): alpha),
//' "fixed_2" (estimated parameter(s): alpha and variance), or "fixed_3" (estimated parameter(s): variance).
//'
//' @returns List of elements including the log-likelihood, the number of iterates,
//' and the optimized parameter values.
//'
// [[Rcpp::export]]
Rcpp::List emstepN(Rcpp::List parmlist, const arma::vec xi, int niter, double epsilon, arma::vec trunc,  std::string type = "free"){
   // Set up progress

   Rcpp::List mint = parmlist;
   // E-step
   Rcpp::List eint = estepN(mint, xi, type, trunc);
   double lnlikeint = llcalcfinalN(eint);

   int count = 0;
   // Iterating E and M steps until convergence is below epsilon
   for(int j = 0; j < (niter + 1); j++){
     Rcpp::checkUserInterrupt();
     //p.increment();
     count = count + 1;

     // M-step
     Rcpp::List moutip = mstepN(eint);


     // E-step
     Rcpp::List eoutip = estepN(moutip, xi, type, trunc);
     double lnlikeip = llcalcfinalN(eoutip);
     double deltalogL = (lnlikeint - lnlikeip);
     deltalogL = std::abs(deltalogL);

     // Set next iter
     eint = eoutip;
     lnlikeint = lnlikeip;

     if(deltalogL<epsilon){
       break;
     }
   }
   double llreturn = llcalcfinalN(eint);
   return  Rcpp::List::create(Rcpp::_["loglikelihood"] = llreturn, Rcpp::_["LL"]=eint["LL"], Rcpp::_["parm.list"] = eint["parm.list"], Rcpp::_["niter.done"] = count);
 }


