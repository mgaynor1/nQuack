// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;
// [[Rcpp::depends(RcppProgress)]]
#include <progress.hpp>
#include "emshared.h"



//' @title E-Step for Expected Maximization - Beta + Beta + Beta Distribution
//' @description This is used in the `Bclean()` function. Here we complete
//' the E-Step and calculate the log-likelihood. Modifications include a correction for
//' the truncated distribution.
//' @param parmlist A list containing initial alpha, mean, and variance.
//' @param xi List of observations, in this case allele frequencies.
//' @param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
//[[Rcpp::export]]
Rcpp::List estepB3(const Rcpp::List parmlist, const arma::vec xi,
                   const arma::vec trunc){

  const arma::vec avec = parmlist["avec"];
  const arma::vec t1vec = parmlist["t1vec"];
  const arma::vec t2vec = parmlist["t2vec"];

  const int m = avec.size();
  const int xl = xi.size();
  arma::mat zprob( xl, m );

  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0){
    for(int i  = 0; i < m; i++) {
      arma::vec ppv = dbetaV(xi, t1vec(i), t2vec(i), 0);
      zprob.col(i) = (avec(i) * ppv);
    }
  } else{
    for(int i  = 0; i < m; i++) {
      double F2F1 = (R::pbeta(trunc(1), t1vec(i), t2vec(i), 1, 0) - R::pbeta(trunc(0), t1vec(i), t2vec(i), 1, 0));
      arma::vec pp = dbetaV(xi, t1vec(i), t2vec(i),0);
      pp = (pp/F2F1);
      zprob.col(i) = (avec(i) * pp);
    }
  }

  arma::vec denom = arma::sum(zprob, 1);

  arma::mat zprobb( xl, m );
  for(int i = 0; i < m; i++){
    zprobb.col(i) = ( zprob.col(i) /denom ) ;
  }

  Rcpp::List P = Rcpp::List::create( Rcpp::Named("avec") = avec, Rcpp::_["t1vec"] = t1vec, Rcpp::_["t2vec"] = t2vec);
  return( Rcpp::List::create( Rcpp::_["zprob"] = zprobb,  Rcpp::_["parm.list"] = P, Rcpp::_["xi"] = xi, Rcpp::_["denom"] = denom, Rcpp::_["trunc"] = trunc));

}

// Augmented Likelihood Calculation - beta Distribution
double llcalcB3(const arma::vec avec, const arma::vec t1vec, const arma::vec t2vec,
                const int zcols, const arma::vec xi,  const arma::mat zprob,
                const arma::vec denom,  const arma::vec trunc) {

  double sumit = 0.0;
  arma::vec lavec = log(avec);

  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
    for(int i = 0; i < zcols; i++){
      arma::vec right = lavec(i) + dbetaV(xi, t1vec(i), t2vec(i), 1);
      arma::vec left = zprob.col(i);
      sumit += sum(left % right);
    }
  } else{
    for(int i = 0; i < zcols; i++){
      double F2F1 = log(R::pbeta(trunc(1), t1vec(i), t2vec(i), 1, 0) - R::pbeta(trunc(0), t1vec(i), t2vec(i), 1, 0));
      arma::vec pp = (dbetaV(xi, t1vec(i), t2vec(i), 1)- F2F1);
      arma::vec right = (lavec(i) + pp);
      arma::vec left = zprob.col(i);
      sumit += sum(left % right);
    }
  }
  return sumit;
}

// Augmented Likelihood Calculation - Beta Distribution
double lnlikecalcB3(const Rcpp::List eout){
  Rcpp::List parmlist = eout["parm.list"];
  arma::vec avec = parmlist["avec"];
  arma::vec t1vec = parmlist["t1vec"];
  arma::vec t2vec = parmlist["t2vec"];
  arma::mat Zprobsmat = eout["zprob"];
  arma::vec xi = eout["xi"];
  arma::vec denom = eout["denom"];
  arma::vec trunc = eout["trunc"];

  const int zcols = avec.size();
  double lnlike = llcalcB3(avec, t1vec, t2vec, zcols, xi, Zprobsmat, denom, trunc);
  return((lnlike*-1.0));
}

// MSTEP
/// Set up for optim
// Defining Ething class
class EthingB3 {
public:
  vec xi;
  int zcols;
  mat zdata;
  vec denom;
  vec trunc;
  EthingB3(vec xi, int zcols, mat zdata, vec denom, vec trunc) : xi(xi), zcols(zcols), zdata(zdata), denom(denom), trunc(trunc){}
};

// Defining the optim function
typedef double optimfn(int n, double *par, void *ex);

// Numerical Optimization - Beta Distribution
//
// description This function is used in nnmin for numeric optimization
// which is necessary for expected maximization to maximize the parameter values.
// param n Number of parameters in par.
// param par The parameters to be optimized.
// param ex Pointer containing all additional information needed.
// return Log likelihood.
double elnlikeB3(int n, double *par, void *ex) {

  EthingB3 *et = (EthingB3 *) ex;

  // zprobsmat
  const int m = et ->zcols;
  arma::mat zprob = et ->zdata;

  // xi
  arma::vec xi = et ->xi;

  // denom
  arma::vec denom = et -> denom;

  // trunc
  arma::vec trunc = et -> trunc;
  double tsum = (trunc(0) + trunc(1));

  // myguess
  arma::vec myguess(n);
  for(int iy = 0; iy<n; iy++){
    myguess(iy) = par[iy];
  }

  // parm
  arma::vec avec(m);
  arma::vec t1vec(m);
  arma::vec t2vec(m);


  if(m == 2){
    avec = invmlogitcD(myguess(0));
  }else{
    avec = invmlogitcU(myguess(span(0, (m-2))));
  }

  t1vec(span(0, (m-3))) = (1.0 + exp(myguess(span((m-1), ((m*2)-4)))));
  t1vec(m-2) = (1.0/(1.0 + exp(-1.0*myguess((m*2)-3))) );
  t1vec(m-1) = (1.0/(1.0 + exp(-1.0*myguess((m*2)-2))) );

  t2vec(span(0, (m-3))) = (1.0 + exp(myguess(span(((m*2)-1), (n-3)))));
  t2vec(m-2) = (1.0/(1.0 + exp(-1.0*myguess(n-2))));
  t2vec(m-1) = (1.0/(1.0 + exp(-1.0*myguess(n-1))));

  double sumit = 0.0;
  arma::vec lavec = log(avec);

  if(tsum == 0){
    for(int i = 0; i < m; i++){
      arma::vec right = (lavec(i) + dbetaV(xi, t1vec(i), t2vec(i), 1));
      arma::vec left = zprob.col(i);
      sumit += sum(left % right);
    }
  } else{
    for(int i = 0; i < m; i++){
      double F2F1 = log(R::pbeta(trunc(1), t1vec(i), t2vec(i), 1, 0) - R::pbeta(trunc(0), t1vec(i), t2vec(i), 1, 0));
      arma::vec pp = (dbetaV(xi, t1vec(i), t2vec(i), 1) - F2F1);
      arma::vec right = (lavec(i) + pp);
      arma::vec left = zprob.col(i);
      sumit += sum(left % right);
    }
  }

  return(sumit*(-1.0));
}




// Defining elnlikeB as an optim function
optimfn elnlikeB3;

// Setting up nmmin
extern "C" {
  void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
             int *fail, double abstol, double intol, void *ex,
             double alpha, double beta, double gamma, int trace,
             int *fncount, int maxit);
}

// M-Step with Numerical Optimization for Expected Maximization - Beta Distribution
//
// description This function is used in expected maximization to maximize the parameter values.
// param eout List with output from the estep
// param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
Rcpp::List mstepB3(Rcpp::List eout){

  const Rcpp::List parmlist = eout["parm.list"];
  const arma::vec avec = parmlist["avec"];
  const arma::vec t1vec = parmlist["t1vec"];
  const arma::vec t2vec = parmlist["t2vec"];
  const arma::mat Zprobsmat = eout["zprob"];
  const arma::mat xi = eout["xi"];
  const arma::vec denom = eout["denom"];
  const arma::vec trunc = eout["trunc"];

  const int m = t1vec.size();

  int gl = ((m-1) + m + m);
  arma::vec guess(gl);

  if(m > 2){
    guess(span(0,(m-2))) = mlogitcU(avec, m);
  } else{
    guess(0) = mlogitcD(avec, m);
  }


  guess(span((m-1), ((m*2)-4))) = log(t1vec(span(0, m-3)) - 1.0);
  guess((m*2)-3) =  (log(t1vec(m-2)) - log(1.0 - t1vec(m-2)) );
  guess((m*2)-2) =  (log(t1vec(m-1)) - log(1.0 - t1vec(m-1)) );

  guess(span(((m*2)-1), ((m*3)-4))) = log(t2vec(span(0, m-3)) - 1.0);
  guess((m*3)-3) =  (log(t2vec(m-2)) - log(1.0 - t2vec(m-2)) );
  guess((m*3)-2) =  (log(t2vec(m-1)) - log(1.0 - t2vec(m-1)) );

  EthingB3 et(xi, m, Zprobsmat, denom, trunc);

  arma::vec vec(gl);
  for(int nn = 0; nn < gl; nn++){
    vec(nn) = guess(nn);
  }

  arma::vec opar(gl);
  double Fmin = 0.0;
  int fail = 0;
  const double abstol = 1.0e-8;
  const double intol = 1.0e-8;
  const double alpha = 1.0;
  const double beta = 0.5;
  const double gamma = 2.0;
  const int trace = 0;
  int fncount = 0;
  const int maxit = 500;

  nmmin(gl, vec.begin(), opar.begin(), &Fmin,
        elnlikeB3, &fail, abstol, intol, &et, alpha, beta,
        gamma, trace, &fncount, maxit);

  // Transform oparback
  arma::vec oparback(gl);
  for(int ou = 0; ou < gl; ou++){
    oparback(ou) = opar[ou];
  }

  arma::vec ahats(m);
  arma::vec t1hats(m);
  arma::vec t2hats(m);

  if(m == 2){
    ahats = invmlogitcD(oparback(0));
  }else{
    ahats = invmlogitcU(oparback(span(0, (m-2))));
  }

  t1hats(span(0, (m-3))) = (1.0 + exp(oparback(span((m-1), ((m*2)-4)))));
  t1hats(m-2) =(1.0/(1.0 + exp(-1.0*oparback((m*2)-3)) ));
  t1hats(m-1) =(1.0/(1.0 + exp(-1.0*oparback((m*2)-2)) ));

  t2hats(span(0, (m-3))) = (1.0 + exp(oparback(span(((m*2)-1), (gl-3)))));
  t2hats(m-2) = (1.0/(1.0 + exp(-1.0*oparback(gl-2)) ));
  t2hats(m-1) = (1.0/(1.0 + exp(-1.0*oparback(gl-1)) ));


  Rcpp::List res = Rcpp::List::create(Rcpp::_["avec"] = ahats, Rcpp::_["t1vec"] = t1hats, Rcpp::_["t2vec"] = t2hats);
  return(res);
}


// Log-Likelihood for BIC
double llcalcfinalB3(Rcpp::List eout){
  Rcpp::List parmlist = eout["parm.list"];
  arma::vec avec = parmlist["avec"];
  arma::vec t1vec = parmlist["t1vec"];
  arma::vec t2vec = parmlist["t2vec"];
  arma::mat xi = eout["xi"];
  arma::vec trunc = eout["trunc"];

  int n = avec.size();
  int xl = xi.size();

  arma::mat mixit(xl, n);
  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
    for(int i = 0; i < n; i++){
      mixit.col(i) =  (avec(i)*(dbetaV(xi, t1vec(i), t2vec(i), 0)));
    }
  }else{
    for(int i  = 0; i < n; i++) {
      double F2F1 = log(R::pbeta(trunc(1), t1vec(i), t2vec(i), 1, 0) - R::pbeta(trunc(0), t1vec(i), t2vec(i), 1, 0));
      arma::vec pp = (avec(i)*(dbetaV(xi, t1vec(i), t2vec(i), 0)));;
      mixit.col(i) =  (pp/F2F1);
    }
  }

  arma::vec sumrows = arma::sum(mixit, 1);
  double sumit = sum(log(sumrows));
  return sumit;
}


//' @title Expected maximization - Beta + Beta + Beta Distribution
//'
//' @description This function is made for the `Bclean()` function and preforms expected maximization with Nelder-Mead
//' numerical optimization for beta distribution.
//'
//' @param parmlist A list containing initial alpha, mean, and variance.
//' @param xi Matrix where the first column is total coverage and the second is the count of base A or B.
//' @param niter Max number of iterates.
//' @param epsilon Epsilon value for convergence tolerance. When the absolute delta log-likelihood is
//'    below this value, convergence is reached.
//' @param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
//'
//' @returns List of elements including the negative log likelihood, the number of iterates,
//'  and the optimized parameter values.
//'
//'
// [[Rcpp::export]]
 Rcpp::List emstepB3(Rcpp::List parmlist, arma::vec xi, int niter, double epsilon, arma::vec trunc){
   Progress p(niter, true);
   Rcpp::List mint = parmlist;
   // E-step
   Rcpp::List eint = estepB3(mint, xi,  trunc);
   double lnlikeint = lnlikecalcB3(eint);

   int count = 0;

   // Iterating E and M steps until convergence is below epsilon
   for(int j = 0; j < (niter + 1); j++){
     count = count+1;
     Rcpp::checkUserInterrupt();
     // M-step
     Rcpp::List moutip = mstepB3(eint);
     p.increment();

     // E-step
     Rcpp::List eoutip = estepB3(moutip, xi, trunc);
     double lnlikeip = lnlikecalcB3(eoutip);
     double deltalogL = (lnlikeint - lnlikeip);
     deltalogL = std::abs(deltalogL);

     // Set next iter
     eint = eoutip;
     lnlikeint = lnlikeip;

     if(deltalogL<epsilon || count > niter){
       break;
     }
   }

   double llreturn = llcalcfinalB3(eint);
   return Rcpp::List::create(Rcpp::_["loglikelihood"] = llreturn, Rcpp::_["negloglikelihood"] = lnlikeint,  Rcpp::_["parm.list"] = eint["parm.list"],  Rcpp::_["niter.done"] = count, Rcpp::_["pir"]=eint["zprob"]);
 }


