// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

#include "emshared.h"
#include "extradistrbbd.h"
#include "betabinomial.h"


// name estepBB
// title E-Step for Expected Maximization - Beta-Binomial Distribution
// description This function is used in expected maximization. Here we complete
//  the E-Step and calculate the log-likelihood. Modifications include a correction for
//  the truncated distribution.
// param parmlist A list containing initial alpha,  $theta_{1}$, and $theta_{2}$ values.
// param xm Matrix where the first column is total coverage and the second is the count of base A or B.
// param type String indicating "Free" or "Fixed".
// param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
// return A list with the following elements:
// \item logL: Log-likelihood of the mixture model.
// \item Zprobs.mat: Matrix where rows represent data points and columns represent each
//   mixture model. Observations are the scaled probability of each data point belonging
//   to each mixture model.
// \item Sj.vec: A vector of column sums for each mixture model.
// \item parm.list: Supplied list containing initial alpha, mu, and sigma values.
// \item  xm: Matrix of coverage and sequencing depth at a site.
Rcpp::List estepBB(const Rcpp::List parmlist, const arma::mat xm, std::string type,
                    const arma::vec trunc){

   const arma::vec avec = parmlist["avec"];
   const arma::vec mvec = parmlist["mvec"];
   const arma::vec svec = parmlist["svec"];

   arma::vec t1vec = t1veccalc(mvec, svec);
   arma::vec t2vec = t2veccalc(mvec, svec);

   const int m = mvec.size();
   const int xl = xm.n_rows;
   arma::mat zprob( xl, m );

   if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0){
     for(int i  = 0; i < m; i++) {
       arma::vec ppv = dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), false);
       zprob.col(i) = (avec(i) * ppv);
     }
   } else{
     for(int i  = 0; i < m; i++) {
       arma::vec F2F1 = pbbinomV((xm.col(0)*trunc(1)), xm.col(0), t1vec(i), t2vec(i)) - pbbinomV((xm.col(0)*trunc(0)), xm.col(0), t1vec(i), t2vec(i));
       arma::vec pp = dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), false);
       pp = (pp/F2F1);
       zprob.col(i) = (avec(i) * pp);
     }
   }

   arma::vec denom = arma::sum(zprob, 1);

   arma::mat zprobb( xl, m );
   for(int i = 0; i < m; i++){
     zprobb.col(i) = ( zprob.col(i) /denom ) ;
   }

   Rcpp::List P = Rcpp::List::create( Rcpp::Named("avec") = avec, Rcpp::_["mvec"] = mvec, Rcpp::_["svec"] = svec);
   return( Rcpp::List::create( Rcpp::_["zprob"] = zprobb,  Rcpp::_["parm.list"] = P, Rcpp::_["xm"] = xm, Rcpp::_["denom"] = denom, Rcpp::_["trunc"] = trunc, Rcpp::_["type"]=type ));

 }

// Augmented Likelihood Calculation - Beta-Binomial Distribution
double llcalcBB(const arma::vec avec, const arma::vec t1vec, const arma::vec t2vec,
                const int zcols, const arma::mat xm,  const arma::mat zprob,
                const arma::vec denom,  const arma::vec trunc) {

   double sumit = 0.0;
   arma::vec lavec = log(avec);

   if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
     for(int i = 0; i < zcols; i++){
       arma::vec right = lavec(i) +(dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), true));
       arma::vec left = zprob.col(i);
       sumit += sum(left % right);
     }
   } else{
     for(int i = 0; i < zcols; i++){
       arma::vec F2F1 = log(pbbinomV((xm.col(0)*trunc(1)), xm.col(0), t1vec(i), t2vec(i)) - pbbinomV((xm.col(0)*trunc(0)), xm.col(0), t1vec(i), t2vec(i)));
       arma::vec pp = (dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), true) - F2F1);
       arma::vec right = (lavec(i) + pp);
       arma::vec left = zprob.col(i);
       sumit += sum(left % right);
      }
    }
   return sumit;
 }

// Augmented Likelihood Calculation - Beta-Binomial Distribution
double lnlikecalcBB(const Rcpp::List eout){
   Rcpp::List parmlist = eout["parm.list"];
   arma::vec avec = parmlist["avec"];
   arma::vec mvec = parmlist["mvec"];
   arma::vec svec = parmlist["svec"];
   arma::mat Zprobsmat = eout["zprob"];
   arma::mat xm = eout["xm"];
   arma::vec denom = eout["denom"];
   std::string type = eout["type"];
   arma::vec trunc = eout["trunc"];

   arma::vec t1vec = t1veccalc(mvec, svec);
   arma::vec t2vec = t2veccalc(mvec, svec);

   const int zcols = avec.size();

   double lnlike = llcalcBB(avec, t1vec, t2vec, zcols, xm, Zprobsmat, denom, trunc);
   //Rprintf("%f \n",lnlike);
   return((lnlike*-1.0));
 }

// MSTEP
/// Set up for optim
// Defining Ething class
class EthingBB {
public:
  mat xm;
  int zcols;
  mat zdata;
  vec denom;
  vec trunc;
  vec parms;
  int pl;
  int flag;

  EthingBB(mat xm, int zcols, mat zdata, vec denom, vec trunc, vec parms, int pl, int flag) : xm(xm), zcols(zcols), zdata(zdata), denom(denom), trunc(trunc), parms(parms), pl(pl), flag(flag) {}
};

// Defining the optim function
typedef double optimfn(int n, double *par, void *ex);

// Numerical Optimization - Beta-Binomial Distribution
//
// description This function is used in nnmin for numeric optimization
// which is necessary for expected maximization to maximize the parameter values.
// param n Number of parameters in par.
// param par The parameters to be optimized.
// param ex Pointer containing all additional information needed.
// return Log likelihood.
double elnlikeBB(int n, double *par, void *ex) {

   EthingBB *et = (EthingBB *) ex;

   // zprobsmat
   const int m = et ->zcols;
   arma::mat zprob = et ->zdata;

   // xi
   arma::mat xm = et ->xm;

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
   arma::vec svec(m);
   arma::vec mvec(m);
   int flag = et ->flag;

   if(flag == 0){
     if(m == 1){
       avec(0) = 1.0;
       mvec(0) = (1.0/(exp(myguess(0)*-1.0)+1.0));
       svec(0) = exp(myguess(1));
     }else {
       if(m == 2){
         avec = invmlogitcD(myguess(0));
       }else{
         avec = invmlogitcU(myguess(span(0, (m-2))));
       }
       mvec = mutransformA(myguess(span((m-1), ((m*2)-2))));
       svec = stransformA(myguess(span(((m*2)-1), (n-1))));
     }
   }else if(flag == 1){ //alpha-free
       if(m == 2){
         avec = invmlogitcD(myguess(0));
       }else{
         avec = invmlogitcU(myguess(span(0, (m-2))));
       }

     arma::vec parms = et -> parms;
     int pl = et -> pl;
     mvec = parms(span(0, (m-1)));
     svec = parms(span(m, (pl-1)));
   }else if(flag == 2){ //alpha-variance-free
     if(m == 1){
       svec(0) = exp(myguess(0));
       arma::vec parms = et -> parms;
       avec(0) = 1.0;
       mvec(0) = parms(1);
     } else{
       if (m == 2){
         avec = invmlogitcD(myguess(0));
       } else{
         avec = invmlogitcU(myguess(span(0, (m-2))));
       }
       svec = stransformA(myguess(span((m-1), (n-1))));
       arma::vec parms = et -> parms;
       int pl = et -> pl;
       mvec = parms(span(0, pl-1));
     }
   } else if(flag == 3){ //variance free
     svec = stransformA(myguess(span(0, (m-1))));
     arma::vec parms = et -> parms;
     int pl = et -> pl;
     avec(span(0, (m-1))) = parms(span(0, (m-1)));
     mvec(span(0, (m-1))) = parms(span(m, (pl-1)));
   }

   arma::vec t1vec = t1veccalc(mvec, svec);
   arma::vec t2vec = t2veccalc(mvec, svec);

   double sumit = 0.0;
   arma::vec lavec = log(avec);

   if(tsum == 0){
     for(int i = 0; i < m; i++){
       arma::vec right = lavec(i) +(dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), true));
       arma::vec left = zprob.col(i);
       sumit += sum(left % right);
     }
   } else{
     for(int i = 0; i < m; i++){
       arma::vec F2F1 = log(pbbinomV((xm.col(0)*trunc(1)), xm.col(0), t1vec(i), t2vec(i)) - pbbinomV((xm.col(0)*trunc(0)), xm.col(0), t1vec(i), t2vec(i)));
       arma::vec pp = (dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), true) - F2F1);
       arma::vec right = (lavec(i) + pp);
       arma::vec left = zprob.col(i);
       sumit += sum(left % right);

     }
   }

   return(sumit*(-1.0));
 }



// Defining elnlikeB as an optim function
optimfn elnlikeBB;

// Setting up nmmin
extern "C" {
   void nmmin(int n, double *xin, double *x, double *Fmin, optimfn fn,
              int *fail, double abstol, double intol, void *ex,
              double alpha, double beta, double gamma, int trace,
              int *fncount, int maxit);
 }

// M-Step with Numerical Optimization for Expected Maximization - Beta-Binomial Distribution
//
// description This function is used in expected maximization to maximize the parameter values.
// param eout List with output from the estep
// param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
Rcpp::List mstepBB(Rcpp::List eout){

   const Rcpp::List parmlist = eout["parm.list"];
   const arma::vec avec = parmlist["avec"];
   const arma::vec mvec = parmlist["mvec"];
   const arma::vec svec = parmlist["svec"];
   const arma::mat Zprobsmat = eout["zprob"];
   const arma::mat xm = eout["xm"];
   const arma::vec denom = eout["denom"];
   std::string type = eout["type"];
   const arma::vec trunc = eout["trunc"];

   const int nmixt = avec.size();
   const int m = mvec.size();

   arma::vec fgp =  sampFT(type, m);
   int flag = fgp(0);
   int gl = fgp(1);
   int pl = fgp(2);

   arma::vec guess(gl);
   arma::vec parms(pl);

   if(flag == 0){
     if(m == 1){
       parms(0) = 1.0;
       guess(0) = (log(mvec(0)) - log(1.0 - mvec(0)));
       guess(1) = log(svec(0));
     }else{
       parms(0) = 0;
       if(m == 2){
         guess(0) = mlogitcD(avec, nmixt);
       } else{
         guess(span(0,(m-2))) = mlogitcU(avec, nmixt);
       }
       guess(span(m-1, ((m*2)-2))) = mvecinA(mvec); // 4-7
       guess(span((m*2)-1, (gl-1))) = svecinA(svec); // 8 -12
     }
   }else if(flag == 1){ //alpha
     if(m == 2){
       guess(0) = mlogitcD(avec, nmixt);
     } else{
       guess(span(0, (m-2))) = mlogitcU(avec, nmixt);
     }
     parms(span(0, (m-1))) = mvec;
     parms(span(m, (pl-1))) = svec;
   }else if(flag == 2){ //both
     if(m == 1){
       guess(0) = log(svec(0));
       parms(0) = 1.0;
       parms(1) = mvec(0);
     } else{
       if(m == 2){
         guess(0) = mlogitcD(avec, nmixt);
       } else{
         guess(span(0,(m-2))) = mlogitcU(avec, nmixt);

       }
       guess(span((m-1), (gl-1))) =  svecinA(svec);
       parms(span(0, (m-1))) = mvec;
      }
   }else if(flag == 3){
     guess(span(0,(m-1))) = svecinA(svec);
     parms(span(0, (m-1))) = avec;
     parms(span(m, (pl-1))) = mvec;
   }


   EthingBB et(xm, nmixt, Zprobsmat, denom, trunc, parms, pl, flag);

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
         elnlikeBB, &fail, abstol, intol, &et, alpha, beta,
         gamma, trace, &fncount, maxit);


   // Transform oparback
   arma::vec oparback(gl);
   for(int ou = 0; ou < gl; ou++){
     oparback(ou) = opar[ou];
   }


   arma::vec ahats(m);
   arma::vec mhats(m);
   arma::vec shats(m);

   if(flag == 0){
     if(m == 1){
       ahats(0) = 1.0;
       mhats(0) = (1.0/(exp(-1.0*oparback(0)) + 1.0));
       shats(0) = exp(oparback(1));
     }else{
       if(m == 2){
         ahats = invmlogitcD(oparback(0));
       }else{
         ahats = invmlogitcU(oparback(span(0, (m-2))));
       }
       mhats = mutransformA(oparback(span((m-1), ((m*2)-2))));
       shats = stransformA(oparback(span(((m*2)-1), (gl-1))));
     }
   }else if(flag == 1){
     if(m == 2){
       ahats = invmlogitcD(oparback(0));
     }else{
       ahats = invmlogitcU(oparback(span(0, (m-2))));
     }
     mhats = mvec;
     shats = svec;
   }else if(flag == 2){
     if(m == 1){
       ahats(0) = 1.0;
       mhats(0) = mvec(0);
       shats(0) =  exp(oparback(0));
     } else {

       if (m == 2){
         ahats = invmlogitcD(oparback(0));
       }else{
         ahats = invmlogitcU(oparback(span(0, (m-2))));
       }
       shats = stransformA(oparback(span((m-1), (m*2)-2)));
       mhats = mvec;

     }
   }else if(flag == 3){
     shats = stransformA(oparback(span(0, (m-1))));
     ahats = avec;
     mhats = mvec;
   }

   // Logic check
   double acheck = ::Rf_fround(arma::min(ahats), 6);
   if(ahats.has_nan() == true || all(acheck) == false ){ //if all are non.zero = false
    ahats = avec;
   }


   Rcpp::List res = Rcpp::List::create(Rcpp::_["avec"] = ahats, Rcpp::_["mvec"] = mhats, Rcpp::_["svec"] = shats);
   return(res);
 }


// Log-Likelihood for BIC
double llcalcfinalBB(Rcpp::List eout){
  Rcpp::List parmlist = eout["parm.list"];
  arma::vec avec = parmlist["avec"];
  arma::vec mvec = parmlist["mvec"];
  arma::vec svec = parmlist["svec"];
  arma::mat xm = eout["xm"];
  arma::vec trunc = eout["trunc"];

  arma::vec t1vec = t1veccalc(mvec, svec);
  arma::vec t2vec = t2veccalc(mvec, svec);

  int n = avec.size();
  int xl = xm.n_rows;

  arma::mat mixit(xl, n);
  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) == 0.0){
    for(int i = 0; i < n; i++){
      mixit.col(i) =  (avec(i)*dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i),false));
    }
  }else{
    for(int i  = 0; i < n; i++) {
      arma::vec F2F1 = (pbbinomV((xm.col(0)*trunc(1)), xm.col(0), t1vec(i), t2vec(i)) - pbbinomV((xm.col(0)*trunc(0)), xm.col(0), t1vec(i), t2vec(i)));
      arma::vec pp = (avec(i)*dbbinomV(xm.col(1), xm.col(0), t1vec(i), t2vec(i), false));
      mixit.col(i) =  (pp/F2F1);
    }
  }

  arma::vec sumrows = arma::sum(mixit, 1);
  double sumit = sum(log(sumrows));
  return sumit;
}


//' @title Expected maximization - Beta-Binomial Distribution
//'
//' @description This function calculates the negative log-likelihood using
//'  the expected maximization algorithm with Nelder-Mead numerical optimization
//'  and beta-binomial distribution.
//'
//' @param parmlist A list containing initial alpha, mean, and variance.
//' @param xm Matrix where the first column is total coverage and the second is the count of base A or B.
//' @param niter Max number of iterates.
//' @param epsilon Epsilon value for convergence tolerance. When the absolute delta log-likelihood is
//'    below this value, convergence is reached.
//' @param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
//' @param type String indicating "Free" or "Fixed".
//'
//' @returns List of elements including the negative log likelihood, the number of iterates,
//'  and the optimized parameter values.
//'
//'
// [[Rcpp::export]]
Rcpp::List emstepBB(Rcpp::List parmlist, arma::mat xm, int niter, double epsilon, arma::vec trunc,  std::string type = "free"){
   Rcpp::List mint = parmlist;
   // E-step
   Rcpp::List eint = estepBB(mint, xm, type, trunc);
   double lnlikeint = lnlikecalcBB(eint);
   //Rprintf("%f \n", lnlikeint);

  int count = 0;

   // Iterating E and M steps until convergence is below epsilon
   for(int j = 0; j < (niter + 1); j++){
    count = count+1;
    Rcpp::checkUserInterrupt();
     // M-step
     Rcpp::List moutip = mstepBB(eint);

     // E-step
     Rcpp::List eoutip = estepBB(moutip, xm, type, trunc);
     double lnlikeip = lnlikecalcBB(eoutip);
     double deltalogL = (lnlikeint - lnlikeip);
     deltalogL = std::abs(deltalogL);
    //Rprintf("%f \n", deltalogL);

     // Set next iter
     eint = eoutip;
     lnlikeint = lnlikeip;

     if(deltalogL<epsilon){
       break;
     }
   }
   double llreturn = llcalcfinalBB(eint);
   return Rcpp::List::create(Rcpp::_["loglikelihood"] = llreturn, Rcpp::_["negloglikelihood"] = lnlikeint,  Rcpp::_["parm.list"] = eint["parm.list"], Rcpp::_["niter.done"] = count);
 }





