#include <Rcpp.h>
using namespace Rcpp;

// Counts A, C, G, and T
//
// param x Matrix with five columns: Depth, A, C, G, and T.
//
// return Logical Vector indicating where the depth and total count are equal.
LogicalVector countACGTrcpp(NumericMatrix x) {
  int dfrows = x.nrow();
  int dfcols = x.ncol();
  NumericVector ACGTcount(dfrows);
  LogicalVector ACGTcountlog(dfrows);
  for( int i = 0; i < dfrows; i++){
    double tempcount = 0;
    for(int j = 1; j < dfcols; j++){
      tempcount+= x(i, j);
    }
    ACGTcount(i) = tempcount;
    ACGTcountlog(i) = ACGTcount(i) == x(i, 0);
  }
  return(ACGTcountlog);
}

// Classifies Type
//
// param x Matrix with five columns: Depth, A, C, G, and T.
//
// return Logical Vector indicating where type = 2, or only counts for two nucleotides is detected.
LogicalVector typeACGTrcpp(NumericMatrix x) {
  int dfrows = x.nrow();
  int dfcols = x.ncol();
  NumericVector ACGTtype(dfrows);

  for( int i = 0; i < dfrows; i++){
    double tempcount = 0;
    for(int j = 1; j < dfcols; j++){
      double obs = x(i, j);
      if(obs != 0){
        tempcount++;
      }
    }
    ACGTtype(i) = tempcount;
  }
  LogicalVector ACGTtypelog = ACGTtype == 2;
  return(ACGTtypelog);
}

// Depth filter
//
// param x Matrix with five columns: Depth, A, C, G, and T.
// param mindepth Minimum depth, default = 15.
// param maxprob Maximum depth quantile cut off, default = 0.9.
//
// return Logical Vector indicating where depth fits the cut off standards
LogicalVector depthfilter(NumericMatrix xy, int mindepth, double maxprob){
  NumericVector x = xy( _ , 0 );
  Environment stats("package:stats");
  Function quantile = stats["quantile"];
  double maxdepth = as<double>(quantile(x, maxprob));
  LogicalVector lx = (x >= mindepth & x <= maxdepth);
  return(lx);
}

// Logical Vector Sum
//
// param a Logical Vector indicating where the depth and total count are equal.
// param b Logical Vector indicating where type = 2, or only counts for two nucleotides is detected.
// param c Logical Vector indicating where depth fits the cut off standards.
//
// return Numeric Vector equal to the sum of the three logical vectors.
NumericVector addlogic(LogicalVector a, LogicalVector b, LogicalVector c){
  int n = a.size();
  NumericVector abcsum(n);
  for( int i = 0; i < n; i++){
    abcsum(i) = (a(i) + b(i) + c(i));
  }
  return(abcsum);
}

//  Filters Numeric Matrix - EW Goolsby
//
//  @param X input matrix.
//  @param Logical Vector indicating if each row fulfills the condition.
//
//  @return Numeric matrix filtered to only rows that fulfill the condition.
NumericMatrix submat_rcpp(NumericMatrix X, LogicalVector condition) {
  int n=X.nrow();
  int k=X.ncol();
  NumericMatrix out(sum(condition),k);
  for (int i = 0, j = 0; i < n; i++) {
    if(condition[i]) {
      out(j,_) = X(i,_);
      j = j+1;
    }
  }
  return(out);
}


//  Calculates and sums the three logical vectors
//
//  @param x Matrix with five columns: Depth, A, C, G, and T.
//  @param mindepth Minimum depth, default = 15.
//  @param maxprob Maximum depth quantile cut off, default = 0.9.
//
//  @return Numeric Vector equal to the sum of the three logical vectors.
NumericVector firstfilter(NumericMatrix x, int mindepth, double maxprob){
  LogicalVector log1 = depthfilter(x, mindepth, maxprob);
  LogicalVector log2 = typeACGTrcpp(x);
  LogicalVector log3 = countACGTrcpp(x);
  NumericVector abc = addlogic(log1, log2, log3);
  return(abc);
}


//  Make three columns
//
//  @param x Matrix with total coverage, coverage of A, coverage of T, coverage of C, and coverage of G
NumericMatrix makethreecolumns(NumericMatrix x){
  int dfrows = x.nrow();
  NumericMatrix xy(dfrows, 3);
  for(int i = 0; i < dfrows; i++){
    NumericVector y = x(i,_);
    // Rcout << y << std::endl;
    std::sort(y.begin(), y.end(),  std::greater<int>());
    //Rcout << y << std::endl;
    xy(i, _) = y;
  }
  return(xy);
}

//  Filter by Sequencing Error Rate
//
//  @param x Matrix with total coverage, coverage of allele A, and coverage of allele B.
//  @param error Sequencing error rate.
//
NumericMatrix errorfilter(NumericMatrix x, double error){
  if(error != 0.0){
     int dfrows = x.nrow();

    LogicalMatrix refaltlog(dfrows, 2);
    NumericVector refaltsum(dfrows);
    for(int i = 0; i < dfrows; i++){
      NumericVector y = x(i, _);
      double lowend = (error*y(0));
      double highend = ((1-error)*y(0));
      refaltlog(i, 0) = (y(1) > lowend && y(1) < highend);
      refaltlog(i, 1) = (y(2) > lowend && y(2) < highend) ;
      refaltsum(i) = refaltlog(i, 0) + refaltlog(i, 1);
    }
    LogicalVector logicvec = refaltsum == 2;

    NumericMatrix table = submat_rcpp(x, logicvec);
    return table;

  } else{
    return x;
  }

 }

//  Calculates allele frequencies for ref/alt bases
//
//  @param x Matrix with three columbs: Total Coverage, Allele A, Allele B.
//
//  @return List including a two column matrix of allele frequencies and a logical vector indicating
//    if the allele frequencies are less than 0.9, but greater than 0.1.
NumericMatrix refaltcalc(NumericMatrix x, NumericVector trunc){
  if(std::accumulate(trunc.begin(), trunc.end(), 0.0) != 0.0){
      int dfrows = x.nrow();


      NumericMatrix refalt(dfrows, 2);
      for(int i = 0; i < dfrows; i++){
        NumericVector y = x(i, _);
        refalt(i, 0) = (y(1)/y(0));
        refalt(i, 1) = (y(2)/y(0));
      }

      LogicalMatrix refaltlog(dfrows, 2);
      NumericVector refaltsum(dfrows);
      for(int i = 0; i < dfrows; i++){
        NumericVector y = x(i, _);
        refaltlog(i, 0) = (refalt(i, 0) >= trunc[0] && refalt(i, 0) <= trunc[1]);
        refaltlog(i, 1) = (refalt(i, 1) >= trunc[0] && refalt(i, 1) <= trunc[1]) ;
        refaltsum(i) = refaltlog(i, 0) + refaltlog(i, 1);
      }
      LogicalVector logicvec = refaltsum == 2;
      NumericMatrix ready = submat_rcpp(x, logicvec);
      return(ready);
  } else{
      return(x);
  }

}

//  Samples the reference or alternative allele frequency at every site.
//
//  @param x Matrix with two columns: ref and alt.
//
//  @return Numeric vector with one allele frequency per site
NumericMatrix samplefreq(NumericMatrix x){
  int dfrows = x.nrow();
  NumericMatrix freqsamp(dfrows, 2);
  for(int i = 0; i < dfrows; i++){
    freqsamp(i, 0) = x(i, 0);
    NumericVector y(2);
    y[0] = x(i,1);
    y[1] = x(i,2);
    NumericVector prob = NumericVector::create(0.5, 0.5);
    double samp = as<double>(sample(y, 1, false, prob));
    freqsamp(i, 1)  = samp;
  }
  return(freqsamp);
}


//' @name process_rcpp
//' @title Data Preparation - Matrix Filtering
//' @description Based on supplied matrix with total depth and sequencing coverage for each
//'  nucleotide (A, C, G, and T) this function remove all but single nucelotide polymorphisms.
//'  When supplied, this function will filter on coverage or  allele frequency. Finally, the function
//'  samples a single allele frequency per site to avoid data duplication.
//'
//' @param x Matrix with five columns: Depth, A, C, G, and T.
//' @param mindepth Minimum depth, default = 15.
//' @param maxprob Maximum depth quantile cut off, default = 0.9.
//' @param trunc List of two values representing the lower and upper bounds, $c_{L}$ and $c_{U}$.
//' @param error Sequencing error rate.
//'
//' @return Numeric Matrix with total coverage and coverage for a randomly sampled allele.
//'
// [[Rcpp::export]]
NumericMatrix process_rcpp(NumericMatrix x, int mindepth, double maxprob, NumericVector trunc, double error) {
  NumericVector logicsum = firstfilter(x, mindepth, maxprob);
  LogicalVector condition = logicsum == 3;
  NumericMatrix out = submat_rcpp(x, condition);
  NumericMatrix cout = makethreecolumns(out);
  NumericMatrix cond = errorfilter(cout, error);
  NumericMatrix condount = refaltcalc(cond, trunc);
  NumericMatrix ready = samplefreq(condount);
  return(ready);
}
