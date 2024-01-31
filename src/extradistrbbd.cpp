#include <Rcpp.h>

// [[Rcpp::interfaces(r, cpp)]]
// [[Rcpp::plugins(cpp11)]]


using Rcpp::NumericVector;

// Note, this is not my code! extraDist
// https://github.com/twolodzko/extraDistr/blob/master/src/beta-binomial-distribution.cpp

// https://github.com/twolodzko/extraDistr/blob/master/src/shared.h

// MACROS

#define GETV(x, i)      x[i % x.length()]    // wrapped indexing of vector
#define GETM(x, i, j)   x(i % x.nrow(), j)   // wrapped indexing of matrix
#define VALID_PROB(p)   ((p >= 0.0) && (p <= 1.0))

//
// functions

inline bool isInteger(double x, bool warn = true);

//https://github.com/twolodzko/extraDistr/blob/master/src/shared.cpp

inline bool isInteger(double x, bool warn) {
  if (ISNAN(x))
    return false;
  if (((x < 0.0) ? std::ceil(x) : std::floor(x)) != x) {
    if (warn) {
      char msg[55];
      std::snprintf(msg, sizeof(msg), "non-integer: %f", x);
      Rcpp::warning(msg);
    }
    return false;
  }
  return true;
}

// https://github.com/twolodzko/extraDistr/blob/master/src/shared_inline.h
inline bool is_large_int(double x) {
  if (x > std::numeric_limits<int>::max())
    return true;
  return false;
}

inline int to_pos_int(double x) {
  if (x < 0.0 || ISNAN(x))
    Rcpp::stop("value cannot be coerced to integer");
  if (is_large_int(x))
    Rcpp::stop("value out of integer range");
  return static_cast<int>(x);
}

inline double to_dbl(int x) {
  return static_cast<double>(x);
}


// https://github.com/twolodzko/extraDistr/blob/master/src/beta-binomial-distribution.cpp

// for dbbinom
inline double logpmf_bbinom(double k, double n, double alpha,
                            double beta, bool& throw_warning) {
#ifdef IEEE_754
  if (ISNAN(k) || ISNAN(n) || ISNAN(alpha) || ISNAN(beta))
    return k+n+alpha+beta;
#endif
  if (alpha < 0.0 || beta < 0.0 || n < 0.0 || !isInteger(n, false)) {
    throw_warning = true;
    return NAN;
  }
  if (!isInteger(k) || k < 0.0 || k > n)
    return R_NegInf;
  // R::choose(n, k) * R::beta(k+alpha, n-k+beta) / R::beta(alpha, beta);
  return R::lchoose(n, k) + R::lbeta(k+alpha, n-k+beta) - R::lbeta(alpha, beta);
}

// for pbbinom
inline std::vector<double> cdf_bbinom_table(double k, double n,
                                            double alpha, double beta) {

  if (k < 0.0 || k > n || alpha < 0.0 || beta < 0.0)
    Rcpp::stop("inadmissible values");

  int ik = to_pos_int(k);
  std::vector<double> p_tab(ik+1);
  double nck, bab, gx, gy, gxy;

  bab = R::lbeta(alpha, beta);
  gxy = R::lgammafn(alpha + beta + n);

  // k = 0

  nck = 0.0;
  gx = R::lgammafn(alpha);
  gy = R::lgammafn(beta + n);
  p_tab[0] = std::exp(nck + gx + gy - gxy - bab);

  if (ik < 1)
    return p_tab;

  // k < 2

  nck += std::log(n);
  gx += std::log(alpha);
  gy -= std::log(n + beta - 1.0);
  p_tab[1] = p_tab[0] + std::exp(nck + gx + gy - gxy - bab);

  if (ik < 2)
    return p_tab;

  // k >= 1

  double dj;

  for (int j = 2; j <= ik; j++) {
    if (j % 10000 == 0)
      Rcpp::checkUserInterrupt();
    dj = to_dbl(j);
    nck += std::log((n + 1.0 - dj)/dj);
    gx += std::log(dj + alpha - 1.0);
    gy -= std::log(n + beta - dj);
    p_tab[j] = p_tab[j-1] + std::exp(nck + gx + gy - gxy - bab);
  }

  return p_tab;
}


NumericVector cpp_dbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& log_prob = false
  ) {

  if (std::min({x.length(), size.length(),
                alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  for (int i = 0; i < Nmax; i++)
    p[i] = logpmf_bbinom(GETV(x, i), GETV(size, i),
                         GETV(alpha, i), GETV(beta, i),
                         throw_warning);

  if (!log_prob)
    p = Rcpp::exp(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}


NumericVector cpp_pbbinom(
    const NumericVector& x,
    const NumericVector& size,
    const NumericVector& alpha,
    const NumericVector& beta,
    const bool& lower_tail = true,
    const bool& log_prob = false
  ) {

  if (std::min({x.length(), size.length(),
                alpha.length(), beta.length()}) < 1) {
    return NumericVector(0);
  }

  int Nmax = std::max({
    x.length(),
    size.length(),
    alpha.length(),
    beta.length()
  });
  NumericVector p(Nmax);

  bool throw_warning = false;

  std::map<std::tuple<int, int, int>, std::vector<double>> memo;

  // maximum modulo size.length(), bounded in [0, size]
  int n = x.length();
  int k = size.length();
  NumericVector mx(k, 0.0);
  for (int i = 0; i < std::max(n, k); i++) {
    if (mx[i % k] < GETV(x, i)) {
      mx[i % k] = std::min(GETV(x, i), GETV(size, i));
    }
  }

  for (int i = 0; i < Nmax; i++) {

    if (i % 100 == 0)
      Rcpp::checkUserInterrupt();

#ifdef IEEE_754
    if (ISNAN(GETV(x, i)) || ISNAN(GETV(size, i)) ||
        ISNAN(GETV(alpha, i)) || ISNAN(GETV(beta, i))) {
      p[i] = GETV(x, i) + GETV(size, i) + GETV(alpha, i) + GETV(beta, i);
      continue;
    }
#endif

    if (GETV(alpha, i) <= 0.0 || GETV(beta, i) <= 0.0 ||
        GETV(size, i) < 0.0 || !isInteger(GETV(size, i), false)) {
      throw_warning = true;
      p[i] = NAN;
    } else if (GETV(x, i) < 0.0) {
      p[i] = 0.0;
    } else if (GETV(x, i) >= GETV(size, i)) {
      p[i] = 1.0;
    } else if (is_large_int(GETV(x, i))) {
      p[i] = NA_REAL;
      Rcpp::warning("NAs introduced by coercion to integer range");
    } else {

      std::vector<double>& tmp = memo[std::make_tuple(
        static_cast<int>(i % size.length()),
        static_cast<int>(i % alpha.length()),
        static_cast<int>(i % beta.length())
      )];

      if (!tmp.size()) {
      //  double mxi = std::min(mx[i % size.length()], GETV(size, i));
        tmp = cdf_bbinom_table(mx[i % size.length()], GETV(size, i), GETV(alpha, i), GETV(beta, i));
      }
      p[i] = tmp[to_pos_int(GETV(x, i))];

    }
  }

  if (!lower_tail)
    p = 1.0 - p;

  if (log_prob)
    p = Rcpp::log(p);

  if (throw_warning)
    Rcpp::warning("NaNs produced");

  return p;
}




