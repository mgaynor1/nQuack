# Expectation maximization - Beta-Binomial Distribution

This function calculates the negative log-likelihood using the
expectation maximization algorithm with Nelder-Mead numerical
optimization and beta-binomial distribution.

## Usage

``` r
emstepBB(parmlist, xm, niter, epsilon, trunc, type = "free")
```

## Arguments

- parmlist:

  A list containing initial alpha, mean, and variance.

- xm:

  Matrix where the first column is total coverage and the second is the
  count of base A or B.

- niter:

  Max number of iterates.

- epsilon:

  Epsilon value for convergence tolerance. When the absolute delta
  log-likelihood is below this value, convergence is reached.

- trunc:

  List of two values representing the lower and upper bounds, \\c\_{L}\\
  and \\c\_{U}\\.

- type:

  String indicating "Free" or "Fixed".

## Value

List of elements including the negative log likelihood, the number of
iterates, and the optimized parameter values.
