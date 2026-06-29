# Expectation maximization - Beta + Beta + Beta Distribution

This function is made for the
[`Bclean()`](http://mlgaynor.com/nQuack/reference/Bclean.md) function
and preforms expectation maximization with Nelder-Mead numerical
optimization for beta distribution.

## Usage

``` r
emstepB3(parmlist, xi, niter, epsilon, trunc)
```

## Arguments

- parmlist:

  A list containing initial alpha, mean, and variance.

- xi:

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

## Value

List of elements including the negative log likelihood, the number of
iterates, and the optimized parameter values.
