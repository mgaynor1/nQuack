# Expectation maximization - Beta Distribution

This function calculates the log-likelihood using the expectation
maximization algorithm with Nelder-Mead numerical optimization and a
beta distribution.

## Usage

``` r
emstepB(parmlist, xi, niter, epsilon, trunc, type = "free")
```

## Arguments

- parmlist:

  A list containing initial alpha, mean, and variance values.

- xi:

  List of observations, in this case allele frequencies.

- niter:

  Max number of iterates.

- epsilon:

  Epsilon value for convergence tolerance. When the absolute delta
  log-likelihood is below this value, convergence is reached.

- trunc:

  List of two values representing the lower and upper bounds, \$c_L\$
  and \$c_U\$.

- type:

  String indicating model type. Options: "free" (estimated parameter(s):
  alpha, mean, and variance), "fixed" (estimated parameter(s): alpha),
  "fixed-2" (estimated parameter(s): alpha and variance), or "fixed-3"
  (estimated parameter(s): variance). If avec is length of 1, fixed and
  fixed-3 will not be able to return a log-likelihood.

## Value

List of elements including the log likelihood, the negative log
likelihood, the number of iterates, and the optimized parameter values.
