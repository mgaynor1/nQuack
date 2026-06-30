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

  List of two values representing the lower and upper bounds, \\c\_{L}\\
  and \\c\_{U}\\.

- type:

  String indicating model type. Options: "free" (estimated parameter(s):
  alpha, mean, and variance), "fixed" (estimated parameter(s): alpha),
  "fixed-2" (estimated parameter(s): alpha and variance), or "fixed-3"
  (estimated parameter(s): variance). If avec is length of 1, fixed and
  fixed-3 will not be able to return a log-likelihood.

## Value

List of elements including the log likelihood, the negative log
likelihood, the number of iterates, and the optimized parameter values.

## Examples

``` r
  if(exists("crazy")){
    xi <- (xm[,2]/xm[,1])
    p = list(avec = c(0.11, 0.22, 0.34, 0.22, 0.11),
             mvec = c(0.20, 0.33, 0.50, 0.67, 0.80),
             svec = c(0.01, 0.01, 0.01, 0.01, 0.01));
    mout <- emstepB(p,
                    xi,
                    niter = 100,
                    epsilon = 0.1,
                    trunc = c(0.0,0.0))
}
```
