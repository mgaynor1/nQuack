# Expectation maximization - Normal Distribution

This function calculates the log-likelihood using the expectation
maximization algorithm with the Normal-Uniform Distribution. This code
is not identical to nQuire and uses an augmented likelihood.

## Usage

``` r
emstepNUA(parmlist, xi, niter, epsilon, trunc, type = "free")
```

## Arguments

- parmlist:

  A list containing initial alpha, mean, and variance values. The list
  of alpha must include a proportion for the uniform mixture.

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
  "fixed_2" (estimated parameter(s): alpha and variance), or "fixed_3"
  (estimated parameter(s): variance).

## Value

List of elements including the log-likelihood, the number of iterates,
and the optimized parameter values.
