# E-Step for Expectation Maximization - Beta + Beta + Beta Distribution

This is used in the
[`Bclean()`](http://mlgaynor.com/nQuack/reference/Bclean.md) function.
Here we complete the E-Step and calculate the log-likelihood.
Modifications include a correction for the truncated distribution.

## Usage

``` r
estepB3(parmlist, xi, trunc)
```

## Arguments

- parmlist:

  A list containing initial alpha, mean, and variance.

- xi:

  List of observations, in this case allele frequencies.

- trunc:

  List of two values representing the lower and upper bounds, c_L and
  c_U.
