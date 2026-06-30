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

  List of two values representing the lower and upper bounds, \\c\_{L}\\
  and \\c\_{U}\\.

## Value

List of zprob, parm.list, xi, denom, and trunc.

## Examples

``` r
 if(exists("crazy")){
 xi <- (xm[,2]/xm[,1])
 tcalc <- alphabetacalcvec(mu = c(0.287, 0.50, 0.713),
                           var = c(0.01, 0.01, 0.01))
 set <-  list(avec = c(0.25, 0.25, 0.25, 0.125, 0.125),
              t1vec = c(tcalc[,1], 0.5, 0.33),
              t2vec = c(tcalc[,2], 0.33, 0.5))
 checkB <- estepB3(set,
                   xi,
                   c(0,0))
 }
```
