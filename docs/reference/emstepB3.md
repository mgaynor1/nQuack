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

## Examples

``` r
 if(exists("crazy")){
    xi <- (xm[,2]/xm[,1])
    tcalc <- alphabetacalcvec(mu = c(0.287, 0.50, 0.713),
                              var = c(0.01, 0.01, 0.01))
    set <-  list(avec = c(0.25, 0.25, 0.25, 0.125, 0.125),
                 t1vec = c(tcalc[,1], 0.5, 0.33),
                 t2vec = c(tcalc[,2], 0.33, 0.5))
    checkB <- emstepB3(set,
                       xi,
                       1000,
                       0.1,
                       c(0,0))
}
```
