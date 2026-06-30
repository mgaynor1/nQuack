# Variance calculation from Mean, Tau, and Sequencing Error

This function is used to calculate variance.

## Usage

``` r
muvarcalcvec(mu, tau, error)
```

## Arguments

- mu:

  Vector of means.

- tau:

  Sequence overdispersion parameter for read counts.

- error:

  Sequencing error rate.

## Value

Mean and variance for the associated tau and error.

## Examples

``` r
var <- muvarcalcvec(mu = 0.5, tau = 0.01, error = 0.01)
```
