# Calculate Variance from Mean, Tau, and Sequencing Error

This function is used to replace variance in mixture model sets.

## Usage

``` r
setconvert(set, tau, error)
```

## Arguments

- set:

  A list of lists, each of the lists must contain avec, mvec, and svec.

- tau:

  Sequence overdispersion parameter for read counts.

- error:

  Sequencing error rate.

## Value

Mean and variance for the associated tau and error.

## Examples

``` r
 set <- c()
 set[[1]] =  list(avec = c(1.00), mvec = c(0.50), svec = c(0.01));
 set[[2]] =  list(avec = c(0.50, 0.50), mvec = c(0.67, 0.33), svec = c(0.01, 0.01));
 exset <- setconvert(set, tau = 0.01, error = 0.001)
```
