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
