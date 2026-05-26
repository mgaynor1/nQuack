# Data Preparation - Use nQuire's Data

This function reduce a three column data frame to two columns by
randomly sampling allele A or B for every site. This is used in our
function
[`process_nquire()`](http://mlgaynor.com/nQuack/reference/process_nquire.md)

## Usage

``` r
nQuire_reformat(xm)
```

## Arguments

- xm:

  A matrix with three columns: Total Coverage, Counts for Allele A, and
  Counts for Allele B.

## Value

Numeric Matrix with total coverage and coverage for a randomly sampled
allele.
