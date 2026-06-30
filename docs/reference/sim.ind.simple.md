# Simulate Allele Counts for Single Individual - Simple Approach

This function is used to simulate coverage of each allele at biallelic
heterozygous sites assuming a binomial distribution.

## Usage

``` r
sim.ind.simple(mvec, cover = 100, s.size = 50000, sampled = TRUE)
```

## Arguments

- mvec:

  Vector of means.

- cover:

  Coverage of sites.

- s.size:

  Number of biallelic sites to generate. Defaults to 50000. Warning, the
  number of sites generated will not be the number of sites returned due
  to filtering steps.

- sampled:

  Default as TRUE. Will randomly sample allele A or allele B, then
  return a data frame with total coverage and coverage of a randomly
  sampled allele will be returned.

## Value

If sampled = FALSE, a data frame with total coverage, coverage of allele
A, and coverage of allele B will be returned. If sampled = TRUE, a data
frame with total coverage and coverage of a randomly sampled allele will
be returned.

## Examples

``` r
xm <- sim.ind.simple(mvec = c(0.5))
```
