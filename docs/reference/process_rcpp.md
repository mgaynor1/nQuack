# Data Preparation - Matrix Filtering

Based on supplied matrix with total depth and sequencing coverage for
each nucleotide (A, C, G, and T) this function remove all but single
nucelotide polymorphisms. When supplied, this function will filter on
coverage or allele frequency. Finally, the function samples a single
allele frequency per site to avoid data duplication.

## Usage

``` r
process_rcpp(x, mindepth, maxprob, trunc, error)
```

## Arguments

- x:

  Matrix with five columns: Depth, A, C, G, and T.

- mindepth:

  Minimum depth, default = 15.

- maxprob:

  Maximum depth quantile cut off, default = 0.9.

- trunc:

  List of two values representing the lower and upper bounds,\\c\_{L}\\
  and \\c\_{U}\\.

- error:

  Sequencing error rate.

## Value

Numeric Matrix with total coverage and coverage for a randomly sampled
allele.
