# Process data - Step 2

Based on the file generated with
[`prepare_data()`](http://mlgaynor.com/nQuack/reference/prepare_data.md),
which contains the total depth and sequencing coverage for each
nucleotide (A, C, G, and T), this function remove all but single
nucelotide polymorphisms. When supplied, this function will filter on
coverage or allele frequency. To filter by total coverage, a user must
supply the `min.depth` and `max.depth.quantile.prob`. If an `error` is
provided, sites will be retained where allele coverage is greater than
the sequencing error rate times the total coverage, but less than one
minus the sequencing error rate times the total coverage. Lastly, based
on `trunc`, allele frequencies will be filtered based on a provided
lower and upper bound. Finally, the function samples a single allele
frequency per site to avoid data duplication.

## Usage

``` r
process_data(
  file,
  min.depth = 2,
  max.depth.quantile.prob = 0.9,
  error = 0.01,
  trunc = c(0, 0)
)
```

## Arguments

- file:

  Output txt file created with
  [`prepare_data()`](http://mlgaynor.com/nQuack/reference/prepare_data.md).

- min.depth:

  Minimum sequencing depth, default as 2.

- max.depth.quantile.prob:

  Maximum sequencing depth quantile cut off, default = 0.9.

- error:

  Sequencing error rate. If an `error` is provided, sites will be
  retained where allele coverage is greater than the sequencing error
  rate times the total coverage, but less than one minus the sequencing
  error rate times the total coverage.

- trunc:

  List of two values representing the lower and upper bounds, \\c\_{L}\\
  and \\c\_{U}\\ which are used to filter allele frequencies.

## Value

Numeric matrix with total coverage and coverage for a randomly sampled
allele.
