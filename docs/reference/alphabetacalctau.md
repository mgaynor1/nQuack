# Calculate Alpha and Beta from Mean, Tau, and Error rate.

Calculate Alpha and Beta from Mean, Tau, and Error rate.

## Usage

``` r
alphabetacalctau(mu, tau, error)
```

## Arguments

- mu:

  Mean.

- tau:

  Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less
  overdispersion and 1 indicates high overdispersion. Here tau must be
  greater than 0.

- error:

  Sequencing error rate.

## Value

Numeric vector of alpha and beta.
