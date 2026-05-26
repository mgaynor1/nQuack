# Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error rate.

Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error rate.

## Usage

``` r
alphabetacalctauvec(mu, tau, error)
```

## Arguments

- mu:

  Vector of mean.

- tau:

  Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less
  overdispersion and 1 indicates high overdispersion. Here tau must be
  greater than 0.

- error:

  Sequencing error rate. Ranges from 0 to 1.

## Value

Numeric matrix of alpha and beta.
