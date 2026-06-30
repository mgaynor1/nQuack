# Model Selection - Expectation Maximization - Beta-Binomial Mixture

This function uses the expectation maximization of both the
beta-binomial and beta-binomial-uniform mixture models for model
selection. Here we can run up to 32 mixture models.

## Usage

``` r
quackBetaBinom(
  xm,
  samplename,
  cores,
  parallel = FALSE,
  trunc = c(0, 0),
  lowvar = FALSE,
  tau = NA,
  error = NA,
  free = FALSE
)
```

## Arguments

- xm:

  Matrix with two columns with total coverage and coverage for a
  randomly sampled allele.

- samplename:

  Name of sample to be included in output.

- cores:

  Threads available to run process in parallel.

- parallel:

  default = FALSE, set to true if cores \> 1.

- trunc:

  List of two values representing the lower and upper bounds for allele
  frequency truncation , c_L and c_U. If allele frequency truncation was
  done to remove error, then you do not need to truncate the expected.
  If no truncation has been done, this should be set to c(0,0), which is
  the default.

- lowvar:

  Default to FALSE. When false, variance is equal to 0.01. If set to
  TRUE and tau and error are not provided, the variance will be set as
  0.001.

- tau:

  Sequencing overdispersion parameter. If tau and error are provided,
  the variance of each mixture will be inferred from these values. If
  not, the variance by default is equal to 0.01 or 0.001.

- error:

  Sequencing error rate. If tau and error are provided, the variance of
  each mixture will be inferred from these values. If not, the variance
  by default is equal to 0.01 or 0.001.

- free:

  default = FALSE, skip the free model calculation and does not
  calculate delta log-likelihood.

## Value

BIC scores and log-likelihood (LL) mixture models including diploid,
triploid, tetraploid, pentaploid, and hexaploid. When free = TRUE, the
delta log-likelihood (dLL) is calculated based on the associated free
model (without or with a uniform mixture). For BIC or delta-log
likelihood, the smallest score is the most likely model. For LL, the
largest score is the most likely model. The type indicates which
parameters are estimated. This function allows all parameters
(`type = 'free'`), only alpha (`type = 'fixed'`), only alpha and
variance (`type = 'fixed_2'`), and only variance (`type ='fixed_3`) to
be estimated for each mixture.

## Examples

``` r
if(exists("crazy")){
  out <- quackBetaBinom(xm[1:100,], samplename = "sample1", cores = 1)
}
```
