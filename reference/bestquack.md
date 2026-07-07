# Model Selection - Expectation Maximization - Optimal Distribution and Type

This function was made to run a subset of models based on a selected
distribution and type. There are many limitations to this function to
make this tractable, as there are 128 models that could be run with our
package. Here we do not include models or comparisons we found
unhelpful, this includes the nQuire implementation and log-likelihood
ratio tests.

## Usage

``` r
bestquack(
  xm,
  distribution,
  type,
  uniform,
  mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid"),
  samplename,
  trunc = c(0, 0),
  lowvar = FALSE,
  tau = NA,
  error = NA,
  verbose = TRUE
)
```

## Arguments

- xm:

  Matrix with two columns with total coverage and coverage for a
  randomly sampled allele.

- distribution:

  May be set to normal, beta, or beta-binomial. We do not include the
  implementation with nQuire.

- type:

  May be equal to fixed, fixed_2, or fixed_3.

- uniform:

  If equal to 1, a uniform mixture is included. If equal to 0, no
  uniform mixture is included.

- mixtures:

  Defaults to
  `c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")`.

- samplename:

  Name of sample to be included in output.

- trunc:

  List of two values representing the lower and upper bounds for allele
  frequency truncation ,c_L and c_U. If allele frequency truncation was
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

- verbose:

  Default to TRUE. If TRUE, progress messages will be printed to the
  console.

## Value

BIC scores and log-likelihood (LL) for the included mixture models. For
BIC, the smallest score is the most likely model. For LL, the largest
score is the most likely model.

## Examples

``` r
 out <- bestquack(xm[1:100,],
                  distribution = "normal",
                  type = "fixed",
                  uniform = 1,
                  samplename = "sample1")
#>           <(.)__ <(.)__ <(-)__
#>            (___/  (___/  (___/  nQuack-in-progress
```
