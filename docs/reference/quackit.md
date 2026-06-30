# Model Selection - Based on BIC or Log-Likelihood

This function is for model interpretation.

## Usage

``` r
quackit(
  model_out,
  summary_statistic = "BIC",
  mixtures = c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")
)
```

## Arguments

- model_out:

  Data frame containing, at minimum, columns labeled LL, type, mixture,
  distribution, and BIC.

- summary_statistic:

  May be equal to BIC or LL.

- mixtures:

  Defaults to
  `c("diploid", "triploid", "tetraploid", "hexaploid", "pentaploid")`.

## Value

Returns data frame with the most likely model for each set of mixtures.
Includes the best and second best mixtures, as well as the difference
between the two. We only use BIC or LL to compare within each
distribution and type. To identify the most accurate model, you will
need to compare accuracy across distributions and types using a set of
known samples. The distributions include Normal, Beta, and
Beta-Binomial - each with and without a uniform mixture. The type
indicates which parameters are estimated for the mixtures: all
parameters (`type = 'free'`, only used to calculate delta
log-likelihood), only alpha (`type = 'fixed'`), only alpha and variance
(`type = 'fixed_2'`), and only variance (`type ='fixed_3`) to be
estimated for each mixture.

## Examples

``` r
out <- quackNormal(xm[1:100,], samplename = "sample1", cores = 1)
#>           <(-)__ <(.)__ <(.)__
#>            (___/  (___/  (___/  nQuack-in-progress
#> parallel set to FALSE
#> Free Model Skipped. Log-likelihood ratio will not be included
#> Calculating likelihood of each mixture with a normal distibution.
#> Calculating likelihood of each mixture with a normal+uniform distibution.
goose <- quackit(out)
```
