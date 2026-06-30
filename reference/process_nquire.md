# Use nQuire's Data

If you happen to like nQuire's data preparation more than ours, uses
their data in our program. After processing samples with nQuire's
`create` and `view` functions, the resulting text file can be read into
R. To prepare the data frame for nQuack, we reduce the three column data
frame to two columns by randomly sampling allele A or B for every site.

## Usage

``` r
process_nquire(file)
```

## Arguments

- file:

  Output text file created with nQuire.

## Value

Numeric matrix with total coverage and coverage for a randomly sampled
allele.

## Examples

``` r
if(file.exists("mynQuirefile.bin")){
  cleaned_data <- process_nquire(file = "mynQuirefile.bin")
}
```
