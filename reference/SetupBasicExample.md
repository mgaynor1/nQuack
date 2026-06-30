# Setup Basic Example

This function was made to download all files needed to run the Basic
Example vignette. It downloads a zipped file from Zenodo and unzips it
to the "inst/extdata/" directory within the nQuack package.

## Usage

``` r
SetupBasicExample(overwrite = FALSE)
```

## Arguments

- overwrite:

  Logical. If TRUE, the function will overwrite any existing files in
  the data directory. Default is FALSE.

## Value

This function downloads files to the "inst/extdata/" directory within
the nQuack package and does not return any value. The function will
print messages as it downloads to keep you updated on the progress.

## Examples

``` r
if(exists("crazy")){
   SetupBasicExample(overwrite = FALSE)
}
```
