# Prepare Data - Step 1

This function transforms a BAM file into a text file. Specifically, this
function uses [samtools
mpileup](http://www.htslib.org/doc/samtools-mpileup.md) to translate
your BAM into a tab-separated file. We then filter this file to remove
indels and deletions. When running this function, a temporary folder
will be created (named 'temp/'), however this folder will be removed
once the process is complete.

## Usage

``` r
prepare_data(name, inpath, outpath, tempfolder = "temp")
```

## Arguments

- name:

  File name without the suffix. For example, if your file is called
  "frog.bam", this input should be "frog".

- inpath:

  Location of input file.

- outpath:

  Location for output file.

- tempfolder:

  Location for temp folder.

## Value

Writes text file with the following columns: chromosome, position,
depth, A, C, G, and T.

## Details

Warning, due to the processing time needed for samtools mpileup, this
step may take some time. This function also requires samtools to be
located locally. Please see our [Data
Preparation](https://mlgaynor.com/nQuack/articles/DataPreparation.html)
article for more information. Warning, this writes a temporary folder
titled 'temp'. If you want to run multiple samples at once, we suggest
you set the working directory to separate locations to ensure that your
temp folder/files are not overwritten.

## Examples

``` r
if(exists("crazy")){
## Prepare many samples
  inpath <- "filtered/"
  outpath <- "Processed/"
  filelist <- list.files(path = inpath, pattern = "*.bam" )
  filelist <- gsub(".bam", "", filelist)
  for( i in 1:length(filelist)){
    prepare_data(filelist[i], inpath, outpath)
  }
}
```
