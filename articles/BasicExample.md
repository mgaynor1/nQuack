# Basic Example

Here we provide an example to show you how to apply this method on real
samples.

Time when reported is based on a MacBook Pro with Apple M2 Max, 12
core-cpu, and 32 GB memory. Not all functions in this package are able
to use multiple CPUs, however we note when functions are able to use
multiple cores. To speed up these analysis, we suggest array SLURM
submissions.

## Download All Data Files

All associated data files need to be downloaded to “inst/extdata/”
folder.

``` r

# Download all data files
SetupBasicExample(overwrite = FALSE)
```

## Data description

After running
[`SetupBasicExample()`](http://mlgaynor.com/nQuack/reference/SetupBasicExample.md),
all data files should be avaliable in “inst/extdata/” folder of this
package. I provide three samples of *Galax urceolata* that were
collected in 2021 as part of my dissertation. These samples have not
been published, but will be submitted to NCBI-SRA in the near future.
This data was generated through target enrichment with species-specific
probes. These samples include a diploid (MLG013), a triploid (MLG015),
and a tetraploid (MLG014).

The raw files were created by following the preprocessing steps outlined
in the [Data
Preparation](https://mlgaynor.com/nQuack/articles/DataPreparation.html).

## Load packages

This tutorial requires `nQuack` and `dplyr`.

``` r

library(nQuack)
library(dplyr)
library(kableExtra)
```

## Data preparation

Please see [Data
Preparation](https://mlgaynor.com/nQuack/articles/DataPreparation.html)
for more information on filtering your data. We also provide a function
to input data from
[nQuire](https://mlgaynor.com/nQuack/reference/process_nquire.html). In
addition, you can use a VCF file for input by following the instructions
found in the [VCF to nQuack
tutorial](https://mlgaynor.com/nQuack/articles/VCF2nQuack.html). Lastly,
you can also [use the output from standardization with
Qploidy2](https://mlgaynor.com/nQuack/articles/Qploidy2nQuack.html) with
nQuack.

### 01. Prepare data

Warning, this step takes a while and results in a file slightly larger
than the input BAM. The time taken in this step can likely be reduce in
the future, however, we keep this file large by conducting no true
filtering steps during this stage. I suggest running this on one CPU
over night on a cluster, though it may take multiple days to finish
depending on the data.

``` r

# Set in and out paths of files
inpath <- "../inst/extdata/01_raw/"
outpath <- "../inst/extdata/02_prepared/"

# List files in the inpath and remove their ending
filelist <- list.files(path = inpath, pattern = "*.bam" )
filelist <- gsub(".bam", "", filelist)

for( i in 1:length(filelist)){
  prepare_data(filelist[i], inpath, outpath)
}
```

### 02. Process data.

Next up, we filter the data file.

This step is fast (2 - 3 seconds per sample) and can be run locally on a
single CPU. Here we are following the filtering approach found in
nQuire: a minimum depth of 10, allele trunctation minimum of 0.15, and
allele truncation maximum of 0.85.

Why are we filtering like this? The most accurate model for this sample
set was under this filtering approach using the normal distribution
implementation with alpha free and a uniform mixture. This model and
filtering approach led to a 97% accuracy with 186 samples ploidal level
correctly assigned.

``` r

inpathtext <- "../inst/extdata/02_prepared/"
newfilelist <- list.files(path = inpathtext, pattern = "*.txt" )

for(i in 1:length(newfilelist)){
    samp <- newfilelist[i]
    temp <- process_data(paste0(inpathtext, samp), 
                         min.depth = 10, 
                         max.depth.quantile.prob = 1, 
                         error = 0.01, 
                         trunc = c(0.15,0.85))
  
  
    write.csv(temp, 
              file = paste0("../inst/extdata/03_processed/", gsub(".txt", "", samp), ".csv"),
              row.names = FALSE)
}
```

## Model inference

### Explore all models

Now we are ready to predict the ploidal level of these samples. When you
are using this method on an unexplored sample set, we suggest to examine
the data under at least 18 model types for all three distributions, a
total of 54 models. These functions can be run on multiple cores.

``` r

samples <- c("MLG013", "MLG014", "MLG015")

for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("../inst/extdata/03_processed/", samples[i], ".csv")))
  out1 <- quackNormal(xm = temp, 
                      samplename = samples[i], 
                      cores = 10, 
                      parallel = TRUE)
  out2 <- quackBeta(xm = temp, 
                    samplename = samples[i], 
                    cores = 10,
                    parallel = TRUE)
  out3 <- quackBetaBinom(xm = temp, 
                         samplename = samples[i],
                         cores = 10, 
                         parallel = TRUE)
  allout <- rbind(out1, out2, out3)
  write.csv(allout, 
            file = paste0("../inst/extdata/04_output/", samples[i], ".csv"),
            row.names = FALSE)
}
```

To run these examples, it took us between 1.46 - 2.09 seconds to run
[`quackNormal()`](http://mlgaynor.com/nQuack/reference/quackNormal.md),
6.41 - 23.16 min to run
[`quackBeta()`](http://mlgaynor.com/nQuack/reference/quackBeta.md), and
3.12 - 27.85 min to run
[`quackBetaBinom()`](http://mlgaynor.com/nQuack/reference/quackBetaBinom.md).
In total, it took 9.54 to 46.15 min to run all models for a sample.

| sample            |     total | normal |    beta | beta.binomial |
|:------------------|----------:|-------:|--------:|--------------:|
| MLG013            |  572.4812 |   1.46 |  384.56 |        186.45 |
| MLG014            | 2769.2193 |   2.09 | 1095.89 |       1671.24 |
| MLG015            | 2380.6057 |   1.81 | 1389.82 |        988.98 |
| Note:             |           |        |         |               |
|  Time in seconds. |           |        |         |               |

### Model interpretation

Using our function
[`quackit()`](http://mlgaynor.com/nQuack/reference/quackit.md), you can
easily interpret model output. Here we are selecting models based on the
BIC score and only considering diploid, triploid, and tetraploid
mixtures.

``` r

inpathtext <- "../inst/extdata/04_output/"
samples <- c("MLG013", "MLG014", "MLG015")

for(i in 1:length(samples)){
  temp <- read.csv(paste0(inpathtext, samples[i], ".csv"))
  summary <- quackit(model_out =  temp, 
                     summary_statistic = "BIC", 
                     mixtures = c("diploid", "triploid", "tetraploid"))
  write.csv(summary, 
            file = paste0("../inst/extdata/05_interpret/", samples[i], ".csv"),
            row.names = FALSE)
}
```

Once we have this output for all our sample, we can pair the outputs
with a key that contains sample names and their ploidal level. To
identify the most accurate model for your data, we tally the accuracy
with some handy `dplyr` functions.

``` r

 # Create key
key <- data.frame(sample = c("MLG013", "MLG014", "MLG015"), 
           ploidal.level = c("diploid", "tetraploid", "triploid"))

# Read in quackit() output
dfs <- lapply(list.files("../inst/extdata/05_interpret/", full.names = TRUE  ), read.csv)
alloutput <- do.call(rbind, dfs)

# Combined
alloutputcombo <- dplyr::left_join(alloutput, key)

# Check the accuracy
alloutputcombo <- alloutputcombo %>%
                  dplyr::mutate(accuracy = ifelse(winnerBIC == ploidal.level, 1, 0))

## What distribution and model type should we use?
sumcheck <- alloutputcombo %>% 
            group_by(Distribution, Type) %>% 
            summarize(total = n(), correct = sum(accuracy))

kbl(sumcheck) %>%
  kable_paper("hover", full_width = F) 
```

| Distribution          | Type    | total | correct |
|:----------------------|:--------|------:|--------:|
| beta                  | fixed   |     3 |       2 |
| beta                  | fixed_2 |     3 |       1 |
| beta                  | fixed_3 |     3 |       2 |
| beta-binomial         | fixed   |     3 |       1 |
| beta-binomial         | fixed_2 |     3 |       2 |
| beta-binomial         | fixed_3 |     3 |       3 |
| beta-binomial-uniform | fixed   |     3 |       1 |
| beta-binomial-uniform | fixed_2 |     3 |       2 |
| beta-binomial-uniform | fixed_3 |     3 |       3 |
| beta-uniform          | fixed   |     3 |       2 |
| beta-uniform          | fixed_2 |     3 |       2 |
| beta-uniform          | fixed_3 |     3 |       2 |
| normal                | fixed   |     3 |       1 |
| normal                | fixed_2 |     3 |       1 |
| normal                | fixed_3 |     3 |       2 |
| normal-uniform        | fixed   |     3 |       3 |
| normal-uniform        | fixed_2 |     3 |       1 |
| normal-uniform        | fixed_3 |     3 |       2 |

For this sample set, we know the normal distribution with alpha free and
a uniform mixture is the most accurate model. So now what?

### Running only the best model

For unknown samples, you can use the best model to predict their ploidal
level with our
[`bestquack()`](http://mlgaynor.com/nQuack/reference/bestquack.md)
function.

``` r

samples <- c("MLG013", "MLG014", "MLG015")
out <- c()

for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("../inst/extdata/03_processed/", samples[i], ".csv")))
  out[[i]] <- bestquack(temp, 
                      distribution = "normal",
                      type = "fixed",
                      uniform = 1, 
                      mixtures = c("diploid", "triploid", "tetraploid"), 
                      samplename = samples[i])
}
```

### Bootstrap replicates

We also provide a function to run bootstrap replicates with the best
model. Warning, this will print a lot of ducks.

``` r

samples <- c("MLG013", "MLG014", "MLG015")
bout <- c()

for(i in 1:length(samples)){
  temp <- as.matrix(read.csv(paste0("../inst/extdata/03_processed/", samples[i], ".csv")))
  bout[[i]] <- quackNboots(temp, 
                        nboots = 100,
                      distribution = "normal",
                      type = "fixed",
                      uniform = 1, 
                      mixtures = c("diploid", "triploid", "tetraploid"), 
                      samplename = samples[i])
}
write.csv(bout[[1]], file = "../inst/extdata/06_boots/MLG013-boots.csv", row.names = FALSE)
write.csv(bout[[2]], file = "../inst/extdata/06_boots/MLG014-boots.csv", row.names = FALSE)
write.csv(bout[[3]], file = "../inst/extdata/06_boots/MLG015-boots.csv", row.names = FALSE)
```

The output of this function includes two rows, the first show the best
model for the original data set, and the second tallies the bootstrap
replicates best model. Here for our diploid, we see the best model for
all replicates is a diploid! This means 100% bootstrap support.

``` r

MLG013boot <- read.csv("../inst/extdata/06_boots/MLG013-boots.csv")
MLG013boot
```

    ##   diploid triploid tetraploid sample
    ## 1       1       NA         NA MLG013
    ## 2     100       NA         NA MLG013

We can use these replicates to identify when a model shouldn’t be
trusted. For example, this model is known to missassign ploidal level to
MLG129, a tetraploid.

``` r

  temp <- as.matrix(read.csv("../inst/extdata/06_boots/MLG129.csv"))
  check  <- quackNboots(temp, 
                        nboots = 1000,
                        distribution = "normal",
                        type = "fixed",
                        uniform = 1, 
                        mixtures = c("diploid", "triploid", "tetraploid"), 
                        samplename = "MLG129")
  
  write.csv(check, file = "../inst/extdata/06_boots/MLG129-boots.csv", row.names = FALSE)
```

Here we found that only 4/1000 bootstrap replicates support the correct
model. Suggesting that any deviation from one mixture/ploidal level may
indicate an untrustworthy model. However, this likely varies across
models and sample sets.

``` r

MLG129boot <- read.csv("../inst/extdata/06_boots/MLG129-boots.csv")
MLG129boot
```

    ##   diploid triploid tetraploid sample
    ## 1       1       NA         NA MLG129
    ## 2     996       NA          4 MLG129
