#' Process data - Step 2
#'
#' @description Based on the file generated with `prepare_data()`, which contains the total depth and sequencing coverage for each
#'  nucleotide (A, C, G, and T), this function remove all but single nucelotide polymorphisms.
#'  When supplied, this function will filter on coverage or  allele frequency. To filter by total coverage,
#'  a user must supply the `min.depth` and `max.depth.quantile.prob`. If an `error` is provided,
#'  sites will be retained where allele coverage is greater than the sequencing error rate times
#'  the total coverage, but less than one minus the sequencing error rate times the total coverage.
#'  Lastly, based on `trunc`, allele frequencies will be filtered based on a provided lower
#'  and upper bound. Finally, the function samples a single allele frequency per site to avoid data duplication.
#'
#' @param file Output txt file created with `prepare_data()`.
#' @param min.depth Minimum sequencing depth, default as 2.
#' @param max.depth.quantile.prob Maximum sequencing depth quantile cut off, default = 0.9.
#' @param trunc List of two values representing the lower and upper bounds, \eqn{c_{L}} and \eqn{c_{U}} which are used to filter allele frequencies.
#' @param error Sequencing error rate. If an `error` is provided,
#'  sites will be retained where allele coverage is greater than the sequencing error rate times
#'  the total coverage, but less than one minus the sequencing error rate times the total coverage.
#'
#' @return Numeric matrix with total coverage and coverage for a randomly sampled allele.
#' @importFrom data.table fread
#' @importFrom stats na.omit

process_data <- function(file, min.depth = 2, max.depth.quantile.prob = 0.9,
                              error = 0.01, trunc = c(0,0)){
  df <- data.table::fread(file, sep = "\t")
  if(ncol(df) != 7){
    message("Did you run `prepare_data()` on this file? Your data is in the wrong format, please correct and try again.")
  } else{
    df <- df[,3:7]
    df <- stats::na.omit(df)
    dfm <- as.matrix(df)

    allelefreq <- process_rcpp(dfm, min.depth, max.depth.quantile.prob, trunc, error )
    return(allelefreq)
  }
}

