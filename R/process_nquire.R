#' Use nQuire's Data
#'
#' @description If you happen to like nQuire's data preparation more than ours,
#' uses their data in our program. After processing samples with nQuire's `create` and `view` functions,
#' the resulting text file can be read into R. To prepare the data frame for nQuack, we reduce the three
#' column data frame to two columns by randomly sampling allele A or B for every site.
#'
#' @param file Output text file created with nQuire.
#'
#' @return Numeric matrix with total coverage and coverage for a randomly sampled allele.
#'
#' @importFrom data.table fread

process_nquire <- function(file){
  df <- data.table::fread(file, sep = "\t")
  dfout <- nQuire_reformat(as.matrix(df))
  return(dfout)
}
