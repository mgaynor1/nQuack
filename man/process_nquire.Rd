% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_nquire.R
\name{process_nquire}
\alias{process_nquire}
\title{Use nQuire's Data}
\usage{
process_nquire(file)
}
\arguments{
\item{file}{Output text file created with nQuire.}
}
\value{
Numeric matrix with total coverage and coverage for a randomly sampled allele.
}
\description{
If you happen to like nQuire's data preparation more than ours,
uses their data in our program. After processing samples with nQuire's \code{create} and \code{view} functions,
the resulting text file can be read into R. To prepare the data frame for nQuack, we reduce the three
column data frame to two columns by randomly sampling allele A or B for every site.
}
