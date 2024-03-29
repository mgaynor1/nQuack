% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{alphabetacalctauvec}
\alias{alphabetacalctauvec}
\title{Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error rate.}
\usage{
alphabetacalctauvec(mu, tau, error)
}
\arguments{
\item{mu}{Vector of mean.}

\item{tau}{Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less overdispersion and 1 indicates high overdispersion.  Here tau must be greater than 0.}

\item{error}{Sequencing error rate. Ranges from 0 to 1.}
}
\value{
Numeric matrix of alpha and beta.
}
\description{
Vector-based - Calculate Alpha and Beta from Mean, Tau, and Error rate.
}
