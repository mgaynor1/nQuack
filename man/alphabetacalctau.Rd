% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{alphabetacalctau}
\alias{alphabetacalctau}
\title{Calculate Alpha and Beta from Mean, Tau, and Error rate.}
\usage{
alphabetacalctau(mu, tau, error)
}
\arguments{
\item{mu}{Mean.}

\item{tau}{Overdispersion parameter. Ranges from 0 to 1, where 0 indicates less overdispersion and 1 indicates high overdispersion.  Here tau must be greater than 0.}

\item{error}{Sequencing error rate.}
}
\value{
Numeric vector of alpha and beta.
}
\description{
Calculate Alpha and Beta from Mean, Tau, and Error rate.
}
