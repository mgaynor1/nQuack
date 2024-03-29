% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sim.ind.BB.tau.R
\name{sim.ind.BB.tau}
\alias{sim.ind.BB.tau}
\title{Simulate Allele Counts for Single Individual - Beta-Binomial Distribution with Overdispersion and Error}
\usage{
sim.ind.BB.tau(
  mvec,
  avec,
  tau = 0.01,
  error = 0.001,
  s.size = 50000,
  lambda = 11,
  max.coverage = 20,
  min.coverage = 2,
  filter.coverage = TRUE,
  max.depth.quantile.prob = 0.9,
  filter.error = TRUE,
  filter.freq = FALSE,
  trunc = c(0, 0),
  sampled = TRUE
)
}
\arguments{
\item{mvec}{Vector of mean values of allele frequency.}

\item{avec}{Vector of alpha values representing the proportion expected of each mean.}

\item{tau}{Overdispersion parameter. Defaults to 0.01.}

\item{error}{Sequencing error rate. Defaults to 0.001.}

\item{s.size}{Number of biallelic sites to generate. Defaults to 50000. Warning,
the number of sites generated will not be the number of sites returned
due to filtering steps.}

\item{lambda}{Set lambda for the truncated poisson distrubtion. Defaults to 11.}

\item{max.coverage}{Maximum sequencing depth per site. Defaults to 20.}

\item{min.coverage}{Minimum sequencing depth per site. Defaults to 2.}

\item{filter.coverage}{Default as TRUE. Filters to only retain sites where total sequencing depth
is greater than the provided minimum coverage and less than the max quantile depth (set with the max.depth.quantile.prob).}

\item{max.depth.quantile.prob}{Maximum depth quantile probability. Defaults to 0.9.}

\item{filter.error}{Default as TRUE. Filter to only retain sites where allele coverage
is greater than the sequencing error rate times the total coverage, but
less than one minus the sequencing error rate times the total coverage.}

\item{filter.freq}{Default as FALSE. When set to true, sites are filtered based on provided \code{trunc}.}

\item{trunc}{List of two values representing the lower and upper bounds, \eqn{c_{L}} and \eqn{c_{U}}.
Defaults as c(0,0) to represent no truncation.}

\item{sampled}{Default as TRUE. Will randomly sample allele A or allele B, then return a data
frame with total coverage and coverage of a randomly sampled allele will be returned.}
}
\value{
If sampled = FALSE, a data frame with total coverage, coverage of allele A,
and coverage of allele B will be returned. If sampled = TRUE, a data frame with total coverage
and coverage of a randomly sampled allele will be returned.
}
\description{
This function is used to simulate the frequency of
biallelic heterozygous sites assuming a beta-binomial distribution. Here we
sample sequence depth from a truncated poisson distribution between a set minimum, maximum,
and lambda. Only heterozygous sites are returned.
Based on input variables, the sites may be filtered based on the total coverage (\code{filter.coverage}), allele
sequencing coverage (\code{filter.error}), or allele frequency (\code{filter.freq}).
}
