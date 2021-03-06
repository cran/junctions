% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/estimate_time_haploid.R
\name{estimate_time_haploid}
\alias{estimate_time_haploid}
\title{estimate time using likelihood for a single chromosome}
\usage{
estimate_time_haploid(
  ancestry_matrix,
  N = 1000,
  freq_ancestor_1 = 0.5,
  lower_lim = 2,
  upper_lim = 1000,
  verbose = FALSE
)
}
\arguments{
\item{ancestry_matrix}{matrix with 3 columns, column 1 = chromosome,
column 2 = location in Morgan, column 3 = ancestry.}

\item{N}{Population Size}

\item{freq_ancestor_1}{Frequency of ancestor 1 at t = 0}

\item{lower_lim}{lower limit of the optimization algorithm. Increase if the
expected admixture time is relatively ancient}

\item{upper_lim}{upper limit of the optimization algorithm. If set too large,
recent admixture events can be overlooked - best to set as low as possible.}

\item{verbose}{return verbose output}
}
\value{
The number of generations passed since the onset of hybridization
}
\description{
Estimate the time since the onset of hybridization, for a
haploid genome
}
