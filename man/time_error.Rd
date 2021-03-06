% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calculate_error.R
\name{time_error}
\alias{time_error}
\title{Estimate the error in the time estimate}
\usage{
time_error(t = NA, N = Inf, R = Inf, H_0 = 0.5, C = 1, relative = TRUE)
}
\arguments{
\item{t}{Inferred time}

\item{N}{Population Size}

\item{R}{Number of genetic markers}

\item{H_0}{Frequency of heterozygosity at t = 0}

\item{C}{Mean number of crossovers per meiosis (e.g. size in Morgan of
the chromosome)}

\item{relative}{Boolean flag, if TRUE: return the relative error, if FALSE:
return error in generations}
}
\value{
Expected error in the time estimate
}
\description{
Calculate the error in the estimate of the onset of
hybridization, following Equations 3 & 4 in the Supplementary information of
Janzen et al. 2018.
}
