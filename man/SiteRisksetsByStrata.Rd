% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02.ParepareDataWithinSites.R
\name{SiteRisksetsByStrata}
\alias{SiteRisksetsByStrata}
\title{Construct risk set data stratifying on a variable}
\usage{
SiteRisksetsByStrata(time = NULL, event, A, W = NULL, strata = NULL)
}
\arguments{
\item{time}{vector of the observed time variable. If omitted, taken to be all 0. Omit when \code{event} is a binary outcome variable without an accompanying time variable.}

\item{event}{vector of the event status binary variable (must be 0, 1). A binary outcome variable can also be used if the outcome of interest is such a variable.}

\item{A}{vector of the exposure status variable (must be 0, 1).}

\item{W}{vector of the weights. Omit if not weighting.}

\item{strata}{vector of stratifying variable}
}
\value{
reduced stratified risk set data frame having one row for each unique evaluation time for each strata.
}
\description{
Gives a stratified risk set data. See \code{\link{SiteRisksetsHelper}} for details
about the risk set data within each stratum. If not \code{strata} vector is supplied,
one stratum with value NA is assumed.
}
\author{
Kazuki Yoshida
}
