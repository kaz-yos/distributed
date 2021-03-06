% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/02.ParepareDataWithinSites.R
\name{SiteRisksetsHelper}
\alias{SiteRisksetsHelper}
\title{Construct risk set data}
\usage{
SiteRisksetsHelper(time = NULL, event, A, W = NULL, eval_times = NULL)
}
\arguments{
\item{time}{vector of the observed time variable. If omitted, taken to be all 0. Omit when \code{event} is a binary outcome variable without an accompanying time variable.}

\item{event}{vector of the event status binary variable (must be 0, 1). A binary outcome variable can also be used if the outcome of interest is such a variable.}

\item{A}{vector of the exposure status variable (must be 0, 1).}

\item{W}{vector of the weights. Omit if not weighting.}

\item{eval_times}{vector of time points to evaluate risk sets and event counts. Time 0 is always included. Required if \code{time} is given. When \code{time} is not given, it is taken to be 0 only.}
}
\value{
reduced risk set data frame having one row for each unique evaluation time.
}
\description{
Construct risk set data classified by exposure status given
the time variable, event status variable, exposure status
variable, as well as time points to evaluate events and risk set
sizes at. The time points for evaluation must include the observed
event times. Otherwise, an interval may contain events and the
riskset size definition at the end of the interval become ambiguous.
As the observed time variable in the simulated dataset is coded in
integer days, specifying all days in the intended follow up duration
will assure inclusion of all event times.
}
\author{
Kazuki Yoshida
}
