% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/04.AnalyzeData.R
\name{AnalyzeSiteRegressionHelper}
\alias{AnalyzeSiteRegressionHelper}
\title{Meta-analyze based on coefficient and variance vector}
\usage{
AnalyzeSiteRegressionHelper(coef, var)
}
\arguments{
\item{coef}{vector of site-specific coefficients}

\item{var}{vector of site-specific variances}
}
\value{
a vector of meta-analysis coefficient and variance.
}
\description{
Perform fixed-effect meta-analysis by inverse variance method, given a coefficient vector and a variance vector.
}
\author{
Kazuki Yoshida
}
