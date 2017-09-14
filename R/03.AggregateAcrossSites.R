################################################################################
### Functions for across-site data aggregation
##
## Created on: 2016-06-17
## Author: Kazuki Yoshida
################################################################################

###
### Undocumented helper functions
################################################################################

## Given list of lists containing binary result df and survival result df,
## returns them stacked up with site name column
NameSitesInResults <- function(lstSiteRes) {

    ## Index for sites
    ii <- seq_along(lstSiteRes)

    ## Add site number as a left most column
    ## Loop over sites
    lstOut <- lapply(ii, function(i) {
        ## No need for stringsAsFactors = FALSE because site is numeric.
        cbind(site = rep(i, nrow(lstSiteRes[[i]])),
              lstSiteRes[[i]])
    })
    ## Stack up all sites into one df
    dplyr::bind_rows(lstOut)
}


###
### Function to request data preparation
################################################################################

##' Request sites perform within-site data preparation.
##'
##' Invokes \code{\link{SitePrepareHelperVariables}} on each \code{ResSite} object within a \code{DistResNet} object. Data preparation is performed within each research site so that further analyses can be conducted. This process includes summary score estimation, matching, stratification, and weighting. other \code{RequestSite*} functions can then ask for sharing of relevant data across sites.
##'
##' @param x \code{DistResNet} object representing a distributed research network.
##'
##' @return \code{DistResNetReady} object that has been prepared for other \code{RequestSite*} functions.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteDataPreparation <- function(x) {
    ## Must be an unready DistResNet object
    assert_that("DistResNet" %in% class(x))
    assert_that(!("DistResNetReady" %in% class(x)))

    ## Prepare each site (for to keep class attributes)
    for (i in seq_along(x)) {
        x[[i]] <- SitePrepareHelperVariables(x[[i]])
    }

    class(x) <- c("DistResNetReady", class(x))
    x
}



###
### REGRESSION RESULTS for meta-analyses
################################################################################

##' Request sites provide site-specific regression results
##'
##' Request sites provide site-specific regression results. Aggregated results from all the participating sites are returned.
##'
##' @param x \code{DistResNetReady} object representing a distributed research network.
##'
##' @return data frame containing results for the binary outcome and survival outcome.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteRegression <- function(x) {
    assert_that("DistResNetReady" %in% class(x))

    NameSitesInResults(lapply(x, SiteRegression))
}


###
### SUMMARY-LEVEL DATASETS for summary-level regressions
################################################################################

##' Request sites provide site-specific summary data.
##'
##' Request sites provide site-specific summary data. Aggregated summary data from all the participating sites are returned.
##'
##' @inheritParams RequestSiteRegression
##'
##' @return data frame containing results for the binary outcome and survival outcome.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteSummary <- function(x) {
    assert_that("DistResNetReady" %in% class(x))

    NameSitesInResults(lapply(x, SiteSummary))
}


###
### RISKSET-LEVEL DATASETS for case-centered logistic
################################################################################

##' Request sites provide site-specific risk set data.
##'
##' Request sites provide site-specific risk set data. Aggregated risk set data from all the participating sites are returned.
##'
##' @inheritParams RequestSiteRegression
##'
##' @return data frame containing results for the binary outcome and survival outcome.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteRisksets <- function(x) {
    assert_that("DistResNetReady" %in% class(x))

    NameSitesInResults(lapply(x, SiteRisksets))
}


###
### INDIVIDUAL-LEVEL DATASETS for regular regressions
################################################################################

##' Request sites provide site-specific individual-level data.
##'
##' Request sites provide site-specific individual-level data. Aggregated individual-level data from all the participating sites are returned.
##'
##' @inheritParams RequestSiteRegression
##'
##' @return data frame containing individual-level observed data.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteDataset <- function(x) {
    assert_that("DistResNetReady" %in% class(x))

    NameSitesInResults(lapply(x, SiteDataset))
}


##' Request sites provide site-specific individual-level data (counterfactuals).
##'
##' Request sites provide site-specific individual-level data. Aggregated individual-level data from all the participating sites are returned. This version will request the latent variables including counterfactuals.
##'
##' @inheritParams RequestSiteRegression
##'
##' @return data frame containing individual-level observed data.
##'
##' @author Kazuki Yoshida
##'
##' @export
RequestSiteTruth <- function(x) {
    assert_that("DistResNetReady" %in% class(x))

    NameSitesInResults(lapply(x, SiteTruth))
}
