################################################################################
### Test across-site data aggregation functions
##
## Created on: 2016-06-17
## Author: Kazuki Yoshida
################################################################################

library(testthat)

###
context("### Test 03.AggregateAcrossSites.R")

set.seed(20160505)

###
### Create multi-site data used for testing
################################################################################

## Sample size
n       <- 10^3
## Randomize treatment
alpha0  <- 0
alphaX  <- rep(0,7)
## Randomize binary outcome
beta0   <- 0
betaX   <- rep(0,7)
betaA   <- 0
betaXA  <- rep(0,7)
## Mean event time of 1
lambda   <- 1
lambda_c <- 1
Tmax     <- 1
survParams <- c(lambda, lambda_c, Tmax)
## Duplicates
K <- 3
## Generate 3-center list differing in sizes only
lstDf <- GenerateDistResNet(lstN          = as.list(seq_len(K) * n),
                          lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                          lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
                          lstSurvParams = rep(list(survParams), K))
## Prepare helper variables
set.seed(20160701)
lstSites <- lapply(lstDf, SitePrepareHelperVariables)
## Request data preparation
set.seed(20160701)
drnReady <- RequestSiteDataPreparation(lstDf)


###
### Function to request data preparation
################################################################################

test_that("data preparation is conducted correctly", {

    ## Correct class assignment
    expect_true("DistResNetReady" %in% class(drnReady))
    expect_true("DistResNet" %in% class(drnReady))
    expect_true("list" %in% class(drnReady))

    ## With the same seed, they should match up.
    for (i in seq_along(drnReady)) {
        expect_equal(drnReady[[i]],
                     lstSites[[i]])
    }

})


###
### REGRESSION RESULTS for meta-analyses
################################################################################

test_that("regression results are stacked up correctly", {

    ## Run site specific regression results
    ans <- lapply(lstSites, SiteRegression)
    ## Name sites
    ans <- lapply(seq_along(ans), function(i) {
        cbind(site = rep(i, nrow(ans[[i]])),
              ans[[i]])
    })
    ## Stack up
    ans <- dplyr::bind_rows(ans)

    ## Result df
    expect_equal(RequestSiteRegression(drnReady),
                 ans)

})


###
### SUMMARY-LEVEL DATASETS for summary-level regressions
################################################################################

test_that("summary datasets are stacked up correctly", {

    ## Run site specific regression results
    ans <- lapply(lstSites, SiteSummary)
    ## Name sites
    ans <- lapply(seq_along(ans), function(i) {
        cbind(site = rep(i, nrow(ans[[i]])),
              ans[[i]])
    })
    ## Stack up
    ans <- dplyr::bind_rows(ans)

    ## Result df
    expect_equal(RequestSiteSummary(drnReady),
                 ans)

})


###
### RISKSET-LEVEL DATASETS for case-centered logistic
################################################################################

test_that("riskset datasets are stacked up correctly", {

    ## Run site specific regression results
    ans <- lapply(lstSites, SiteRisksets)
    ## Name sites
    ans <- lapply(seq_along(ans), function(i) {
        cbind(site = rep(i, nrow(ans[[i]])),
              ans[[i]])
    })
    ## Stack up
    ans <- dplyr::bind_rows(ans)

    ## Result df
    resRisk <- RequestSiteRisksets(drnReady)
    expect_equal(resRisk, ans)


    ## Further checking for binary data
    ## Check binary data construction
    ## eDrsBMatch
    expect_equal(as.numeric(subset(resRisk, site == 1 & outcome == "binary" & method == "eDrsBMatch")[,c("events_A0", "events_A1", "riskset_A0", "riskset_A1")]),
                 as.numeric(c(with(subset(drnReady[[1]], eDrsBMatch == 1),
                                   tapply(Y, A, sum)),
                              with(subset(drnReady[[1]], eDrsBMatch == 1),
                                   tapply(Y, A, length)))))
    ## ePsMatch
    expect_equal(as.numeric(subset(resRisk, site == 1 & outcome == "binary" & method == "ePsMatch")[,c("events_A0", "events_A1", "riskset_A0", "riskset_A1")]),
                 as.numeric(c(with(subset(drnReady[[1]], ePsMatch == 1),
                                   tapply(Y, A, sum)),
                              with(subset(drnReady[[1]], ePsMatch == 1),
                                   tapply(Y, A, length)))))

})

###
### INDIVIDUAL-LEVEL DATASETS for regular regressions
################################################################################

test_that("individual-level datasets are stacked up correctly", {

    ## Run site results
    ans <- lapply(lstSites, SiteDataset)
    ## Name sites
    ans <- lapply(seq_along(ans), function(i) {
        cbind(site = rep(i, nrow(ans[[i]])),
              ans[[i]])
    })
    ## Stack up
    ans <- dplyr::bind_rows(ans)

    ## Result df
    expect_equal(RequestSiteDataset(drnReady),
                 ans)

})


test_that("individual-level counterfactual datasets are stacked up correctly", {

    ## Run site results
    ans <- lapply(lstSites, SiteTruth)
    ## Name sites
    ans <- lapply(seq_along(ans), function(i) {
        cbind(site = rep(i, nrow(ans[[i]])),
              ans[[i]])
    })
    ## Stack up
    ans <- dplyr::bind_rows(ans)

    ## Result df
    expect_equal(RequestSiteTruth(drnReady),
                 ans)

})
