################################################################################
### Test data analysis functions
##
## Created on: 2016-06-17
## Author: Kazuki Yoshida
################################################################################

library(testthat)
library(survey)
library(gnm)
library(dplyr)


###
context("### Test 04.AnalyzeData.R")

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
drn <- GenerateDistResNet(lstN          = as.list(seq_len(K) * n),
                        lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                        lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
                        lstSurvParams = rep(list(survParams), K))
## Prepare helper variables
drnReady <- RequestSiteDataPreparation(drn)
drnReady


###
### Tests for meta-analysis on site-specific regression results
################################################################################

test_that("site-specific regression results are meta-analyzed ok", {

    ## Obtain site specific results
    resReg <- RequestSiteRegression(drnReady)
    ## Inverse variance weighting
    ## https://en.wikipedia.org/wiki/Inverse-variance_weighting
    ans <- resReg %>%
        dplyr::group_by(outcome, method) %>%
        dplyr::summarize(coef = sum(coef * (1/var)) / sum(1/var),
                         var  = 1                   / sum(1/var)) %>%
        as.data.frame

    ## Using split
    ans2 <- resReg %>%
        ## Group by outcome and methods
        split(f = resReg[c("outcome","method")], drop = TRUE) %>%
        lapply(function(df) {
            ## Obtain meta-analyzed results
            res <- AnalyzeSiteRegressionHelper(coef = df$coef,
                                               var  = df$var)
            ## Add labels
            data.frame(outcome = df$outcome[1],
                       method  = df$method[1],
                       coef    = res[1],
                       var     = res[2],
                       stringsAsFactors = FALSE)
        }) %>%
        do.call(what = rbind)
    rownames(ans2) <- NULL
    ans2 <- arrange(ans2, outcome, method)

    ## Check equivalence of dplyr and split method
    expect_equal(ans,
                 ans2)

    ## Add data type
    ans <- cbind(data = rep("meta", nrow(ans)),
                 ans,
                 stringsAsFactors = FALSE)
    ##
    ## Site-specific estimated coefficient augmentation
    ## Spread after dropping variance column
    ans_sites <- tidyr::spread(resReg[,-(which(names(resReg) == "var"))],
                               key = site, value = coef)
    names(ans_sites)[-c(1,2)] <- paste0("coef_site", names(ans_sites)[-c(1,2)])
    ## Augment with site-specific results
    ans_aug <- merge(x = ans,
                     y = ans_sites,
                     by = c("outcome","method"),
                     all.x = TRUE,  all.y = FALSE)
    ##
    ## Site-specific estimated variance augmentation
    ## Spread after dropping var column
    ans_sites2 <- tidyr::spread(resReg[,-(which(names(resReg) == "coef"))],
                                key = site, value = var)
    names(ans_sites2)[-c(1,2)] <- paste0("var_site", names(ans_sites2)[-c(1,2)])
    ## Augment with site-specific results
    ans_aug2 <- merge(x = ans_aug,
                      y = ans_sites2,
                      by = c("outcome","method"),
                      all.x = TRUE,  all.y = FALSE)
    ##
    ## Expectations
    expect_equal(AnalyzeSiteRegression(resReg),
                 ans_aug2)

    ## NA handling by exclusion
    ## df with NA coef/var
    df <- structure(list(
        site = c(1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L,
                 1L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 3L, 3L,
                 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 3L, 4L, 4L, 4L, 4L, 4L, 4L,
                 4L, 4L, 4L, 4L, 4L, 4L),
        outcome = c("binary", "binary", "binary",
                    "binary", "binary", "binary", "survival", "survival", "survival",
                    "survival", "survival", "survival", "binary", "binary", "binary",
                    "binary", "binary", "binary", "survival", "survival", "survival",
                    "survival", "survival", "survival", "binary", "binary", "binary",
                    "binary", "binary", "binary", "survival", "survival", "survival",
                    "survival", "survival", "survival", "binary", "binary", "binary",
                    "binary", "binary", "binary", "survival", "survival", "survival",
                    "survival", "survival", "survival"),
        method = c("ePsMatch", "eDrsBMatch",
                   "ePsStrata", "eDrsBStrata", "ePsSIptw", "ePsMw", "ePsMatch",
                   "eDrsSMatch", "ePsStrata", "eDrsSStrata", "ePsSIptw", "ePsMw",
                   "ePsMatch", "eDrsBMatch", "ePsStrata", "eDrsBStrata", "ePsSIptw",
                   "ePsMw", "ePsMatch", "eDrsSMatch", "ePsStrata", "eDrsSStrata",
                   "ePsSIptw", "ePsMw", "ePsMatch", "eDrsBMatch", "ePsStrata", "eDrsBStrata",
                   "ePsSIptw", "ePsMw", "ePsMatch", "eDrsSMatch", "ePsStrata", "eDrsSStrata",
                   "ePsSIptw", "ePsMw", "ePsMatch", "eDrsBMatch", "ePsStrata", "eDrsBStrata",
                   "ePsSIptw", "ePsMw", "ePsMatch", "eDrsSMatch", "ePsStrata", "eDrsSStrata",
                   "ePsSIptw", "ePsMw"),
        coef = c(-0.0206670211235778, -0.0359256282253753,
                 NA, NA, -0.0324363322777292, -0.0322850842183134, 0.0220560126279101,
                 0.0140049071631618, 0.0167938810788718, 0.0160771752521929, 0.016114372965551,
                 0.0159248054508719, -0.0740122570647937, -0.063164687950672,
                 -0.068026816460471, -0.0681755066950125, -0.0672431267731625,
                 -0.067948177256375, -0.0496525275658786, -0.029953971037688,
                 -0.0330645364343326, -0.0322730605798652, -0.0320033946690058,
                 -0.0315832553432548, 0.0401524896479436, 0.0193015251437017,
                 0.0187398115189868, 0.0183815771183701, 0.0202039288419552, 0.0174674110215562,
                 0.0252265993210464, 0.0166996695087391, 0.0165600159403483, 0.0184443875492624,
                 0.0195358022819071, 0.0217928042991456, -0.0953583925039326,
                 -0.176367486516478, -0.094863780177608, -0.107194261821463, -0.0985059412500863,
                 -0.0933366797533911, -0.0424521823636991, 0.0489053579419892,
                 0.0143387529288193, 0.00386208995080197, 0.00266573694146372,
                 0.000229282078056276),
        var = c(0.000510305479697766, 0.000509610557964662,
                NA, NA, 0.000487305303929011, 0.000487358220388539, 0.000295417046743524,
                0.000295626364462884, 0.00028269790404288, 0.000282732841717295,
                0.000282573934940181, 0.00028259346361012, 0.00274181200078064,
                0.00274673400490608, 0.00247879963494072, 0.00247929757943473,
                0.00247976002267762, 0.00248242918113736, 0.00158551271054998,
                0.00157515971152444, 0.0014214415587498, 0.0014232274146712,
                0.00141395224389523, 0.00141874739338038, 0.00297447268100481,
                0.00296951460328587, 0.00255383536180751, 0.00255569479401723,
                0.00255374619543534, 0.00255587972823593, 0.00166406892071795,
                0.00167095012028587, 0.00144101004948219, 0.00144156168479951,
                0.00143753253512358, 0.00143917897723988, 0.0127198156204243,
                0.0126155486600761, 0.0105985149548199, 0.0105368635094055, 0.0105968213384742,
                0.0106326171090129, 0.00747048341237477, 0.007173658929895, 0.00618275814683362,
                0.00613120550805594, 0.00613153196387081, 0.00614683006386217
                )), .Names = c("site", "outcome", "method", "coef", "var"),
        row.names = c(NA, -48L), class = "data.frame")

    ## df without NA's
    df_complete_rows <- df[!is.na(df$coef) & !is.na(df$var), ]

    expect_equal(AnalyzeSiteRegression(df),
                 AnalyzeSiteRegression(df_complete_rows))
})


test_that("meta-analysis handle NAs gracefully and all methods still return results (which may be NA)", {

    ## Test helper function
    ## Complete
    coefs1 <- c(0.18473589, 0.09948962, 0.06935778)
    vars1  <- c(0.017619107, 0.008654851, 0.005549735)
    expect_equal(AnalyzeSiteRegressionHelper(coefs1, vars1),
                 c(sum(coefs1 * (1/vars1)) / sum(1/vars1),
                   1                       / sum(1/vars1)))
    ## Broken
    coefs2 <- c(0.18473589, NA, 0.06935778)
    vars2  <- c(0.017619107, 0.008654851, NA)
    expect_equal(AnalyzeSiteRegressionHelper(coefs2, vars2),
                 c(coefs2[1], vars2[1]))
    ## Empty
    coefs3 <- c(NA, NA, NA)
    vars3  <- c(NA, NA, NA)
    expect_equal(AnalyzeSiteRegressionHelper(coefs3, vars3),
                 c(NA, NA))


    ## Obtain normal site-specific regression results
    resReg <- RequestSiteRegression(drnReady)

    ## Introduce NAs to binary outcomes for examination
    ## Drop ePsMatch result from site 1 (keep sites 2 and 3)
    resReg[resReg$outcome == "binary" & resReg$site == 1 & resReg$method == "ePsMatch",
           c("coef","var")] <- NA
    ## Drop ePsStrataB coef from site 1 and variance from site 2
    resReg[resReg$outcome == "binary" & resReg$site == 1 & resReg$method == "ePsStrataB",
           "coef"] <- NA
    resReg[resReg$outcome == "binary" & resReg$site == 2 & resReg$method == "ePsStrataB",
           "var"] <- NA
    ## Drop eDrsBMatch results from all sites
    resReg[resReg$outcome == "binary" & resReg$method == "eDrsBMatch",
           c("coef","var")] <- NA

    ## Perform meta-analysis
    res <- AnalyzeSiteRegression(resReg)

    ## Expectations
    ## Site 1 results are NA for ePsMatch
    expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "ePsMatch",
                                c("coef_site1","var_site1")]),
                 as.numeric(c(NA,NA)))
    ## Sites 2 and 3 are ok
    expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "ePsMatch",
                                c("coef_site2","var_site2")]),
                 as.numeric(resReg[resReg$outcome == "binary" & resReg$site == 2 & resReg$method == "ePsMatch",
                                   c("coef","var")]))
    expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "ePsMatch",
                                c("coef_site3","var_site3")]),
                 as.numeric(resReg[resReg$outcome == "binary" & resReg$site == 3 & resReg$method == "ePsMatch",
                                   c("coef","var")]))
    ## Summary results match meta-analysis of sites 2 and 3 excluding invalid site 1
    ans <- resReg %>%
        dplyr::filter(resReg$outcome == "binary",
                      resReg$site %in% c(2,3),
                      resReg$method == "ePsMatch") %>%
        dplyr::group_by(outcome, method) %>%
        dplyr::summarize(coef = sum(coef * (1/var)) / sum(1/var),
                         var  = 1                   / sum(1/var)) %>%
        as.data.frame
        expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "ePsMatch",
                                    c("coef","var")]),
                     as.numeric(ans[,c("coef","var")]))


    ## Summary results match site 3 (only valid site) for ePsStrataB
    expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "ePsStrataB",
                                c("coef","var")]),
                 as.numeric(resReg[resReg$outcome == "binary" & resReg$site == 3 & resReg$method == "ePsStrataB",
                                   c("coef","var")]))


    ## Even without any valid site-specific regression results, results must be returned with NA's
    expect_equal(as.numeric(res[res$outcome == "binary" & res$method == "eDrsBMatch",
                                c("coef","var")]),
                 as.numeric(c(NA,NA)))

})


###
### Tests for regression on site-specific summary data
################################################################################

test_that("summary-to-IPD expander for binary works correctly", {

    ## Data frame of summary data aggregated across sites
    df <- data.frame(site = c(1,1),
                     A = c(0,1),
                     events = c(2,3),
                     denom = c(5, 10))

    ans <- data.frame(site = rep(1, 15),
                      A = rep(c(0,1), c(5,10)),
                      Y = rep(c(1,0,1,0), c(2,3,3,7)))

    ## Expectations
    expect_equal(ExpandToIpd(df),
                 ans)

    ## Using simulated data
    resSummary <- RequestSiteSummary(drnReady)

    ## IPD for binary outco
    binIpd <- ExpandToIpd(subset(resSummary, outcome == "binary"))

    ## Total number check
    ## Matched
    expect_equal(length(subset(binIpd, method == "ePsMatch")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "ePsMatch")$denom))
    expect_equal(sum(subset(binIpd, method == "ePsMatch")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "ePsMatch")$events))
    expect_equal(length(subset(binIpd, method == "eDrsBMatch")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "eDrsBMatch")$denom))
    expect_equal(sum(subset(binIpd, method == "eDrsBMatch")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "eDrsBMatch")$events))
    ## Stratified
    expect_equal(length(subset(binIpd, method == "ePsStrata")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "ePsStrata")$denom))
    expect_equal(sum(subset(binIpd, method == "ePsStrata")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "ePsStrata")$events))
    expect_equal(length(subset(binIpd, method == "eDrsBStrata")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "eDrsBStrata")$denom))
    expect_equal(sum(subset(binIpd, method == "eDrsBStrata")$Y),
                 sum(subset(resSummary, outcome == "binary" & method == "eDrsBStrata")$events))


    ## Test for degenerate data (should expand to degenerate data)
    df_corrupt <- data.frame(site = c(3L, 3L, 3L, 3L, 3L),
                             outcome = c("binary", "binary", "binary", "binary", "binary"),
                             method = c("unadj", "unadj", "ePsMatch", "ePsMatch", "eDrsBMatch"),
                             strata = as.numeric(c(NA, NA, NA, NA, NA)),
                             A = c(0L, 1L, 0L, 1L, NA),
                             events = c(0, 2, 0, 0, NA), denom = c(9667, 10333, 5829, 5829, NA),
                             stringsAsFactors = FALSE)
    df_corrupt[5,]
    expect_equal(ExpandToIpd(df_corrupt[5,]),
                 data.frame(site = 3,
                            outcome = "binary",
                            method = "eDrsBMatch",
                            strata = as.numeric(NA),
                            A = as.numeric(NA),
                            Y = as.numeric(NA),
                            stringsAsFactors = FALSE))

})

test_that("summary data regression is performed correctly", {

    resSummary <- RequestSiteSummary(drnReady)

    ## IPD for binary outco
    binIpd <- ExpandToIpd(subset(resSummary, outcome == "binary"))

    ## Fit models
    ## Unadjusted
    bin_unadjE       <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "unadj"),
                               method = c("exact", "approximate", "efron", "breslow")[3])
    ## Matched
    bin_ePsMatchE    <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "ePsMatch"),
                               method = c("exact", "approximate", "efron", "breslow")[3])
    bin_eDrsBMatchE  <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "eDrsBMatch"),
                               method = c("exact", "approximate", "efron", "breslow")[3])
    ## Stratified
    bin_ePsStrataE   <- clogit(formula = Y ~ A + strata(site, strata),
                               data = subset(binIpd, method == "ePsStrata"),
                               method = c("exact", "approximate", "efron", "breslow")[3])
    bin_eDrsBStrataE <- clogit(formula = Y ~ A + strata(site, strata),
                               data = subset(binIpd, method == "eDrsBStrata"),
                               method = c("exact", "approximate", "efron", "breslow")[3])
    ## Breslow
    ## Unadjusted
    bin_unadjB       <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "unadj"),
                               method = c("exact", "approximate", "efron", "breslow")[4])
    ## Matched
    bin_ePsMatchB    <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "ePsMatch"),
                               method = c("exact", "approximate", "efron", "breslow")[4])
    bin_eDrsBMatchB  <- clogit(formula = Y ~ A + strata(site),
                               data = subset(binIpd, method == "eDrsBMatch"),
                               method = c("exact", "approximate", "efron", "breslow")[4])
    ## Stratified
    bin_ePsStrataB   <- clogit(formula = Y ~ A + strata(site, strata),
                               data = subset(binIpd, method == "ePsStrata"),
                               method = c("exact", "approximate", "efron", "breslow")[4])
    bin_eDrsBStrataB <- clogit(formula = Y ~ A + strata(site, strata),
                               data = subset(binIpd, method == "eDrsBStrata"),
                               method = c("exact", "approximate", "efron", "breslow")[4])

    ## Combined
    ans_bin <- rbind(CoefVar(bin_unadjE, "A"),
                     CoefVar(bin_ePsMatchE, "A"),
                     CoefVar(bin_eDrsBMatchE, "A"),
                     CoefVar(bin_ePsStrataE, "A"),
                     CoefVar(bin_eDrsBStrataE, "A"),
                     ## Breslow
                     CoefVar(bin_unadjB, "A"),
                     CoefVar(bin_ePsMatchB, "A"),
                     CoefVar(bin_eDrsBMatchB, "A"),
                     CoefVar(bin_ePsStrataB, "A"),
                     CoefVar(bin_eDrsBStrataB, "A"))

    ## Rate data (survival data)
    resSurv <- subset(resSummary, outcome == "survival")

    ## Conditional Poisson with gnm (generalized non-linear model)
    ## https://cran.r-project.org/web/packages/gnm/index.html
    ## https://cran.r-project.org/web/packages/gnm/vignettes/gnmOverview.pdf
    ## http://bmcmedresmethodol.biomedcentral.com/articles/10.1186/1471-2288-14-122
    ## Matched without intercept estimation
    surv_unadj       <- gnm(formula = events ~ A,
                            family = poisson(link = "log"),
                            data = subset(resSurv, method == "unadj"),
                            eliminate = factor(site),
                            offset = log(denom))
    surv_ePsMatch    <- gnm(formula = events ~ A,
                            family = poisson(link = "log"),
                            data = subset(resSurv, method == "ePsMatch"),
                            eliminate = factor(site),
                            offset = log(denom))
    surv_eDrsSMatch  <- gnm(formula = events ~ A,
                            family = poisson(link = "log"),
                            data = subset(resSurv, method == "eDrsSMatch"),
                            eliminate = factor(site),
                            offset = log(denom))
    ## Stratified without intercept estimation
    surv_ePsStrata   <- gnm(formula = events ~ A,
                            family = poisson(link = "log"),
                            data = subset(resSurv, method == "ePsStrata"),
                            eliminate = interaction(site, strata),
                            offset = log(denom))
    surv_eDrsSStrata <- gnm(formula = events ~ A,
                            family = poisson(link = "log"),
                            data = subset(resSurv, method == "eDrsSStrata"),
                            eliminate = interaction(site, strata),
                            offset = log(denom))
    ## Combined
    ans_surv <- rbind(CoefVar(surv_unadj, "A"),
                      CoefVar(surv_ePsMatch, "A"),
                      CoefVar(surv_eDrsSMatch, "A"),
                      CoefVar(surv_ePsStrata, "A"),
                      CoefVar(surv_eDrsSStrata, "A"))

    ## Unadjusted with intercept estimation
    surv2_unadj       <- glm(formula = events ~ A + factor(site) + offset(log(denom)),
                             family  = poisson(link = "log"),
                             data    = subset(resSurv, method == "unadj"))
    ## Matched with intercept estimation
    surv2_ePsMatch    <- glm(formula = events ~ A + factor(site) + offset(log(denom)),
                             family  = poisson(link = "log"),
                             data    = subset(resSurv, method == "ePsMatch"))
    surv2_eDrsSMatch  <- glm(formula = events ~ A + factor(site) + offset(log(denom)),
                             family  = poisson(link = "log"),
                             data    = subset(resSurv, method == "eDrsSMatch"))
    ## Statified with intercept estimation
    surv2_ePsStrata   <- glm(formula = events ~ A + interaction(site, strata) + offset(log(denom)),
                             family  = poisson(link = "log"),
                             data    = subset(resSurv, method == "ePsStrata"))
    surv2_eDrsSStrata <- glm(formula = events ~ A + interaction(site, strata) + offset(log(denom)),
                             family  = poisson(link = "log"),
                             data    = subset(resSurv, method == "eDrsSStrata"))

    ## Check equality of point estimates
    expect_equal(coef(surv_unadj)["A"],
                 coef(surv2_unadj)["A"])
    expect_equal(coef(surv_ePsMatch)["A"],
                 coef(surv2_ePsMatch)["A"])
    expect_equal(coef(surv_eDrsSMatch)["A"],
                 coef(surv2_eDrsSMatch)["A"])
    expect_equal(coef(surv_ePsStrata)["A"],
                 coef(surv2_ePsStrata)["A"])
    expect_equal(coef(surv_eDrsSStrata)["A"],
                 coef(surv2_eDrsSStrata)["A"])
    ## Check equality of models
    expect_equal(deviance(surv_unadj)["A"],
                 deviance(surv2_unadj)["A"])
    expect_equal(deviance(surv_ePsMatch),
                 deviance(surv2_ePsMatch))
    expect_equal(deviance(surv_eDrsSMatch),
                 deviance(surv2_eDrsSMatch))
    expect_equal(deviance(surv_ePsStrata),
                 deviance(surv2_ePsStrata))
    expect_equal(deviance(surv_eDrsSStrata),
                 deviance(surv2_eDrsSStrata))

    ## Target object
    ans <- rbind(ans_bin,
                 ans_surv)
    ans0 <- data.frame(data = "summary",
                       outcome = c(rep("binary", nrow(ans_bin)),
                                   rep("survival", nrow(ans_surv))),
                       method = c(paste0(c("unadj",
                                           "ePsMatch", "eDrsBMatch",
                                           "ePsStrata", "eDrsBStrata"),
                                         "E"),
                                  ## Breslow
                                  paste0(c("unadj",
                                           "ePsMatch", "eDrsBMatch",
                                           "ePsStrata", "eDrsBStrata"),
                                         "B"),
                                  ## Survival
                                  "unadj",
                                  "ePsMatch", "eDrsSMatch",
                                  "ePsStrata", "eDrsSStrata"),
                       stringsAsFactors = FALSE)
    ans <- cbind(ans0, ans)

    ## Expectations
    expect_equal(AnalyzeSiteSummary(resSummary),
                 ans)

})


###
### Tests for case-centered conditional logistic regression
################################################################################

## Model specification
## SAS: https://support.sas.com/documentation/cdl/en/statug/63033/HTML/default/viewer.htm#statug_logistic_sect010.htm
##  exposed_events / total_events =  / offset = logodds_exposure
## R: https://stat.ethz.ch/pipermail/r-help/2008-August/170545.html
##  cbind(exposed_events, unexposed_events) ~ 1 + offset(logodds_exposure)

resRisk <- RequestSiteRisksets(drnReady)


test_that("weighted-risk-set-data-to-IPD expanders function correctly", {

    ## Extract simplest one
    df0 <- subset(resRisk, outcome == "binary" & method == "ePsSIptw" & site == 1)
    ## Expand
    out0 <- RisksetsToIpd(df0)
    ## Appropriate columns exist
    expect_true(all(c("site","outcome","method","strata",
                      "start_time", "stop_time", "event", "A", "W") %in% names(out0)))
    ## Weights must be positive
    expect_true(min(out0$W) > 0)
    ## Total size is the sum of risk set sizes
    expect_equal(nrow(out0), sum(df0$riskset_A0, df0$riskset_A1))
    expect_equal(nrow(subset(out0, A == 0)),
                 df0$riskset_A0)
    expect_equal(nrow(subset(out0, A == 1)),
                 df0$riskset_A1)
    ## Event count match up
    expect_equal(sum(subset(out0, A == 0)$event),
                 df0$events_A0)
    expect_equal(sum(subset(out0, A == 1)$event),
                 df0$events_A1)
    ## Mean of weights match up
    expect_equal(mean(subset(out0, A == 0 & event == 0)$W),
                 with(df0, (w_riskset_A0 - w_events_A0) / (riskset_A0 - events_A0)))
    expect_equal(mean(subset(out0, A == 1 & event == 0)$W),
                 with(df0, (w_riskset_A1 - w_events_A1) / (riskset_A1 - events_A1)))
    expect_equal(mean(subset(out0, A == 0 & event == 1)$W),
                 with(df0, w_events_A0 / events_A0))
    expect_equal(mean(subset(out0, A == 1 & event == 1)$W),
                 with(df0, w_events_A1 / events_A1))
    ## Variance of weights match up (may not be exact)
    expect_equal(safe_var(subset(out0, A == 0 & event == 0)$W),
                 with(df0, varw_nonevents_A0))
    expect_equal(safe_var(subset(out0, A == 1 & event == 0)$W),
                 with(df0, varw_nonevents_A1))
    expect_equal(safe_var(subset(out0, A == 0 & event == 1)$W),
                 with(df0, varw_events_A0))
    expect_equal(safe_var(subset(out0, A == 1 & event == 1)$W),
                 with(df0, varw_events_A1))
    ## Compressed format comes out ok
    out0_compress <- out0 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n()) %>%
        as.data.frame
    rownames(out0_compress) <- NULL
    expect_equal(RisksetsToIpd(df0, compress = TRUE),
                 out0_compress)
    ## Skewed distribution version
    expect_equal(RisksetsToIpdSkewed(df0, compress = TRUE),
                 out0_compress)


###  SkewedTwoPointSample: Check for desired skewed two-point (three-point) sample distributions
    ## Call order: RisksetsToIpdSkewed -> RisksetsToIpdSkewedHelper -> SkewedTwoPointSample

    expect_equal(SkewedTwoPointSample(n = df0$events_A0,
                                      m = df0$w_events_A0 / df0$events_A0,
                                      v = df0$varw_events_A0),
                 data.frame(subset(out0_compress, A == 0 & event == 1)[c("W","count")],
                            row.names = NULL))
    expect_equal(SkewedTwoPointSample(n = df0$events_A1,
                                      m = df0$w_events_A1 / df0$events_A1,
                                      v = df0$varw_events_A1),
                 data.frame(subset(out0_compress, A == 1 & event == 1)[c("W","count")],
                            row.names = NULL))
    expect_equal(SkewedTwoPointSample(n = with(df0, (riskset_A0 - events_A0)),
                                      m = with(df0, (w_riskset_A0 - w_events_A0) /
                                                    (riskset_A0 - events_A0)),
                                      v = df0$varw_nonevents_A0),
                 data.frame(subset(out0_compress, A == 0 & event == 0)[c("W","count")],
                            row.names = NULL))
    expect_equal(SkewedTwoPointSample(n = with(df0, (riskset_A1 - events_A1)),
                                      m = with(df0, (w_riskset_A1 - w_events_A1) /
                                                    (riskset_A1 - events_A1)),
                                      v = df0$varw_nonevents_A1),
                 data.frame(subset(out0_compress, A == 1 & event == 0)[c("W","count")],
                            row.names = NULL))


    ## Matching weights (smaller weights; potential for negative regenerated weights)
    df1 <- subset(resRisk, outcome == "binary" & method == "ePsMw" & site == 1)
    ## Expand
    out1 <- RisksetsToIpd(df1)
    ## Appropriate columns exist
    expect_true(all(c("site","outcome","method","strata",
                      "start_time", "stop_time", "event", "A", "W") %in% names(out1)))
    ## Weights must be positive
    expect_true(min(out1$W) > 0)
    ## Total size is the sum of risk set sizes
    expect_equal(nrow(out1), sum(df1$riskset_A0, df1$riskset_A1))
    expect_equal(nrow(subset(out1, A == 0)),
                 df1$riskset_A0)
    expect_equal(nrow(subset(out1, A == 1)),
                 df1$riskset_A1)
    ## Event count match up
    expect_equal(sum(subset(out1, A == 0)$event),
                 df1$events_A0)
    expect_equal(sum(subset(out1, A == 1)$event),
                 df1$events_A1)
    ## Mean of weights match up
    expect_equal(mean(subset(out1, A == 0 & event == 0)$W),
                 with(df1, (w_riskset_A0 - w_events_A0) / (riskset_A0 - events_A0)))
    expect_equal(mean(subset(out1, A == 1 & event == 0)$W),
                 with(df1, (w_riskset_A1 - w_events_A1) / (riskset_A1 - events_A1)))
    expect_equal(mean(subset(out1, A == 0 & event == 1)$W),
                 with(df1, w_events_A0 / events_A0))
    expect_equal(mean(subset(out1, A == 1 & event == 1)$W),
                 with(df1, w_events_A1 / events_A1))
    ## Variance of weights match up (may not be exact)
    expect_equal(safe_var(subset(out1, A == 0 & event == 0)$W),
                 with(df1, varw_nonevents_A0))
    expect_equal(safe_var(subset(out1, A == 1 & event == 0)$W),
                 with(df1, varw_nonevents_A1))
    expect_equal(safe_var(subset(out1, A == 0 & event == 1)$W),
                 with(df1, varw_events_A0))
    expect_equal(safe_var(subset(out1, A == 1 & event == 1)$W),
                 with(df1, varw_events_A1))
    ## Compressed format comes out ok
    out1_compress <- out1 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n()) %>%
        as.data.frame
    expect_equal(RisksetsToIpd(df1, compress = TRUE),
                 out1_compress)
    ## Skewed version
    expect_equal(RisksetsToIpdSkewed(df1, compress = TRUE),
                 out1_compress)


    ## Survival data with matching weights
    ## Risk set data to be expanded
    df2 <- subset(resRisk, outcome == "survival" & method == "ePsMw" & site == 1)
    ## Generate IPD
    out2 <- RisksetsToIpd(df2)
    ## Appropriate columns exist
    expect_true(all(c("site","outcome","method","strata",
                      "start_time", "stop_time", "event", "A", "W") %in% names(out2)))
    ## Weights must be positive
    expect_true(min(out2$W) > 0)
    ## Total size is the sum of risk set sizes
    expect_equal(nrow(out2), sum(df2$riskset_A0, df2$riskset_A1))
    expect_equal(nrow(subset(out2, A == 0)),
                 sum(df2$riskset_A0))
    expect_equal(nrow(subset(out2, A == 1)),
                 sum(df2$riskset_A1))
    ## Event count match up
    expect_equal(sum(subset(out2, A == 0)$event),
                 sum(df2$events_A0))
    expect_equal(sum(subset(out2, A == 1)$event),
                 sum(df2$events_A1))
    ## Compressed format comes out ok
    out2_compress <- out2 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n()) %>%
        as.data.frame
    expect_equal(RisksetsToIpd(df2, compress = TRUE),
                 out2_compress)
    ## Skewed version
    expect_equal(RisksetsToIpdSkewed(df2, compress = FALSE),
                 ## The orders are different. Skewed version is sorted.
                 arrange(out2, start_time, A, event, W))
    ##
    ## Mean of weights match up (first moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out2, A == 0 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df2, (w_riskset_A0 - w_events_A0) / (riskset_A0 - events_A0)))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out2, A == 1 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df2, (w_riskset_A1 - w_events_A1) / (riskset_A1 - events_A1)))
    ##  Events among unexposed (only meaningful among cells where events exist)
    expect_equal(as.numeric(with(subset(out2, A == 0 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df2, (w_events_A0 / events_A0)[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out2, A == 1 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df2, (w_events_A1 / events_A1)[events_A1 > 0]))
    ##
    ## Variance of weights match up (second moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out2, A == 0 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df2, varw_nonevents_A0))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out2, A == 1 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df2, varw_nonevents_A1))
    ##  Events among unexposed
    expect_equal(as.numeric(with(subset(out2, A == 0 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df2, varw_events_A0[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out2, A == 1 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df2, varw_events_A1[events_A1 > 0]))


###  RisksetsToIpdSkewedHelper: Helper for RisksetsToIpdSkewed
    ## No row
    expect_equal(RisksetsToIpdSkewedHelper(n = 0, m = 1, v = 0, A = 0, event = 1),
                 data.frame(A = as.numeric(), event = as.numeric(),
                            W = as.numeric(), count = as.numeric()))
    ## One row
    expect_equal(RisksetsToIpdSkewedHelper(n = 1, m = 1, v = 0, A = 0, event = 1),
                 data.frame(A = 0, event = 1,
                            W = 1, count = 1))
    expect_equal(RisksetsToIpdSkewedHelper(n = 100, m = 1, v = 0, A = 1, event = 0),
                 data.frame(A = 1, event = 0,
                            W = 1, count = 100))
    ## More rows
    ansMoreRows <- cbind(data.frame(A = 1, event = 0)[rep(1,3),],
                         SkewedTwoPointSample(n = 3, m = 1, v = 0.3^2))
    rownames(ansMoreRows) <- NULL
    expect_equal(RisksetsToIpdSkewedHelper(n = 3, m = 1, v = 0.3^2, A = 1, event = 0),
                 ansMoreRows)


###  Use of RisksetsToIpdSkewed
    ## Survival data with matching weights
    ## Risk set data to be expanded
    df3 <- subset(resRisk, outcome == "survival" & method == "ePsMw" & site == 1)
    ## Generate IPD
    out3 <- RisksetsToIpdSkewed(df3, compress = FALSE)
    ## Appropriate columns exist
    expect_true(all(c("site","outcome","method","strata",
                      "start_time", "stop_time", "event", "A", "W") %in% names(out3)))
    ## Weights must be positive
    expect_true(min(out3$W) > 0)
    ## Total size is the sum of risk set sizes
    expect_equal(nrow(out3), sum(df3$riskset_A0, df3$riskset_A1))
    expect_equal(nrow(subset(out3, A == 0)),
                 sum(df3$riskset_A0))
    expect_equal(nrow(subset(out3, A == 1)),
                 sum(df3$riskset_A1))
    ## Event count match up
    expect_equal(sum(subset(out3, A == 0)$event),
                 sum(df3$events_A0))
    expect_equal(sum(subset(out3, A == 1)$event),
                 sum(df3$events_A1))
    ## Recompressing result in the same thing in this simple case not requiring skew
    out3_compress <- out3 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n()) %>%
        as.data.frame
    expect_equal(RisksetsToIpdSkewed(df3, compress = FALSE),
                 arrange(out3, start_time, A, event, W))
    ##
    ## Mean of weights match up (first moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_riskset_A0 - w_events_A0) / (riskset_A0 - events_A0)))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_riskset_A1 - w_events_A1) / (riskset_A1 - events_A1)))
    ##  Events among unexposed (only meaningful among cells where events exist)
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_events_A0 / events_A0)[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_events_A1 / events_A1)[events_A1 > 0]))
    ##
    ## Variance of weights match up (second moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_nonevents_A0))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_nonevents_A1))
    ##  Events among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_events_A0[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_events_A1[events_A1 > 0]))


###  Use of RisksetsToIpdExtreme (Bruce's 1 observation above mean version))
    ## Survival data with matching weights
    ## Risk set data to be expanded
    df3 <- subset(resRisk, outcome == "survival" & method == "ePsMw" & site == 1)
    ## Generate IPD
    out3 <- RisksetsToIpdSkewed(df3, compress = FALSE, helper_fun = RisksetsToIpdExtremeHelper)
    ## Appropriate columns exist
    expect_true(all(c("site","outcome","method","strata",
                      "start_time", "stop_time", "event", "A", "W") %in% names(out3)))
    ## Weights must be positive
    expect_true(min(out3$W) > 0)
    ## Total size is the sum of risk set sizes
    expect_equal(nrow(out3), sum(df3$riskset_A0, df3$riskset_A1))
    expect_equal(nrow(subset(out3, A == 0)),
                 sum(df3$riskset_A0))
    expect_equal(nrow(subset(out3, A == 1)),
                 sum(df3$riskset_A1))
    ## Event count match up
    expect_equal(sum(subset(out3, A == 0)$event),
                 sum(df3$events_A0))
    expect_equal(sum(subset(out3, A == 1)$event),
                 sum(df3$events_A1))
    ## Recompressing result in the same thing in this simple case not requiring skew
    out3_compress <- out3 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n()) %>%
        as.data.frame
    expect_equal(RisksetsToIpdSkewed(df3, compress = FALSE, helper_fun = RisksetsToIpdExtremeHelper),
                 arrange(out3, start_time, A, event, W))
    ##
    ## Mean of weights match up (first moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_riskset_A0 - w_events_A0) / (riskset_A0 - events_A0)))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 0),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_riskset_A1 - w_events_A1) / (riskset_A1 - events_A1)))
    ##  Events among unexposed (only meaningful among cells where events exist)
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_events_A0 / events_A0)[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 1),
                                 tapply(W, start_time, mean))),
                 with(df3, (w_events_A1 / events_A1)[events_A1 > 0]))
    ##
    ## Variance of weights match up (second moment preservation)
    ##  Nonevents among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_nonevents_A0))
    ##  Nonevents among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 0),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_nonevents_A1))
    ##  Events among unexposed
    expect_equal(as.numeric(with(subset(out3, A == 0 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_events_A0[events_A0 > 0]))
    ##  Events among exposed
    expect_equal(as.numeric(with(subset(out3, A == 1 & event == 1),
                                 tapply(W, start_time, safe_var))),
                 with(df3, varw_events_A1[events_A1 > 0]))

})


test_that("compressed weighted IPD work (THEY DO NOT IN R)", {

    ## Extract simplest one
    df0 <- subset(resRisk, outcome == "binary" & method == "ePsSIptw" & site == 1)
    ## Expand
    out0 <- RisksetsToIpd(df0)
    ## Compressed data set with additional frequency weighting
    out0_comp <- out0 %>%
        group_by(site, outcome, method, strata, start_time, stop_time, A, event, W) %>%
        summarize(count = n())

    ## Using expanded data (standard by svycoxph)
    svycoxph0 <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ W, data = out0))
    ## Using expanded data (coxph and robust)
    coxph0 <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                    weights = W,
                    data = out0,
                    robust = TRUE)
    ## Point estimates match up
    expect_equal(CoefVar(svycoxph0, "A")[1],
                 CoefVar(coxph0, "A")[1])
    ## SE fairly close
    expect_equal(CoefVar(svycoxph0, "A")[2],
                 ## Gives robust SE to CoefVar()
                 CoefVar(coxph0, "A")[2],
                 tolerance = 10^-5)

    ## Using compressed data with product weights
    svycoxph1 <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ I(W * count), data = out0_comp))
    ## Even the point estimates differ for some reason
    expect_equal(CoefVar(svycoxph1, "A")[1],
                 CoefVar(svycoxph0, "A")[1],
                 tolerance = 10^-2)
    ## Variance is overestimated as the "sample size" is small
    expect_equal(CoefVar(svycoxph1, "A")[2],
                 CoefVar(svycoxph0, "A")[2],
                 tolerance = 1)

    ## Doing the same using coxph
    coxph1 <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                    weights = I(W * count),
                    data = out0_comp,
                    robust = TRUE)
    ## Point estimates match up
    expect_equal(CoefVar(svycoxph1, "A")[1],
                 CoefVar(coxph1, "A")[1])
    ## SE are very different
    expect_equal(CoefVar(svycoxph1, "A")[2],
                 ## Gives robust SE to CoefVar()
                 CoefVar(coxph1, "A")[2],
                 tolerance = 1)

    ## Trying to use additional weights does not work
    svycoxph2 <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                          ## Use count as frequency weights (?) -> Doesn't work
                          weights = count,
                          ## Use W as sampling weightis
                          design = svydesign(ids = ~ 1, weights = ~ W, data = out0_comp))
    ## Even the point estimates differ for some reason
    expect_equal(CoefVar(svycoxph2, "A")[1],
                 CoefVar(svycoxph0, "A")[1],
                 tolerance = 10^-2)
    ## Variance is overestimated as the "sample size" is small
    expect_equal(CoefVar(svycoxph2, "A")[2],
                 CoefVar(svycoxph0, "A")[2],
                 tolerance = 1)
    ## weights = count seems to be multiplied with W to result in the same
    expect_equal(CoefVar(svycoxph2, "A")[1],
                 CoefVar(svycoxph1, "A")[1],
                 tolerance = 10^-15)
    ## same for variance
    expect_equal(CoefVar(svycoxph2, "A")[2],
                 CoefVar(svycoxph1, "A")[2],
                 tolerance = 10^-15)

})


test_that("SAS-powered functions work", {

    ## Generate a SAS script string
    out_string <- GenerateSasProcPhreg(datafile = "./data/input.csv",
                                       method   = "BRESLOW",
                                       strata   = "site age_cat",
                                       outfile  = "./data/output.csv")
    cat(out_string)
    expect_true(grepl("proc import datafile = './data/input.csv'",
                      out_string))
    expect_true(grepl("/ ties = BRESLOW;",
                      out_string))
    expect_true(grepl("strata site age_cat;",
                      out_string))
    expect_true(grepl("outfile = './data/output.csv'",
                      out_string))


    ## Check methods for sasfit object (manually constructed)
    model1 <- list(coef = c(A = 1.1),
                   vcov = matrix(c(2.2), nrow = 1, ncol = 1,
                                 dimnames = list("A", "A")))
    class(model1) <- c("sasfit", class(model1))
    expect_equal(as.numeric(CoefVar(model1, "A")), c(1.1, 2.2))


    ## Extract simplest one
    df0 <- subset(resRisk, outcome == "survival" & method == "ePsSIptw" & site == 1)
    ## Expand
    out0 <- RisksetsToIpd(df0, compress = TRUE)

    ## Expectations
    if (length(suppressWarnings(system("type sas", intern = TRUE, ignore.stderr = TRUE))) == 0) {
        ## If there is no sas command
        expect_error(CaseCenteredLogisticWeightedSas(df = df0, method = "breslow", strata = "site"),
                     regexp = "sas command not available")

    } else {

        ## If sas command is available, it should give a sasfit object
        out <- CaseCenteredLogisticWeightedSas(df = df0, method = "breslow", strata = "site")
        expect_equal(class(out)[1], "sasfit")

    }
})


test_that("risk sets are constructed correctly (retests for binary)", {

    ## Check binary data construction
    ## Construct risksets for binary data
    site1RisksetsBin <- SiteRisksetsBin(drnReady[[1]])
    ## ePsMatch
    site1ePsMatch <- drnReady[[1]][drnReady[[1]]$ePsMatch == 1,]
    expect_equal(as.numeric(site1RisksetsBin[site1RisksetsBin$method == "ePsMatch",
                                             c("events_A0", "events_A1",
                                               "riskset_A0", "riskset_A1")]),
                 as.numeric(c(tapply(site1ePsMatch$Y, site1ePsMatch$A, sum),
                              tapply(site1ePsMatch$Y, site1ePsMatch$A, length))))
    expect_equal(as.numeric(SiteRisksetsByStrata(
        event  = site1ePsMatch[site1ePsMatch$ePsMatch == 1, "Y"],
        A      = site1ePsMatch[site1ePsMatch$ePsMatch == 1, "A"])[,c("events_A0", "events_A1",
                                                                     "riskset_A0", "riskset_A1")]),
        as.numeric(c(tapply(site1ePsMatch$Y, site1ePsMatch$A, sum),
                     tapply(site1ePsMatch$Y, site1ePsMatch$A, length))))
    expect_equal(as.numeric(SiteRisksetsHelper(
        event  = site1ePsMatch[site1ePsMatch$ePsMatch == 1, "Y"],
        A      = site1ePsMatch[site1ePsMatch$ePsMatch == 1, "A"])[,c("events_A0", "events_A1",
                                                                     "riskset_A0", "riskset_A1")]),
        as.numeric(c(tapply(site1ePsMatch$Y, site1ePsMatch$A, sum),
                     tapply(site1ePsMatch$Y, site1ePsMatch$A, length))))
    ## eDrsMatch
    site1eDrsBMatch <- drnReady[[1]][drnReady[[1]]$eDrsBMatch == 1,]
    expect_equal(as.numeric(site1RisksetsBin[site1RisksetsBin$method == "eDrsBMatch",
                                             c("events_A0", "events_A1",
                                               "riskset_A0", "riskset_A1")]),
                 as.numeric(c(tapply(site1eDrsBMatch$Y, site1eDrsBMatch$A, sum),
                              tapply(site1eDrsBMatch$Y, site1eDrsBMatch$A, length))))
    expect_equal(as.numeric(SiteRisksetsByStrata(
        event  = site1eDrsBMatch[site1eDrsBMatch$eDrsBMatch == 1, "Y"],
        A      = site1eDrsBMatch[site1eDrsBMatch$eDrsBMatch == 1, "A"])[,c("events_A0", "events_A1",
                                                                           "riskset_A0", "riskset_A1")]),
        as.numeric(c(tapply(site1eDrsBMatch$Y, site1eDrsBMatch$A, sum),
                     tapply(site1eDrsBMatch$Y, site1eDrsBMatch$A, length))))
    ## ePsStrata
    expect_equal(as.numeric(as.matrix(
        site1RisksetsBin[site1RisksetsBin$method == "ePsStrata",
                         c("events_A0", "events_A1",
                           "riskset_A0", "riskset_A1")])),
        as.numeric(
            cbind(t(tapply(drnReady[[1]]$Y, drnReady[[1]][,c("A","ePsStrata")], sum)),
                  t(tapply(drnReady[[1]]$Y, drnReady[[1]][,c("A","ePsStrata")], length)))))
    ## eDrsBStrata
    expect_equal(as.numeric(as.matrix(
        site1RisksetsBin[site1RisksetsBin$method == "eDrsBStrata",
                         c("events_A0", "events_A1",
                           "riskset_A0", "riskset_A1")])),
        as.numeric(
            cbind(t(tapply(drnReady[[1]]$Y, drnReady[[1]][,c("A","eDrsBStrata")], sum)),
                  t(tapply(drnReady[[1]]$Y, drnReady[[1]][,c("A","eDrsBStrata")], length)))))
})


test_that("Breslow coxph and case-centered logistic are equivalent", {

    ## Equivalence of case-centered logistic and Cox
    ## Data with no tied times (not even in censored)
    data0 <- as.data.frame(drnReady[[3]][,c("time","event","A")])
    data0 <- data0[!duplicated(data0$time),]
    ## Fit models
    cox0_efron   <- coxph(formula = Surv(time, event) ~ A, data = data0, ties = "efron")
    cox0_breslow <- coxph(formula = Surv(time, event) ~ A, data = data0, ties = "breslow")
    cox0_exact   <- coxph(formula = Surv(time, event) ~ A, data = data0, ties = "exact")
    ## Efron is affected by tied times even if no tied event times exist
    expect_equal(coef(cox0_efron), coef(cox0_breslow), tolerance = 10^(-10))
    expect_equal(coef(cox0_efron), coef(cox0_exact),   tolerance = 10^(-10))
    expect_equal(coef(cox0_exact), coef(cox0_breslow), tolerance = 10^(-10))
    ## Risk set construction
    data0_riskset <- with(data0, SiteRisksetsByStrata(time = time, event = event, A = A))
    ## No event row at day 0 does not count
    ccl0 <- CaseCenteredLogistic(data0_riskset)
    expect_equal(CoefVar(ccl0, "(Intercept)"),
                 CoefVar(CaseCenteredLogistic(data0_riskset[-1,]), "(Intercept)"),
                 tolerance = 10^(-10))
    ## Match up to efron when without ties
    expect_equal(as.numeric(CoefVar(ccl0, "(Intercept)")[1]),
                 as.numeric(coef(cox0_efron)["A"]),     tolerance = 10^(-10))
    expect_equal(as.numeric(CoefVar(ccl0, "(Intercept)")[2]),
                 as.numeric(vcov(cox0_efron)["A","A"]), tolerance = 10^(-9))

    ## Data with no tied events (censored can ties)
    data1 <- as.data.frame(drnReady[[3]][,c("time","event","A")])
    data1_A0 <- data1[data1$A == 0,]
    data1_A1 <- data1[data1$A == 1,]
    ## Drop duplicated event times
    data1_A1 <- data1_A1[!duplicated(data1_A1$time),]
    data1 <- rbind(data1_A0, data1_A1)
    ## Fit models
    cox1_efron   <- coxph(formula = Surv(time, event) ~ A, data = data1, ties = "efron")
    cox1_breslow <- coxph(formula = Surv(time, event) ~ A, data = data1, ties = "breslow")
    cox1_exact   <- coxph(formula = Surv(time, event) ~ A, data = data1, ties = "exact")
    ## Efron is affected by tied times even if no tied event times exist
    expect_true(abs(coef(cox1_efron) - coef(cox1_breslow)) > 10^(-4))
    expect_true(abs(coef(cox1_efron) - coef(cox1_exact))   > 10^(-4))
    expect_equal(coef(cox1_exact), coef(cox1_breslow), tolerance = 10^(-5))
    ## Risk set construction
    data1_riskset <- with(data1, SiteRisksetsByStrata(time = time, event = event, A = A))
    ## No event row at day 0 does not count
    ccl1 <- CaseCenteredLogistic(data1_riskset)
    expect_equal(CoefVar(ccl1, "(Intercept)"),
                 CoefVar(CaseCenteredLogistic(data1_riskset[-1,]), "(Intercept)"),
                 tolerance = 10^(-10))
    ## Do not match with Efron
    expect_true(abs(CoefVar(ccl1, "(Intercept)")[1] - as.numeric(coef(cox1_efron)["A"]))     > 10^(-5))
    expect_true(abs(CoefVar(ccl1, "(Intercept)")[2] - as.numeric(vcov(cox1_efron)["A","A"])) > 10^(-7))
    ## Match with Breslow with some numerical tolerance
    expect_equal(as.numeric(CoefVar(ccl1, "(Intercept)")[1]),
                 as.numeric(coef(cox1_breslow)["A"]),     tolerance = 10^(-10))
    expect_equal(as.numeric(CoefVar(ccl1, "(Intercept)")[2]),
                 as.numeric(vcov(cox1_breslow)["A","A"]), tolerance = 10^(-7))
    ## Do not match with exact
    expect_true(abs(CoefVar(ccl1, "(Intercept)")[1] - as.numeric(coef(cox1_exact)["A"]))     > 10^(-6))
    expect_true(abs(CoefVar(ccl1, "(Intercept)")[2] - as.numeric(vcov(cox1_exact)["A","A"])) > 10^(-7))

    ## Data with tied events
    data2 <- data.frame(time  = c(1,2,2,3,4,5,6,7,8,9,10,10),
                        event = c(0,1,1,0,1,0,1,0,1,0,1,0),
                        A     = c(0,1,1,1,1,1,0,0,0,0,1,0))
    cox2_efron   <- coxph(formula = Surv(time, event) ~ A, data = data2, ties = "efron")
    cox2_breslow <- coxph(formula = Surv(time, event) ~ A, data = data2, ties = "breslow")
    cox2_exact   <- coxph(formula = Surv(time, event) ~ A, data = data2, ties = "exact")
    ## With ties, approximation methods do not agree
    expect_true(abs(coef(cox2_efron)   - coef(cox2_breslow)) > 10^(-5))
    expect_true(abs(coef(cox2_efron)   - coef(cox2_exact))   > 10^(-5))
    expect_true(abs(coef(cox2_breslow) - coef(cox2_exact))   > 10^(-5))
    ## Risk set construction
    data2_riskset <- with(data2, SiteRisksetsByStrata(time = time, event = event, A = A))
    ## No event row at day 0 does not count
    ccl2 <- CaseCenteredLogistic(data2_riskset)
    expect_equal(CoefVar(ccl2, "(Intercept)"),
                 CoefVar(CaseCenteredLogistic(data2_riskset[-1,]), "(Intercept)"))
    ## Do not match with Efron
    expect_true(abs(CoefVar(ccl2, "(Intercept)")[1] - as.numeric(coef(cox2_efron)["A"]))     > 10^(-3))
    expect_true(abs(CoefVar(ccl2, "(Intercept)")[2] - as.numeric(vcov(cox2_efron)["A","A"])) > 10^(-3))
    ## Match with Breslow with some numerical tolerance
    expect_equal(as.numeric(CoefVar(ccl2, "(Intercept)")[1]),
                 as.numeric(coef(cox2_breslow)["A"]),     tolerance = 10^(-10))
    expect_equal(as.numeric(CoefVar(ccl2, "(Intercept)")[2]),
                 as.numeric(vcov(cox2_breslow)["A","A"]), tolerance = 10^(-7))
    ## Do not match with exact
    expect_true(abs(CoefVar(ccl2, "(Intercept)")[1] - as.numeric(coef(cox2_exact)["A"]))     > 10^(-3))
    expect_true(abs(CoefVar(ccl2, "(Intercept)")[2] - as.numeric(vcov(cox2_exact)["A","A"])) > 10^(-3))

    ## Data with many ties
    data3 <- as.data.frame(drnReady[[3]][,c("time","event","A")])
    ## Check there are ties
    expect_true(any(table(subset(data3, event == 1)$time) > 1))
    ## With ties, approximation methods do not agree
    cox3_efron   <- coxph(formula = Surv(time, event) ~ A, data = data3, ties = "efron")
    cox3_breslow <- coxph(formula = Surv(time, event) ~ A, data = data3, ties = "breslow")
    cox3_exact   <- coxph(formula = Surv(time, event) ~ A, data = data3, ties = "exact")
    ## With ties, approximation methods do not agree
    expect_true(abs(coef(cox3_efron)   - coef(cox3_breslow)) > 10^(-5))
    expect_true(abs(coef(cox3_efron)   - coef(cox3_exact))   > 10^(-5))
    expect_true(abs(coef(cox3_breslow) - coef(cox3_exact))   > 10^(-5))
    ## Risk set construction
    data3_riskset <- with(data3, SiteRisksetsByStrata(time = time, event = event, A = A))
    ## No event row at day 0 does not count
    ccl3 <- CaseCenteredLogistic(data3_riskset)
    expect_equal(CoefVar(ccl3, "(Intercept)"),
                 CoefVar(CaseCenteredLogistic(data3_riskset[-1,]), "(Intercept)"),
                 tolerance = 10^(-10))
    ## Do not match with Efron
    expect_true(abs(CoefVar(ccl3, "(Intercept)")[1] - as.numeric(coef(cox3_efron)["A"]))     > 10^(-5))
    expect_true(abs(CoefVar(ccl3, "(Intercept)")[2] - as.numeric(vcov(cox3_efron)["A","A"])) > 10^(-9))
    ## Match with Breslow with some numerical tolerance
    expect_equal(as.numeric(CoefVar(ccl3, "(Intercept)")[1]),
                 as.numeric(coef(cox3_breslow)["A"]),     tolerance = 10^(-10))
    expect_equal(as.numeric(CoefVar(ccl3, "(Intercept)")[2]),
                 as.numeric(vcov(cox3_breslow)["A","A"]), tolerance = 10^(-10))
    ## Do not match with Efron
    expect_true(abs(CoefVar(ccl3, "(Intercept)")[1] - as.numeric(coef(cox3_exact)["A"])) > 10^(-5))
    expect_true(abs(CoefVar(ccl3, "(Intercept)")[2] - as.numeric(vcov(cox3_exact)["A","A"])) > 10^(-6))

    ## Data with all ties (binary data)
    data4 <- as.data.frame(drnReady[[3]][,c("Y","A")])
    ## With ties, approximation methods do not agree (exact is not feasible)
    cox4_efron   <- coxph(formula = Surv(rep(1, length(Y)), Y) ~ A, data = data4, ties = "efron")
    cox4_breslow <- coxph(formula = Surv(rep(1, length(Y)), Y) ~ A, data = data4, ties = "breslow")
    ## With ties, approximation methods do not agree
    expect_true(abs(coef(cox4_efron)   - coef(cox4_breslow)) > 10^(-2))
    ## Risk set construction
    data4_riskset <- with(data4, SiteRisksetsByStrata(event = Y, A = A))
    expect_equal(as.numeric(tapply(data4$Y, data4$A, sum)),
                 as.numeric(data4_riskset[,c("events_A0","events_A1")]))
    expect_equal(as.numeric(tapply(data4$Y, data4$A, length)),
                 as.numeric(data4_riskset[,c("riskset_A0","riskset_A1")]))
    ## Case-centered approach
    ccl4 <- CaseCenteredLogistic(data4_riskset)
    ## Do not match with Efron
    expect_true(abs(as.numeric(CoefVar(ccl4, "(Intercept)"))[1] - as.numeric(coef(cox4_efron)["A"]))     > 10^(-2))
    expect_true(abs(as.numeric(CoefVar(ccl4, "(Intercept)"))[2] - as.numeric(vcov(cox4_efron)["A","A"])) > 10^(-9))
    ## Match with Breslow with some numerical tolerance
    expect_equal(as.numeric(CoefVar(ccl4, "(Intercept)"))[1], as.numeric(coef(cox4_breslow)["A"]),     tolerance = 10^(-10))
    expect_equal(as.numeric(CoefVar(ccl4, "(Intercept)"))[2], as.numeric(vcov(cox4_breslow)["A","A"]), tolerance = 10^(-10))

})


test_that("Breslow svycoxph and weighted riskset analyses are equivalent", {

###  MW without ties
    ## Weighted individual data (drop duplicated times)
    data0 <- as.data.frame(drnReady[[3]][,c("time","event","A","ePsMw")])
    data0 <- data0[!duplicated(data0$time),]
    dim(data0)

    ## IPD weighted analyses
    svycoxph0 <- svycoxph(formula = Surv(time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ ePsMw, data = data0),
                          method = "breslow")
    coxph0 <- coxph(formula = Surv(time, event) ~ A,
                    data = data0,
                    method = "breslow",
                    weights = ePsMw,
                    robust = TRUE)

    ## Riskset construction with weight summaries
    risksets0 <- SiteRisksetsByStrata(time = data0[, "time"],
                                      event = data0[, "event"],
                                      A = data0[, "A"],
                                      W = data0[, "ePsMw"])
    ## IPD reconstruction from risksets with weight summaries
    risksets0$site <- 1
    risksets0$outcome <- "survival"
    risksets0$method <- "ePsMw"
    risksets0$strata <- NA
    risksets0_expanded <- RisksetsToIpdSkewed(risksets0, compress = FALSE)

    ## Reconstructed-IPD weighted analyses
    svycoxph0_expanded <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                                   design = svydesign(ids = ~ 1, weights = ~ W, data = risksets0_expanded),
                                   method = "breslow")
    coxph0_expanded <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                             data = risksets0_expanded,
                             method = "breslow",
                             weights = W,
                             robust = TRUE)

    ## Coefficients should match up exactly
    ## Two fitting methods
    expect_equal(coef(svycoxph0), coef(coxph0))
    ## original weighted IPD vs expaneded weighted IPD
    expect_equal(coef(svycoxph0), coef(svycoxph0_expanded))
    expect_equal(coef(coxph0), coef(coxph0_expanded))

    ## Variances
    ## Two robust variance methods
    ## Finite-population correction results in larger variance with MW, but close
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(coxph0)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^4)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^5)
    ## Variance: original weighted IPD > expaneded weighted IPD
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(svycoxph0_expanded)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^3)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^4)
    ##
    expect_true(as.numeric(vcov(coxph0)) > as.numeric(vcov(coxph0_expanded)))
    expect_true(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^3)
    expect_false(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^4)


###  MW with ties
    ## Weighted individual data (drop duplicated times)
    data0 <- as.data.frame(drnReady[[3]][,c("time","event","A","ePsMw")])

    ## IPD weighted analyses
    svycoxph0 <- svycoxph(formula = Surv(time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ ePsMw, data = data0),
                          method = "breslow")
    coxph0 <- coxph(formula = Surv(time, event) ~ A,
                    data = data0,
                    method = "breslow",
                    weights = ePsMw,
                    robust = TRUE)

    ## Riskset construction with weight summaries
    risksets0 <- SiteRisksetsByStrata(time = data0[, "time"],
                                      event = data0[, "event"],
                                      A = data0[, "A"],
                                      W = data0[, "ePsMw"])
    ## IPD reconstruction from risksets with weight summaries
    risksets0$site <- 1
    risksets0$outcome <- "survival"
    risksets0$method <- "ePsMw"
    risksets0$strata <- NA
    risksets0_expanded <- RisksetsToIpdSkewed(risksets0, compress = FALSE)

    ## Reconstructed-IPD weighted analyses
    svycoxph0_expanded <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                                   design = svydesign(ids = ~ 1, weights = ~ W, data = risksets0_expanded),
                                   method = "breslow")
    coxph0_expanded <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                             data = risksets0_expanded,
                             method = "breslow",
                             weights = W,
                             robust = TRUE)

    ## Coefficients should match up exactly
    ## Two fitting methods
    expect_equal(coef(svycoxph0), coef(coxph0))
    expect_equal(coef(svycoxph0_expanded), coef(coxph0_expanded))
    ## original weighted IPD vs expaneded weighted IPD
    expect_equal(coef(svycoxph0), coef(svycoxph0_expanded))
    expect_equal(coef(coxph0), coef(coxph0_expanded))

    ## Variances
    ## Two robust variance methods
    ## Finite-population correction results in larger variance with MW, but close
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(coxph0)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^6)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^7)
    ## Variance: original weighted IPD < expaneded weighted IPD
    expect_true(as.numeric(vcov(svycoxph0)) < as.numeric(vcov(svycoxph0_expanded)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^5)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^6)
    ##
    expect_true(as.numeric(vcov(coxph0)) < as.numeric(vcov(coxph0_expanded)))
    expect_true(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^5)
    expect_false(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^6)


###  Stabilized IPTW without ties
    ## Weighted individual data (drop duplicated times)
    data0 <- as.data.frame(drnReady[[3]][,c("time","event","A","ePsSIptw")])
    data0 <- data0[!duplicated(data0$time),]
    dim(data0)

    ## IPD weighted analyses
    svycoxph0 <- svycoxph(formula = Surv(time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ ePsSIptw, data = data0),
                          method = "breslow")
    coxph0 <- coxph(formula = Surv(time, event) ~ A,
                    data = data0,
                    method = "breslow",
                    weights = ePsSIptw,
                    robust = TRUE)

    ## Riskset construction with weight summaries
    risksets0 <- SiteRisksetsByStrata(time = data0[, "time"],
                                      event = data0[, "event"],
                                      A = data0[, "A"],
                                      W = data0[, "ePsSIptw"])
    ## IPD reconstruction from risksets with weight summaries
    risksets0$site <- 1
    risksets0$outcome <- "survival"
    risksets0$method <- "ePsSIptw"
    risksets0$strata <- NA
    risksets0_expanded <- RisksetsToIpdSkewed(risksets0, compress = FALSE)

    ## Reconstructed-IPD weighted analyses
    svycoxph0_expanded <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                                   design = svydesign(ids = ~ 1, weights = ~ W, data = risksets0_expanded),
                                   method = "breslow")
    coxph0_expanded <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                             data = risksets0_expanded,
                             method = "breslow",
                             weights = W,
                             robust = TRUE)

    ## Coefficients should match up exactly
    ## Two fitting methods
    expect_equal(coef(svycoxph0), coef(coxph0))
    ## original weighted IPD vs expaneded weighted IPD
    expect_equal(coef(svycoxph0), coef(svycoxph0_expanded))
    expect_equal(coef(coxph0), coef(coxph0_expanded))

    ## Variances
    ## Two robust variance methods
    ## Finite-population correction results in larger variance with MW, but close
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(coxph0)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^4)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^5)
    ## Variance: original weighted IPD > expaneded weighted IPD
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(svycoxph0_expanded)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^3)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^4)
    ##
    expect_true(as.numeric(vcov(coxph0)) > as.numeric(vcov(coxph0_expanded)))
    expect_true(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^3)
    expect_false(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^4)


###  Stabilized IPTW with ties
    ## Weighted individual data (drop duplicated times)
    data0 <- as.data.frame(drnReady[[3]][,c("time","event","A","ePsSIptw")])

    ## IPD weighted analyses
    svycoxph0 <- svycoxph(formula = Surv(time, event) ~ A,
                          design = svydesign(ids = ~ 1, weights = ~ ePsSIptw, data = data0),
                          method = "breslow")
    coxph0 <- coxph(formula = Surv(time, event) ~ A,
                    data = data0,
                    method = "breslow",
                    weights = ePsSIptw,
                    robust = TRUE)

    ## Riskset construction with weight summaries
    risksets0 <- SiteRisksetsByStrata(time = data0[, "time"],
                                      event = data0[, "event"],
                                      A = data0[, "A"],
                                      W = data0[, "ePsSIptw"])
    ## IPD reconstruction from risksets with weight summaries
    risksets0$site <- 1
    risksets0$outcome <- "survival"
    risksets0$method <- "ePsSIptw"
    risksets0$strata <- NA
    risksets0_expanded <- RisksetsToIpdSkewed(risksets0, compress = FALSE)

    ## Reconstructed-IPD weighted analyses
    svycoxph0_expanded <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                                   design = svydesign(ids = ~ 1, weights = ~ W, data = risksets0_expanded),
                                   method = "breslow")
    coxph0_expanded <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                             data = risksets0_expanded,
                             method = "breslow",
                             weights = W,
                             robust = TRUE)

    ## Coefficients should match up exactly
    ## Two fitting methods
    expect_equal(coef(svycoxph0), coef(coxph0))
    expect_equal(coef(svycoxph0_expanded), coef(coxph0_expanded))
    ## original weighted IPD vs expaneded weighted IPD
    expect_equal(coef(svycoxph0), coef(svycoxph0_expanded))
    expect_equal(coef(coxph0), coef(coxph0_expanded))

    ## Variances
    ## Two robust variance methods
    ## Finite-population correction results in larger variance with MW, but close
    expect_true(as.numeric(vcov(svycoxph0)) > as.numeric(vcov(coxph0)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^6)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(coxph0))) < 1/10^7)
    ## Variance: original weighted IPD < expaneded weighted IPD
    expect_true(as.numeric(vcov(svycoxph0)) < as.numeric(vcov(svycoxph0_expanded)))
    expect_true(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^5)
    expect_false(abs(as.numeric(vcov(svycoxph0)) - as.numeric(vcov(svycoxph0_expanded))) < 1/10^6)
    ##
    expect_true(as.numeric(vcov(coxph0)) < as.numeric(vcov(coxph0_expanded)))
    expect_true(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^5)
    expect_false(abs(as.numeric(vcov(coxph0)) - as.numeric(vcov(coxph0_expanded))) < 1/10^6)


###  Comparison of skewed and extreme weight reconstruction
    ## Examine if robust variance is invariant to weight rearrangement within event and treatment
    risksets0_expanded_ex <- RisksetsToIpdSkewed(risksets0, compress = FALSE,
                                                 helper_fun = RisksetsToIpdExtremeHelper)
    ## NOT the same weights data
    summary(risksets0_expanded_ex$W)
    summary(risksets0_expanded$W)
    expect_false(identical(summary(risksets0_expanded_ex$W),
                           summary(risksets0_expanded$W)))
    ## Means match
    expect_equal(mean(risksets0_expanded_ex$W),
                 mean(risksets0_expanded$W))
    ## Variances match
    expect_equal(var(risksets0_expanded_ex$W),
                 var(risksets0_expanded$W))

    ## Reconstructed-IPD weighted analyses with extreme weights
    svycoxph0_expanded_ex <- svycoxph(formula = Surv(start_time, stop_time, event) ~ A,
                                      design = svydesign(ids = ~ 1, weights = ~ W, data = risksets0_expanded_ex),
                                      method = "breslow")
    coxph0_expanded_ex <- coxph(formula = Surv(start_time, stop_time, event) ~ A,
                                data = risksets0_expanded_ex,
                                method = "breslow",
                                weights = W,
                                robust = TRUE)
    ## Coefficients match
    expect_equal(coef(svycoxph0_expanded), coef(svycoxph0_expanded_ex))
    expect_equal(coef(coxph0_expanded), coef(coxph0_expanded_ex))
    ## Variances match
    expect_equal(vcov(svycoxph0_expanded), vcov(svycoxph0_expanded_ex))
    expect_equal(vcov(coxph0_expanded), vcov(coxph0_expanded_ex))

})

test_that("case-centered logistic regression is performed correctly", {

    ## Limit time points for speed
    resRisk <- subset(resRisk, eval_time < 10)

    ## Split but drop empty levels
    lstDf <- split(resRisk, f = resRisk[,c("outcome","method")],
                   drop = TRUE)


    ## Check if case-centered logistic regression is working for unweighted data
    lstDfUnwt <- lstDf[!grepl(pattern = "ePsSIptw|ePsMw", names(lstDf))]

    ## Run case-centered conditional logistic regression
    ans <- lapply(lstDfUnwt, function(df) {
        ## If unweighted, run regular case-centered approach
        ## Stratification by site is implicit by constructing risk sets at each site.
        out <- glm(formula = cbind(events_A1, events_A0) ~ 1,
                   family  = binomial(link = "logit"),
                   offset  = log(riskset_A1 / riskset_A0),
                   data    = df)
    })
    ## Expectation
    expect_equal(lapply(lstDfUnwt, CaseCenteredLogistic),
                 ans)


    ## Check depleted riskset handling
    ## If either exposed or unexposed are depleted (zero), time points should be dropped
    ## Untreated are depleted
    zeroRisksetA0 <- subset(resRisk, method == "ePsStrata" & outcome == "survival")
    zeroRisksetA0[zeroRisksetA0$eval_time > 8, c("events_A0","riskset_A0")] <- 0
    expect_equal(CaseCenteredLogistic(zeroRisksetA0),
                 CaseCenteredLogistic(subset(zeroRisksetA0, riskset_A0 > 0 & riskset_A1 > 0)))
    ## Treated are depleted
    zeroRisksetA1 <- subset(resRisk, method == "ePsStrata" & outcome == "survival")
    zeroRisksetA1[zeroRisksetA1$eval_time > 8, c("events_A1","riskset_A1")] <- 0
    CaseCenteredLogistic(zeroRisksetA1)
    expect_equal(CaseCenteredLogistic(zeroRisksetA1),
                 CaseCenteredLogistic(subset(zeroRisksetA1, riskset_A0 > 0 & riskset_A1 > 0)))


    ## Using AnalyzeSiteRisksets for unweighted
    lstBin  <- AnalyzeSiteRisksetsBin(resRisk[resRisk$outcome == "binary",])
    lstSurv <- AnalyzeSiteRisksetsSurv(resRisk[resRisk$outcome == "survival",])
    out2 <- do.call(rbind, c(lstBin, lstSurv))
    rownames(out2) <- NULL
    df2  <- data.frame(data = rep("risksets", nrow(out2)),
                       outcome = c(rep("binary", length(lstBin)),
                                   rep("survival", length(lstSurv))),
                       method = c(names(lstBin),
                                  names(lstSurv)),
                       stringsAsFactors = FALSE)
    ans2 <- cbind(df2, out2)

    ## Expectation
    expect_equal(AnalyzeSiteRisksets(resRisk),
                 ans2)


    ## Check if weighted analysis is working explicitly
    ## Efron
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedEfron(lstDf[["survival.ePsSIptw"]]),
        "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd((lstDf[["survival.ePsSIptw"]]))),
                         method = "efron"),
                "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedEfron((lstDf[["survival.ePsMw"]])), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd((lstDf[["survival.ePsMw"]]))),
                         method = "efron"), "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedEfron(lstDf[["binary.ePsSIptw"]]), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd(lstDf[["binary.ePsSIptw"]])),
                         method = "efron"), "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedEfron(lstDf[["binary.ePsMw"]]), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd(lstDf[["binary.ePsMw"]])),
                         method = "efron"), "A"))
    ## Breslow
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedBreslow(lstDf[["survival.ePsSIptw"]]), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd((lstDf[["survival.ePsSIptw"]]))),
                         method = "breslow"), "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedBreslow((lstDf[["survival.ePsMw"]])), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd((lstDf[["survival.ePsMw"]]))),
                         method = "breslow"), "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedBreslow(lstDf[["binary.ePsSIptw"]]), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd(lstDf[["binary.ePsSIptw"]])),
                         method = "breslow"), "A"))
    expect_equal(CoefVar(
        CaseCenteredLogisticWeightedBreslow(lstDf[["binary.ePsMw"]]), "A"),
        CoefVar(svycoxph(formula = Surv(start_time, stop_time, event) ~ A + strata(site),
                         design = svydesign(ids = ~ 1, weights = ~ W,
                                            data = RisksetsToIpd(lstDf[["binary.ePsMw"]])),
                         method = "breslow"), "A"))

})


###
### Tests for individual-level data regression
################################################################################

test_that("individual-level data regression is performed correctly", {

    df <- RequestSiteDataset(drnReady)

    ## Efron
    ## Binary
    ## Unadjusted stratified by site
    bin_unadjE       <- clogit(formula = Y ~ A + strata(site),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Matched cohort stratified by sites
    bin_ePsMatchE    <- clogit(formula = Y ~ A + strata(site),
                               data    = df[df$ePsMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    bin_eDrsBMatchE  <- clogit(formula = Y ~ A + strata(site),
                               data    = df[df$eDrsBMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Stratified analysis stratified by summary score strata and sites
    bin_ePsStrataE   <- clogit(formula = Y ~ A + strata(site, ePsStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    bin_eDrsBStrataE <- clogit(formula = Y ~ A + strata(site, eDrsBStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Weighted analysis stratified by sites
    bin_ePsSIptwE    <- svycoxph(formula = Surv(rep(1, length(Y)), Y) ~ A + strata(site),
                                 design = svydesign(ids = ~ 1, data = df, weights = ~ ePsSIptw),
                                 method  = c("exact", "approximate", "efron", "breslow")[3])
    bin_ePsMwE       <- svycoxph(formula = Surv(rep(1, length(Y)), Y) ~ A + strata(site),
                                 design = svydesign(ids = ~ 1, data = df, weights = ~ ePsMw),
                                 method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Survival
    ## Unadjusted stratified by site
    surv_unadjE       <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Matched cohort stratified by sites
    surv_ePsMatchE    <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df[df$ePsMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    surv_eDrsSMatchE  <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df[df$eDrsSMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Stratified analysis stratified by summary score strata and sites
    surv_ePsStrataE   <- coxph(formula = Surv(time, event) ~ A + strata(site, ePsStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    surv_eDrsSStrataE <- coxph(formula = Surv(time, event) ~ A + strata(site, eDrsSStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Weighted analysis stratified by sites
    surv_ePsSIptwE    <- svycoxph(formula = Surv(time, event) ~ A + strata(site),
                                  design = svydesign(ids = ~ 1, data = df, weights = ~ ePsSIptw),
                                  method  = c("exact", "approximate", "efron", "breslow")[3])
    surv_ePsMwE       <- svycoxph(formula = Surv(time, event) ~ A + strata(site),
                                  design = svydesign(ids = ~ 1, data = df, weights = ~ ePsMw),
                                  method  = c("exact", "approximate", "efron", "breslow")[3])

    ## Breslow
    ## Binary
    ## Unadjusted stratified by site
    bin_unadjB       <- clogit(formula = Y ~ A + strata(site),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Matched cohort stratified by sites
    bin_ePsMatchB    <- clogit(formula = Y ~ A + strata(site),
                               data    = df[df$ePsMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    bin_eDrsBMatchB  <- clogit(formula = Y ~ A + strata(site),
                               data    = df[df$eDrsBMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Stratified analysis stratified by summary score strata and sites
    bin_ePsStrataB   <- clogit(formula = Y ~ A + strata(site, ePsStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    bin_eDrsBStrataB <- clogit(formula = Y ~ A + strata(site, eDrsBStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Weighted analysis stratified by sites
    bin_ePsSIptwB    <- svycoxph(formula = Surv(rep(1, length(Y)), Y) ~ A + strata(site),
                                 design = svydesign(ids = ~ 1, data = df, weights = ~ ePsSIptw),
                                 method  = c("exact", "approximate", "efron", "breslow")[4])
    bin_ePsMwB       <- svycoxph(formula = Surv(rep(1, length(Y)), Y) ~ A + strata(site),
                                 design = svydesign(ids = ~ 1, data = df, weights = ~ ePsMw),
                                 method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Survival
    ## Unadjusted stratified by site
    surv_unadjB       <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Matched cohort stratified by sites
    surv_ePsMatchB    <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df[df$ePsMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    surv_eDrsSMatchB  <- coxph(formula = Surv(time, event) ~ A + strata(site),
                               data    = df[df$eDrsSMatch == 1,],
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Stratified analysis stratified by summary score strata and sites
    surv_ePsStrataB   <- coxph(formula = Surv(time, event) ~ A + strata(site, ePsStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    surv_eDrsSStrataB <- coxph(formula = Surv(time, event) ~ A + strata(site, eDrsSStrata),
                               data    = df,
                               method  = c("exact", "approximate", "efron", "breslow")[4])
    ## Weighted analysis stratified by sites
    surv_ePsSIptwB    <- svycoxph(formula = Surv(time, event) ~ A + strata(site),
                                  design = svydesign(ids = ~ 1, data = df, weights = ~ ePsSIptw),
                                  method  = c("exact", "approximate", "efron", "breslow")[4])
    surv_ePsMwB       <- svycoxph(formula = Surv(time, event) ~ A + strata(site),
                                  design = svydesign(ids = ~ 1, data = df, weights = ~ ePsMw),
                                  method  = c("exact", "approximate", "efron", "breslow")[4])

    ## Pile up binary (both Efron and Breslow) first
    mat <- rbind(CoefVar(bin_unadjE, "A"),
                 CoefVar(bin_ePsMatchE, "A"),
                 CoefVar(bin_eDrsBMatchE, "A"),
                 CoefVar(bin_ePsStrataE, "A"),
                 CoefVar(bin_eDrsBStrataE, "A"),
                 CoefVar(bin_ePsSIptwE, "A"),
                 CoefVar(bin_ePsMwE, "A"),
                 ## Breslow
                 CoefVar(bin_unadjB, "A"),
                 CoefVar(bin_ePsMatchB, "A"),
                 CoefVar(bin_eDrsBMatchB, "A"),
                 CoefVar(bin_ePsStrataB, "A"),
                 CoefVar(bin_eDrsBStrataB, "A"),
                 CoefVar(bin_ePsSIptwB, "A"),
                 CoefVar(bin_ePsMwB, "A"),
                 ## Survival
                 ## Efron
                 CoefVar(surv_unadjE, "A"),
                 CoefVar(surv_ePsMatchE, "A"),
                 CoefVar(surv_eDrsSMatchE, "A"),
                 CoefVar(surv_ePsStrataE, "A"),
                 CoefVar(surv_eDrsSStrataE, "A"),
                 CoefVar(surv_ePsSIptwE, "A"),
                 CoefVar(surv_ePsMwE, "A"),
                 ## Breslow
                 CoefVar(surv_unadjB, "A"),
                 CoefVar(surv_ePsMatchB, "A"),
                 CoefVar(surv_eDrsSMatchB, "A"),
                 CoefVar(surv_ePsStrataB, "A"),
                 CoefVar(surv_eDrsSStrataB, "A"),
                 CoefVar(surv_ePsSIptwB, "A"),
                 CoefVar(surv_ePsMwB, "A"))

    ans <- data.frame(data = rep("dataset", nrow(mat)),
                      outcome = rep(c("binary","survival"), each = nrow(mat) / 2),
                      ##         ## Binary (clogit by both Efron and Breslow)
                      method = c(paste0(c("unadj",
                                          "ePsMatch", "eDrsBMatch",
                                          "ePsStrata", "eDrsBStrata",
                                          "ePsSIptw", "ePsMw"),
                                        "E"),
                                 paste0(c("unadj",
                                          "ePsMatch", "eDrsBMatch",
                                          "ePsStrata", "eDrsBStrata",
                                          "ePsSIptw", "ePsMw"),
                                        "B"),
                                 ## Survival (both Efron and Breslow)
                                 paste0(c("unadj",
                                          "ePsMatch", "eDrsSMatch",
                                          "ePsStrata", "eDrsSStrata",
                                          "ePsSIptw", "ePsMw"),
                                        "E"),
                                 paste0(c("unadj",
                                          "ePsMatch", "eDrsSMatch",
                                          "ePsStrata", "eDrsSStrata",
                                          "ePsSIptw", "ePsMw"),
                                        "B")),
                      stringsAsFactors = FALSE)

    ans <- cbind(ans, mat)

    ## Expectations
    expect_equal(AnalyzeSiteDataset(df),
                 ans)
})


test_that("individual-level data COUNTERFACTUAL regression is performed correctly", {

    df <- RequestSiteTruth(drnReady)

    ## Need to set seed for randomness control
    set.seed(113)
    ## Duplicate for counterfactuals
    df_A0 <- df
    df_A1 <- df
    ## Manipulate treatment
    df_A0$A_ <- 0
    df_A1$A_ <- 1
    ## Realize binary outcome under each manipulated treatments
    ## FIXME: Stochastic and not ideal, but should work over iterations
    df_A0$Y_ <- rbinom(n = seq_along(df_A0$pY), size = 1, prob = df_A0$pY0)
    df_A1$Y_ <- rbinom(n = seq_along(df_A1$pY), size = 1, prob = df_A1$pY1)
    ## Realize survival outcome under each manipulated treatment
    df_A0$T_ <- rexp(n = seq_along(df_A0$rate0), rate = df_A0$rate0)
    df_A1$T_ <- rexp(n = seq_along(df_A1$rate1), rate = df_A1$rate1)
    ## Everybody has observed survival time
    df_A0$event_ <- 1
    df_A1$event_ <- 1
    ## Create a counterfactual dataset containing everybody twice under both treatment
    ## All manipulated or counterfactuals are marked with *_ variable names.
    df_cf <- rbind(df_A0, df_A1)

    ## Binary
    ## Efron
    bin_tAllE        <- clogit(formula = Y_ ~ A_ + strata(site),
                             data    = df_cf,
                             method  = c("exact", "approximate", "efron", "breslow")[3])
    bin_tTreatedE    <- clogit(formula = Y_ ~ A_ + strata(site),
                             data    = df_cf[df_cf$A == 1,],
                             method  = c("exact", "approximate", "efron", "breslow")[3])
    bin_tUntreatedE  <- clogit(formula = Y_ ~ A_ + strata(site),
                              data    = df_cf[df_cf$A == 0,],
                              method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Breslow
    bin_tAllB        <- clogit(formula = Y_ ~ A_ + strata(site),
                             data    = df_cf,
                             method  = c("exact", "approximate", "efron", "breslow")[4])
    bin_tTreatedB    <- clogit(formula = Y_ ~ A_ + strata(site),
                             data    = df_cf[df_cf$A == 1,],
                             method  = c("exact", "approximate", "efron", "breslow")[4])
    bin_tUntreatedB  <- clogit(formula = Y_ ~ A_ + strata(site),
                              data    = df_cf[df_cf$A == 0,],
                              method  = c("exact", "approximate", "efron", "breslow")[4])

    ## Survival
    ## Efron
    surv_tAllE       <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                             data    = df_cf,
                             method  = c("exact", "approximate", "efron", "breslow")[3])
    surv_tTreatedE   <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                             data    = df_cf[df_cf$A == 1,],
                             method  = c("exact", "approximate", "efron", "breslow")[3])
    surv_tUntreatedE <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                              data    = df_cf[df_cf$A == 0,],
                              method  = c("exact", "approximate", "efron", "breslow")[3])
    ## Breslow
    surv_tAllB       <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                             data    = df_cf,
                             method  = c("exact", "approximate", "efron", "breslow")[4])
    surv_tTreatedB   <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                             data    = df_cf[df_cf$A == 1,],
                             method  = c("exact", "approximate", "efron", "breslow")[4])
    surv_tUntreatedB <- coxph(formula = Surv(T_, event_) ~ A_ + strata(site),
                              data    = df_cf[df_cf$A == 0,],
                              method  = c("exact", "approximate", "efron", "breslow")[4])


    ## Efron
    mat <- rbind(CoefVar(bin_tAllE, "A_"),
                 CoefVar(bin_tTreatedE, "A_"),
                 CoefVar(bin_tUntreatedE, "A_"),
                 ##
                 CoefVar(bin_tAllB, "A_"),
                 CoefVar(bin_tTreatedB, "A_"),
                 CoefVar(bin_tUntreatedB, "A_"),
                 ## Survival
                 CoefVar(surv_tAllE, "A_"),
                 CoefVar(surv_tTreatedE, "A_"),
                 CoefVar(surv_tUntreatedE, "A_"),
                 ##
                 CoefVar(surv_tAllB, "A_"),
                 CoefVar(surv_tTreatedB, "A_"),
                 CoefVar(surv_tUntreatedB, "A_"))

    ans <- data.frame(data = rep("truth", 6*2),
                      outcome = rep(c("binary","survival"), each = 6),
                      ##         ## Binary (clogit by both Efron and Breslow)
                      method = rep(c(paste0(c("tAll",
                                              "tTreated",
                                              "tUntreated"),
                                            "E"),
                                     paste0(c("tAll",
                                              "tTreated",
                                              "tUntreated"),
                                            "B")), 2),
                      stringsAsFactors = FALSE)

    ans <- cbind(ans, mat)

    ## Expectations
    expect_equal(data.frame(do.call(rbind, c(AnalyzeSiteTruthBin(df_cf),
                                             AnalyzeSiteTruthSurv(df_cf))),
                            row.names = NULL),
                 ans[c("coef","var")])

    set.seed(113)
    expect_equal(AnalyzeSiteTruth(df),
                 ans)
})


###
### All analyses together
################################################################################

test_that("all analyses can be conducted together", {

    cat("### Data being used\n")
    print(drnReady)

    cat("### Time required for data preparation\n")
    cat("### RequestSiteRegression \n")
    print(system.time(reqSiteRegression <- RequestSiteRegression(drnReady)))
    cat("### RequestSiteSummary \n")
    print(system.time(reqSiteSummary    <- RequestSiteSummary(drnReady)))
    cat("### RequestSiteRisksets \n")
    print(system.time(reqSiteRisksets   <- RequestSiteRisksets(drnReady)))
    cat("### RequestSiteDataset \n")
    print(system.time(reqSiteDataset    <- RequestSiteDataset(drnReady)))
    cat("### RequestSiteTruth \n")
    print(system.time(reqSiteTruth      <- RequestSiteTruth(drnReady)))

    cat("### Time required for data analysis\n")
    cat("### AnalyzeSiteRegression \n")
    print(system.time(resSiteRegression <- AnalyzeSiteRegression(reqSiteRegression)))
    cat("### AnalyzeSiteSummary \n")
    print(system.time(resSiteSummary    <- AnalyzeSiteSummary(reqSiteSummary)))
    cat("### AnalyzeSiteRisksets \n")
    print(system.time(resSiteRisksets   <- AnalyzeSiteRisksets(reqSiteRisksets)))
    cat("### AnalyzeSiteDataset \n")
    print(system.time(resSiteDataset    <- AnalyzeSiteDataset(reqSiteDataset)))
    cat("### AnalyzeSiteTruth \n")
    ## Need seed for counterfactual realization
    set.seed(131)
    print(system.time(resSiteTruth      <- AnalyzeSiteTruth(reqSiteTruth)))

    ## All data types except truth has unadjusted analyses
    expect_true(any(grepl("unadj", resSiteRegression$method)))
    expect_true(any(grepl("unadj", resSiteSummary$method)))
    expect_true(any(grepl("unadj", resSiteRisksets$method)))
    expect_true(any(grepl("unadj", resSiteDataset$method)))
    ## Breslow and Efron are all paired
    expect_equal(sum(grepl("unadjE", resSiteRegression$method)),
                 sum(grepl("unadjB", resSiteRegression$method)))
    expect_equal(sum(grepl("unadjE", resSiteSummary$method)),
                 sum(grepl("unadjB", resSiteSummary$method)))
    expect_equal(sum(grepl("unadjE", resSiteRisksets$method)),
                 sum(grepl("unadjB", resSiteRisksets$method)))
    expect_equal(sum(grepl("unadjE", resSiteDataset$method)),
                 sum(grepl("unadjB", resSiteDataset$method)))

    ## Combine as an answer to the wrapper function
    ans <- dplyr::bind_rows(resSiteRegression,
                            resSiteSummary,
                            resSiteRisksets,
                            resSiteDataset,
                            resSiteTruth)

    ## Data types
    data_types <- c("meta", "summary", "risksets", "dataset", "truth")
    expect_true(all(ans$data %in% data_types))
    expect_true(all(data_types %in% ans$data))
    ans$data    <- factor(ans$data,
                          levels = data_types)
    ## Outcome types
    outcome_types <- c("binary", "survival")
    expect_true(all(ans$outcome %in% outcome_types))
    expect_true(all(outcome_types %in% ans$outcome))
    ans$outcome <- factor(ans$outcome,
                          levels = outcome_types)
    ## Do not specify levels for method as they are subject to change
    ans$method  <- factor(ans$method)

    ## Reorder for readability
    ans <- dplyr::arrange_(ans, "data", "outcome", "method")

    ## Expectations
    if (length(suppressWarnings(system("type sas", intern = TRUE, ignore.stderr = TRUE))) == 0) {
        ## If sas is not available, drop ones with expected NA's
        expect_true(all(!is.na(ans[!(ans$outcome == "survival" & ans$method %in% c("ePsMwB","ePsMwE","ePsSIptwB","ePsSIptwE") & ans$data == "risksets"),
                                   c("outcome","method","data","coef","var")])))
    } else {
        ## If sas is available
        expect_true(all(!is.na(ans[c("outcome","method","data","coef","var")])))
    }


    ## Conduct all analyses to gether
    cat("### Time required for entire analysis together\n")
    ## Need seed for counterfactual realization
    set.seed(131)
    print(system.time(res <- Analyze(drnReady)))

    ## These columns should not have NA's
    if (length(suppressWarnings(system("type sas", intern = TRUE, ignore.stderr = TRUE))) == 0) {
        ## If sas is not available, drop ones with expected NA's
        expect_true(all(!is.na(res[!(res$outcome == "survival" & res$method %in% c("ePsMwB","ePsMwE","ePsSIptwB","ePsSIptwE") & res$data == "risksets"),
                                   c("outcome","method","data","coef","var")])))
    } else {
        ## If sas is available
        expect_true(all(!is.na(res[c("outcome","method","data","coef","var")])))
    }

    ## Site result columns included
    expect_true(any(grepl("coef_site", names(res))))
    expect_true(any(grepl("var_site", names(res))))
    ## All columns should not have NA's for meta-analysis (site coef/var included).
    expect_true(all(!is.na(subset(res, method == "meta"))))
    ## Match up
    expect_equal(res,
                 ans)

})
