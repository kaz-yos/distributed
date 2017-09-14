################################################################################
### Test data preparation functions
##
## Created on: 2016-05-17
## Author: Kazuki Yoshida
################################################################################

library(testthat)

### Context (1 for each file)
context("### Test 02.ParepareDataWithinSites.R")


set.seed(20160505)


### Create a very simple example data

## Load real data example for PS estimation
library(MatchIt)
data(lalonde)
## Add PS
lalonde_formula <- treat ~ age + educ + black + hispan + married + nodegree
glm1 <- glm(formula = lalonde_formula,
            family  = binomial(link = "logit"),
            data    = lalonde)
lalonde$PS1 <- predict(glm1, type = "response")


## Sample size
n        <- 10^4
## Randomize treatment
alpha0   <- 0
alphaX   <- rep(0,7)
## Randomize binary outcome
beta0    <- 0
betaX    <- rep(0,7)
betaA    <- 0
betaXA   <- rep(0,7)
## Mean event time of 1
lambda   <- 1
lambda_c <- 0.5
Tmax     <- 1
## Generate dataset
df1 <- GenerateOneCenter(n = n,
                         alphas = c(alpha0, alphaX),
                         betas = c(beta0, betaX, betaA, betaXA),
                         survParams = c(lambda, lambda_c, Tmax))
## Prepare within site with helper variables
df_prep1 <- SitePrepareHelperVariables(df1)


## Create a pathological data with NA adjustment helper variables
df_prep1_NA <- df_prep1
df_prep1_NA[,c("ePsSIptw", "ePsMw",
               "ePsStrata", "eDrsBStrata", "eDrsSStrata")] <- NA
## when matching is impossible, match indicator is all zeros
df_prep1_NA[,c("ePsMatch", "eDrsBMatch", "eDrsSMatch")] <- 0


## Generate data where model fitting fail
df_rare_tx <- GenerateOneCenter(n = n,
                                alphas = c(-10, alphaX),
                                betas = c(beta0, betaX, betaA, betaXA),
                                survParams = c(lambda, lambda_c, Tmax))
summary(df_rare_tx)
## Only keep the untreated
df_rare_tx <- subset(df_rare_tx, A == 0)
## df_prep_rare_tx <- SitePrepareHelperVariables(df_rare_tx)


## Rare outcome
df_rare_dis <- GenerateOneCenter(n = 2*10^4,
                                 alphas = c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                                 betas = c(beta0 = -10.0,
                                           betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                           betaA = 0,
                                           betaXA = rep(0,7)),
                                 survParams = c(-log(0.99995), -log(0.9), 1))
## Drop events in the untreated
df_rare_dis <- subset(df_rare_dis, !(A == 0 & Y == 1))
df_rare_dis <- subset(df_rare_dis, !(A == 0 & event == 1))
summary(df_rare_dis)
## df_prep_rare_dis <- SitePrepareHelperVariables(df_rare_dis)


test_that("estimated weights are reasonable", {

    ## Max 1 for matching weight
    expect_equal(max(df_prep1$ePsMw), 1)
    ## Mean stabilized IPTW should be around 1
    expect_equal(round(mean(df_prep1$ePsSIptw), 5), 1)

})


###
### Functions for helper variable generation
################################################################################

test_that("PS are generated correctly", {

    ## Compare to lalonde PS generated manually
    expect_equal(EstimateScoreBin(formula = lalonde_formula, data = lalonde, newdata = lalonde),
                 lalonde$PS1)


    ## Pathological data with only one treatment arm
    ## Warning due to model non-convergence
    expect_warning(EstimateScoreBin(formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                                     data = df_rare_tx, newdata = df_rare_tx))
    ## When there are only untreated, all PS are very small
    expect_true(all(suppressWarnings(
        EstimateScoreBin(formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                         data = df_rare_tx, newdata = df_rare_tx)) < 10^(-10)))
    ## Model fitting error (induced by non-exisiting covariate) gives NA's
    expect_equal(EstimateScoreBin(formula = A ~ X1 + X2 + X3 + X4 + X5 + X6 + X8,
                                  data = df_rare_tx, newdata = df_rare_tx),
                 rep(as.numeric(NA), nrow(df_rare_tx)))
})


test_that("DRS for binary outcome are generated correctly", {

    ## Only use the untreated
    glm1 <- glm(formula = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                family  = binomial(link = "logit"),
                data    = subset(df1, A == 0))
    eDrsStd <- as.numeric(predict(glm1, newdata = df1, type = "response"))

    eDrs <- EstimateScoreBin(formula = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                             data = subset(df1, A == 0),
                             newdata = df1)
    ## Compare to the standard
    expect_equal(eDrs, eDrsStd)
    ## Must be [0,1]
    expect_true(all(eDrs >= 0))
    expect_true(all(eDrs <= 1))


    ## Pathological data with no untreated events
    ## Warning due to model non-convergence
    expect_warning(EstimateScoreBin(formula = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                                    data = subset(df_rare_dis, A == 0),
                                    newdata = df_rare_dis))
    ## When there are no untreated events, all DRS are very small
    expect_true(all(suppressWarnings(
        EstimateScoreBin(formula = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                         data = subset(df_rare_dis, A == 0),
                         newdata = df_rare_dis)) < 10^(-10)))
    ## Model fitting error (induced by non-exisiting covariate) gives NA's
    expect_equal(EstimateScoreBin(formula = Y ~ X1 + X2 + X3 + X4 + X5 + X6 + X8,
                                  data = subset(df_rare_dis, A == 0),
                                  newdata = df_rare_dis),
                 rep(as.numeric(NA), nrow(df_rare_dis)))
})


test_that("DRS for survival outcome are generated correctly", {

    ## Only use the untreated
    coxph1 <- coxph(formula = Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                    data    = subset(df1, A == 0))

    ## predict.coxph() details
    ## https://stat.ethz.ch/R-manual/R-devel/library/survival/html/predict.coxph.html

    ## How to interpret the output of predict.coxph?
    ## http://stats.stackexchange.com/questions/44896/how-to-interpret-the-output-of-predict-coxph

    ## Predict probability from Cox PH model
    ## http://stackoverflow.com/questions/30993197/predict-probability-from-cox-ph-model

    cat("###  predict.coxph\n")
    cbind(lp = predict(coxph1, newdata = df1, type = "lp"),
          risk = predict(coxph1, newdata = df1, type = "risk"),
          expected = predict(coxph1, newdata = df1, type = "expected"),
          predict(coxph1, newdata = df1, type = "terms")) %>%
        head(., n = 10) %>% print

    eDrsStd <- as.numeric(predict(coxph1, newdata = df1, type = "lp"))

    eDrs <- EstimateScoreSurv(formula = Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                              data = subset(df1, A == 0),
                              newdata = df1)
    ## Compare to the standard
    expect_equal(eDrs, eDrsStd)
    ## [-Inf, Inf]. No good way to test?


    ## Pathological data with no untreated events
    ## Warning due to model fit failure
    expect_warning(EstimateScoreSurv(formula = Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                                     data = subset(df_rare_dis, A == 0),
                                     newdata = df_rare_dis))
    ## Model fitting error. coxph gives no model fit here.
    expect_equal(suppressWarnings(
        EstimateScoreSurv(formula = Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7,
                          data = subset(df_rare_dis, A == 0),
                          newdata = df_rare_dis)),
        rep(as.numeric(NA), nrow(df_rare_dis)))
    ## Model fitting error (induced by non-exisiting covariate) gives NA's
    expect_equal(EstimateScoreSurv(formula = Surv(time, event) ~ X1 + X2 + X3 + X4 + X5 + X6 + X8,
                                   data = subset(df_rare_dis, A == 0),
                                   newdata = df_rare_dis),
                 rep(as.numeric(NA), nrow(df_rare_dis)))
})


test_that("IPW are generated correctly", {

    ## Write definitions
    df1 <- data.frame(status = c(0,     0,     1,     1,     1),
                      score  = c(0.2,   0.3,   0.4,   0.5,   0.6),
                      iptw   = c(1/0.8, 1/0.7, 1/0.4, 1/0.5, 1/0.6),
                      siptw  = c(1/0.8 * 2/5,
                                 1/0.7 * 2/5,
                                 1/0.4 * 3/5,
                                 1/0.5 * 3/5,
                                 1/0.6 * 3/5),
                      mw     = c(1/0.8 * 0.2,
                                 1/0.7 * 0.3,
                                 1/0.4 * 0.4,
                                 1/0.5 * 0.5,
                                 1/0.6 * 0.4))

    expect_equal(CreateIpw(score = df1$score, status = df1$status, type = "naive"), df1$iptw)
    expect_equal(CreateIpw(score = df1$score, status = df1$status, type = "stabilized"), df1$siptw)
    expect_equal(CreateIpw(score = df1$score, status = df1$status, type = "mw"), df1$mw)


    ## NA PS gives NA weights
    expect_equal(CreateIpw(score = as.numeric(c(NA,NA)), status = c(0,1), type = "naive"),
                 as.numeric(c(NA,NA)))
    expect_equal(CreateIpw(score = as.numeric(c(NA,NA)), status = c(0,1), type = "stabilized"),
                 as.numeric(c(NA,NA)))
    expect_equal(CreateIpw(score = as.numeric(c(NA,NA)), status = c(0,1), type = "mw"),
                 as.numeric(c(NA,NA)))
})


test_that("Strata are generated correctly", {

    ## 4 categories
    score <- seq(0,1, length.out = 1000)
    quartiles <- quantile(score)
    scoreCat4 <- as.numeric(cut(x = score,
                               breaks = quartiles,
                               right = FALSE,
                               include.lowest = TRUE))
    expect_equal(CreateStrata(score = score, nStrata = 4), scoreCat4)
    expect_true(all(table(CreateStrata(score = score, nStrata = 4)) == 250))
    expect_true(length(table(CreateStrata(score = score, nStrata = 4))) == 4)


    ## 10 categories
    deciles <- quantile(score, probs = seq(0, 1, length.out = 10 + 1))
    scoreCat10 <- as.numeric(cut(x = score,
                                 breaks = deciles,
                                 right = FALSE,
                                 include.lowest = TRUE))
    expect_equal(CreateStrata(score = score, nStrata = 10), scoreCat10)
    expect_true(all(table(CreateStrata(score = score, nStrata = 10)) == 100))
    expect_true(length(table(CreateStrata(score = score, nStrata = 10))) == 10)


    ## If only < nStrata discrete values, just utilize discreteness as is
    score_level1 <- rep(1,10)
    score_level2 <- c(1,1,1, rep(5, 7))
    score_level3 <- c(1,1,1, 3, 3, rep(8, 5))
    expect_equal(CreateStrata(score = score_level1, nStrata = 10),
                 rep(1, 10))
    expect_equal(CreateStrata(score = score_level2, nStrata = 10),
                 c(1,1,1, rep(2, 7)))
    expect_equal(CreateStrata(score = score_level3, nStrata = 10),
                 c(1,1,1, 2, 2, rep(3, 5)))


    ## If nStrata = 1, just assign 1 to everybody
    expect_equal(CreateStrata(1:7, nStrata = 1),
                 rep(1,7))


    ## Avoid threshold value duplication by unique()
    score_few_values <- c(rep(0.1, 4), rep(0.2, 6))
    uniq_cutoffs     <- unique(quantile(score_few_values, probs = seq(0, 1, length.out = 10 + 1)))
    ansScoreCatFew   <- as.numeric(cut(x = score_few_values,
                                       breaks = uniq_cutoffs,
                                       right = FALSE,
                                       include.lowest = TRUE))
    expect_equal(CreateStrata(score_few_values, nStrata = 10),
                 ansScoreCatFew)


    ## NA scores gives NA strata indicators
    expect_equal(CreateStrata(score = rep(as.numeric(NA), 1000), nStrata = 10),
                 rep(as.numeric(NA), 1000))
})


test_that("Matched indicators are generated correctly", {

    ## Using lalonde dataset
    lalonde$PS2 <- EstimateScoreBin(formula = lalonde_formula,
                                    data = lalonde, newdata = lalonde)

    lalonde$matched <- CreateMatches(score = lalonde$PS2, status = lalonde$treat,
                                     logit = TRUE, caliper = 0.05)

    matched_count <- tapply(lalonde$matched, lalonde$treat, FUN = sum)
    expect_equal(length(matched_count), 2)
    expect_equal(as.numeric(matched_count[1]),
                 as.numeric(matched_count[2]))


    ## Simpler data
    score <- c(0.1, 0.11, 0.2, 0.21, 0.3, 0.8)
    status <- rep(c(0,1), 3)

    ## With a tight caliper only first four match
    expect_equal(CreateMatches(score = score, status = status,
                               logit = FALSE, caliper = 0.05),
                 c(1,1,1,1,0,0))

    ## With a very loose caliper all can match (3 * sd(score) > abs(0.3 - 0.8))
    expect_equal(CreateMatches(score = score, status = status,
                               logit = FALSE, caliper = 3),
                 c(1,1,1,1,1,1))


    ## When no matches can occur by data anomaly, the match indicator is all zeros.
    ## NA summary scores gives all-zero match indicators
    expect_equal(CreateMatches(score = rep(as.numeric(NA), 1000), status = rep(c(0,1), 500),
                               logit = FALSE, caliper = 0.05),
                 rep(0, 1000))
    ## Having only one treatment group also gives all-zero match indicators
    expect_equal(CreateMatches(score = rep(as.numeric(0.5), 1000), status = rep(1, 1000),
                               logit = FALSE, caliper = 0.05),
                 rep(0, 1000))
})


test_that("dataset is prepared correctly", {

    expect_true("ResSiteReady" %in% class(df_prep1))

    ## Variables to expect
    vars <- c("ePs",
              "eDrsB",
              "eDrsS",
              ## Weights
              "ePsSIptw",
              "ePsMw",
              ## Strata indicators
              "ePsStrata",
              "eDrsBStrata",
              "eDrsSStrata",
              ## Matching indicators
              "ePsMatch",
              "eDrsBMatch",
              "eDrsSMatch")

    ## Includes correct variables
    expect_true(all(vars %in% names(df_prep1)))

    ## 8 Covariates
    df2 <- GenerateOneCenter(n = n,
                             alphas = c(alpha0, c(alphaX,0)),
                             betas = c(beta0, c(betaX,0), betaA, c(betaXA,0)),
                             survParams = c(lambda, lambda_c, Tmax))
    ## Correct number of covariates
    expect_equal(length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(df2))), 8)
    ## Prepare
    df_prep2 <- SitePrepareHelperVariables(df2)
    ## Includes correct variables
    expect_true(all(vars %in% names(df_prep2)))

    ## 6 Covariates
    df3 <- GenerateOneCenter(n = n,
                             alphas = c(alpha0, alphaX[-1]),
                             betas = c(beta0, betaX[-1], betaA, betaXA[-1]),
                             survParams = c(lambda, lambda_c, Tmax))
    ## Correct number of covariates
    expect_equal(length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(df3))), 6)
    ## Prepare
    df_prep3 <- SitePrepareHelperVariables(df3)
    ## Includes correct variables
    expect_true(all(vars %in% names(df_prep3)))

})


###
### Functions for within-site regression
################################################################################

test_that("site-specific logistic regression is conducted correctly", {

    ## Helper function
    lstModels1 <- FitSiteLogistic(df_prep1)
    ## Check first class elements
    expect_equal(as.character(sapply(lstModels1,
                                     function(model) {
                                         class(model)[1]
                                     })),
                 c("glm",
                   "glm","glm",
                   "clogit","clogit",
                   "svyglm","svyglm",
                   ## Breslow
                   "clogit","clogit"))
    ## Fitting methods
    expect_equal(as.character(sapply(lstModels1, "[", "method")),
                 c("glm.fit",
                   "glm.fit","glm.fit",
                   "efron","efron",
                   "glm.fit","glm.fit",
                   "breslow","breslow"))

    ## Variance examination

    ## IPTW using glm
    ## quasibinomial
    logit_ePsSIptwQuasi <- glm(formula = Y ~ A,
                               family  = quasibinomial(link = "logit"),
                               weights = ePsSIptw,
                               data    = df_prep1)
    ## binomial
    suppressWarnings({
        logit_ePsSIptwBinom <- glm(formula = Y ~ A,
                                   family  = binomial(link = "logit"),
                                   weights = ePsSIptw,
                                   data    = df_prep1)
    })
    ## Their model-based vcov are NOT identical.
    expect_true(!identical(vcov(logit_ePsSIptwQuasi),
                           vcov(logit_ePsSIptwBinom)))
    ## Their robust vcov are nearly identical as the model-based ones are irrelevant here.
    expect_true(identical(round(sandwich::sandwich(logit_ePsSIptwQuasi), 13),
                          round(sandwich::sandwich(logit_ePsSIptwBinom), 13)))

    ## IPTW using svyglm
    suppressWarnings({
        logit_ePsSIptwBinomSvy <- svyglm(formula = Y ~ A,
                                         family  = binomial(link = "logit"),
                                         design  = svydesign(ids = ~ 1,
                                                             weights = ~ ePsSIptw,
                                                             data = df_prep1))
    })
    logit_ePsSIptwQuasiSvy <- svyglm(formula = Y ~ A,
                                     family  = quasibinomial(link = "logit"),
                                     design  = svydesign(ids = ~ 1,
                                                         weights = ~ ePsSIptw,
                                                         data = df_prep1))
    ## Their vcov are the same since they are already robust vcov.
    ## Do not apply sandwich on these.
    expect_true(identical(vcov(logit_ePsSIptwBinomSvy),
                          vcov(logit_ePsSIptwQuasiSvy)))
    ##  The glm-sandwich vcov and svyglm vcov should be very close.
    expect_true(all(abs((sandwich::sandwich(logit_ePsSIptwQuasi) - vcov(logit_ePsSIptwQuasiSvy)) /
                        vcov(logit_ePsSIptwQuasiSvy)) < 0.001))


    ## MW using glm
    ## quasibinomial
    logit_ePsMwQuasi <- glm(formula = Y ~ A,
                            family  = quasibinomial(link = "logit"),
                            weights = ePsMw,
                            data    = df_prep1)
    ## binomial
    suppressWarnings({
        logit_ePsMwBinom <- glm(formula = Y ~ A,
                                family  = binomial(link = "logit"),
                                weights = ePsMw,
                                data    = df_prep1)
    })
    ## Their model-based vcov are NOT identical.
    expect_true(!identical(vcov(logit_ePsMwQuasi),
                           vcov(logit_ePsMwBinom)))
    ## Their robust vcov are nearly identical as the model-based ones are irrelevant here.
    expect_true(identical(round(sandwich::sandwich(logit_ePsMwQuasi), 10),
                          round(sandwich::sandwich(logit_ePsMwBinom), 10)))

    ## MW using svyglm
    suppressWarnings({
        logit_ePsMwBinomSvy <- svyglm(formula = Y ~ A,
                                      family  = binomial(link = "logit"),
                                      design  = svydesign(ids = ~ 1,
                                                          weights = ~ ePsMw,
                                                          data = df_prep1))
    })
    logit_ePsMwQuasiSvy <- svyglm(formula = Y ~ A,
                                  family  = quasibinomial(link = "logit"),
                                  design  = svydesign(ids = ~ 1,
                                                      weights = ~ ePsMw,
                                                      data = df_prep1))
    ## Their vcov are the same since they are already robust vcov.
    ## Do not apply sandwich on these.
    expect_true(identical(vcov(logit_ePsMwBinomSvy),
                          vcov(logit_ePsMwQuasiSvy)))
    ## The glm-sandwich vcov and svyglm vcov should be very close.
    expect_true(all(abs((sandwich::sandwich(logit_ePsMwQuasi) - vcov(logit_ePsMwQuasiSvy)) /
                        vcov(logit_ePsMwQuasiSvy)) < 0.001))


###  Keep only coefficients and variance

    ## Extract coefficients for A
    coefsA <- sapply(lstModels1, function(model) {
        as.numeric(coef(model)["A"])
    })
    ## Extract variances for A
    varsA <- sapply(lstModels1, function(model) {
        as.numeric(vcov(model)["A","A"])
    })
    ## Collect in a data frame
    ans1 <- data.frame(method = names(coefsA),
                       coef = coefsA,
                       var = varsA,
                       stringsAsFactors = FALSE)
    rownames(ans1) <- NULL

    ## Logistic regression
    res1 <- SiteRegressionBin(df_prep1)
    ## Unadjusted x 1, Matched x 2, (Stratified x 2) x 2, Weighted x 2
    expect_equal(nrow(res1), 9)
    ## method, coef, var columns
    expect_equal(ncol(res1), 3)
    ## Check first elements match the coefs
    expect_equal(as.numeric(res1[,"coef"]),
                 as.numeric(coefsA))
    ## Check second elements match the vars
    expect_equal(as.numeric(res1[,"var"]),
                 as.numeric(varsA))
    ## Check weighted analyses specifically
    expect_equal(as.numeric(res1[res1$method %in% c("ePsSIptw","ePsMw"),"var"]),
                 c(as.numeric(vcov(logit_ePsSIptwQuasiSvy)["A","A"]),
                   as.numeric(vcov(logit_ePsMwQuasiSvy)["A","A"])))
    ## Check equality of df
    expect_equal(res1, ans1)

})


test_that("site-specific Cox regression is conducted correctly", {

    ## Helper function
    lstModels1 <- FitSiteCox(df_prep1)
    ## Check first class elements
    expect_equal(as.character(sapply(lstModels1,
                                     function(model) {
                                         class(model)[1]
                                     })),
                 c("coxph",
                   "coxph","coxph",
                   "coxph","coxph",
                   "svycoxph","svycoxph",
                   ## Breslow
                   "coxph",
                   "coxph","coxph",
                   "coxph","coxph",
                   "svycoxph","svycoxph"))
    ## Check approximation methods
    expect_equal(as.character(sapply(lstModels1, "[", "method")),
                 rep(c("efron", "breslow"),
                     each = 7))

    ## IPTW
    ## coxph with robust = TRUE (vcov gives robust vcov)
    cox_ePsSIptw <- coxph(formula = Surv(time, event) ~ A,
                          weights = ePsSIptw,
                          data    = df_prep1,
                          robust  = TRUE)
    cox_ePsSIptwSvy <- svycoxph(formula = Surv(time, event) ~ A,
                                design  = svydesign(ids = ~ 1,
                                                    weights = ~ ePsSIptw,
                                                    data = df_prep1))
    ## coxph-robust vcov and svycoxph vcov should be very close.
    expect_true(all(abs((vcov(cox_ePsSIptw) - vcov(cox_ePsSIptwSvy)) /
                        vcov(cox_ePsSIptwSvy)) < 0.001))

    ## MW
    ## coxph with robust = TRUE (vcov gives robust vcov)
    cox_ePsMw <- coxph(formula = Surv(time, event) ~ A,
                       weights = ePsMw,
                       data    = df_prep1,
                       robust  = TRUE)
    cox_ePsMwSvy <- svycoxph(formula = Surv(time, event) ~ A,
                             design  = svydesign(ids = ~ 1,
                                                 weights = ~ ePsMw,
                                                 data = df_prep1))
    ## coxph-robust vcov and svycoxph vcov should be very close.
    expect_true(all(abs((vcov(cox_ePsMw) - vcov(cox_ePsMwSvy)) /
                        vcov(cox_ePsMwSvy)) < 0.001))


###  Keep only coefficients and variance
    ## Extract coefficients for A
    coefsA <- sapply(lstModels1, function(model) {
        as.numeric(coef(model)["A"])
    })
    ## Extract variances for A
    varsA <- sapply(lstModels1, function(model) {
        as.numeric(vcov(model)["A","A"])
    })
    ## Collect in a data frame
    ans1 <- data.frame(method = names(coefsA),
                       coef = coefsA,
                       var = varsA,
                       stringsAsFactors = FALSE)
    rownames(ans1) <- NULL

    ## Cox PH
    res1 <- SiteRegressionSurv(df_prep1)
    ## (Unadjusted x 1, Matched x 2, Stratified x 2, Weighted x 2) x 2
    expect_equal(nrow(res1), 14)
    ## method, coef, var columns
    expect_equal(ncol(res1), 3)
    ## Check first elements match the coefs
    expect_equal(as.numeric(res1[,"coef"]),
                 as.numeric(coefsA))
    ## Check second elements match the vars
    expect_equal(as.numeric(res1[,"var"]),
                 as.numeric(varsA))
    ## Check weighted analyses specifically
    expect_equal(as.numeric(res1[res1$method %in% c("ePsSIptwE","ePsMwE"),"var"]),
                 c(as.numeric(vcov(cox_ePsSIptwSvy)["A","A"]),
                   as.numeric(vcov(cox_ePsMwSvy)["A","A"])))
    ## Check equality of df
    expect_equal(res1, ans1)

})


test_that("site-specific regressions (binary and survival) are conducted correctly (also NA helper variablehandling)", {

    ## Expected data frame
    binary   <- SiteRegressionBin(df_prep1)
    survival <- SiteRegressionSurv(df_prep1)
    binary$outcome   <- "binary"
    survival$outcome <- "survival"
    ans <- dplyr::bind_rows(binary,
                            survival)[c("outcome","method","coef","var")]

    expect_equal(SiteRegression(df_prep1),
                 ans)


    ## When helper variables are NA,
    ## all values should be NA except for unadjusted analysis results
    ans_NA <- ans
    ans_NA[!grepl("unadj", ans_NA$method), c("coef","var")] <- NA

    expect_equal(suppressWarnings(SiteRegression(df_prep1_NA)),
                 ans_NA)
})


###
### Functions for within-site summary dataset construction
################################################################################

test_that("site-specific summary data construction works for binary outcome", {
    ## Unadjusted
    ans_unadj <- df_prep1 %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(Y),
                         denom  = length(Y)) %>%
        as.data.frame
    ## Test
    res_unadj <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                         data = df_prep1)
    expect_equal(res_unadj,
                 ans_unadj)

    ## PS matching
    ans_ePsMatch <- subset(df_prep1, ePsMatch == 1) %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(Y),
                         denom  = length(Y)) %>%
        as.data.frame
    ## Test
    res_ePsMatch <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                            data = subset(df_prep1, ePsMatch == 1))
    expect_equal(res_ePsMatch,
                 ans_ePsMatch)
    ## Matching is 1:1
    expect_equal(as.numeric(diff(res_ePsMatch$denom)), 0)


    ## Binary DRS matching
    ans_eDrsBMatch <- subset(df_prep1, eDrsBMatch == 1) %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(Y),
                         denom  = length(Y)) %>%
        as.data.frame
    ## Test
    res_eDrsBMatch <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                              data = subset(df_prep1, eDrsBMatch == 1))
    expect_equal(res_eDrsBMatch,
                 ans_eDrsBMatch)
    ## Matching is 1:1
    expect_equal(as.numeric(diff(res_eDrsBMatch$denom)), 0)


    ## Handling empty matched data set
    res_NoMatch <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                           data = subset(df_prep1_NA, eDrsBMatch %in% 1))
    expect_equal(res_NoMatch,
                 dplyr::tibble(A      = as.integer(NA),
                               events = as.integer(NA),
                               denom  = as.integer(NA)))


    ## Stratum-specific by PS
    ans_ePsStrata <- df_prep1 %>%
        dplyr::group_by(ePsStrata, A) %>%
        dplyr::summarize(events = sum(Y),
                         denom  = length(Y)) %>%
        as.data.frame
    names(ans_ePsStrata)[1] <- "strata"
    ## Test
    res_ePsStrata <- SummarizeStratumEventsTotalsBy(var = "Y", by = "A", strata = "ePsStrata",
                                                    data = df_prep1)
    expect_equal(res_ePsStrata,
                 ans_ePsStrata)

    ## Stratum-specific by binary DRS
    ans_eDrsBStrata <- df_prep1 %>%
        dplyr::group_by(eDrsBStrata, A) %>%
        dplyr::summarize(events = sum(Y),
                         denom  = length(Y)) %>%
        as.data.frame
    names(ans_eDrsBStrata)[1] <- "strata"
    ## Test
    res_eDrsBStrata <- SummarizeStratumEventsTotalsBy(var = "Y", by = "A", strata = "eDrsBStrata",
                                                      data = df_prep1)
    expect_equal(res_eDrsBStrata,
                 ans_eDrsBStrata)


    ## Handling empty stratified data
    res_NoStrata <- SummarizeStratumEventsTotalsBy(var = "Y", by = "A", strata = "eDrsBStrata",
                                                   data = subset(df_prep1_NA, !is.na(eDrsBStrata)))
    expect_equal(res_NoStrata,
                 dplyr::tibble(strata = as.integer(NA),
                               A      = as.integer(NA),
                               events = as.integer(NA),
                               denom  = as.integer(NA)))


    ## Combination
    ans_unadj$strata       <- NA
    ans_ePsMatch$strata    <- NA
    ans_eDrsBMatch$strata  <- NA
    ans_unadj$method       <- "unadj"
    ans_ePsMatch$method    <- "ePsMatch"
    ans_eDrsBMatch$method  <- "eDrsBMatch"
    ans_ePsStrata$method   <- "ePsStrata"
    ans_eDrsBStrata$method <- "eDrsBStrata"
    ans_combo <- rbind(ans_unadj,
                       ans_ePsMatch,
                       ans_eDrsBMatch,
                       ans_ePsStrata,
                       ans_eDrsBStrata)[c("method","strata","A","events","denom")]
    ## Check
    expect_equal(SiteSummaryBin(df_prep1),
                 ans_combo)

})


test_that("site-specific summary data construction works for survival outcome", {
    ## Unadjusted
    ans_unadj <- df_prep1 %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(event),
                         denom  = sum(time)) %>%
        as.data.frame
    ## Test
    ## Here PT do not necessarily match
    res_unadj <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                        data = df_prep1)
    expect_equal(res_unadj,
                 ans_unadj)

    ## PS matching
    ans_ePsMatch <- subset(df_prep1, ePsMatch == 1) %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(event),
                         denom  = sum(time)) %>%
        as.data.frame
    ## Test
    ## Here PT do not necessarily match
    res_ePsMatch <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                        data = subset(df_prep1, ePsMatch == 1))
    expect_equal(res_ePsMatch,
                 ans_ePsMatch)


    ## Binary DRS matching
    ans_eDrsSMatch <- subset(df_prep1, eDrsSMatch == 1) %>%
        dplyr::group_by(A) %>%
        dplyr::summarize(events = sum(event),
                         denom  = sum(time)) %>%
        as.data.frame
    ## Test
    ## Here PT do not necessarily match
    res_eDrsSMatch <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                          data = subset(df_prep1, eDrsSMatch == 1))
    expect_equal(res_eDrsSMatch,
                 ans_eDrsSMatch)


    ## Handling empty matched data set
    res_NoMatch <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                       data = subset(df_prep1_NA, eDrsBMatch %in% 1))
    expect_equal(res_NoMatch,
                 dplyr::tibble(A      = as.integer(NA),
                               events = as.integer(NA),
                               denom  = as.integer(NA)))



    ## Stratum-specific by PS
    ans_ePsStrata <- df_prep1 %>%
        dplyr::group_by(ePsStrata, A) %>%
        dplyr::summarize(events = sum(event),
                         denom  = sum(time)) %>%
        as.data.frame
    names(ans_ePsStrata)[1] <- "strata"
    ## Test
    res_ePsStrata <- SummarizeStratumEventsPtBy(event = "event", time = "time",
                                                by = "A", strata = "ePsStrata",
                                                data = df_prep1)
    expect_equal(res_ePsStrata,
                 ans_ePsStrata)


    ## Stratum-specific by survival DRS
    ans_eDrsSStrata <- df_prep1 %>%
        dplyr::group_by(eDrsSStrata, A) %>%
        dplyr::summarize(events = sum(event),
                         denom  = sum(time)) %>%
        as.data.frame
    names(ans_eDrsSStrata)[1] <- "strata"
    ## Test
    res_eDrsSStrata <- SummarizeStratumEventsPtBy(event = "event", time = "time",
                                                  by = "A", strata = "eDrsSStrata",
                                                  data = df_prep1)
    expect_equal(res_eDrsSStrata,
                 ans_eDrsSStrata)


    ## Handling empty stratified data
    res_NoStrata <- SummarizeStratumEventsPtBy(event = "event", time = "time",
                                               by = "A", strata = "eDrsSStrata",
                                               data = subset(df_prep1_NA, !is.na(eDrsBStrata)))
    expect_equal(res_NoStrata,
                 dplyr::tibble(strata = as.integer(NA),
                               A      = as.integer(NA),
                               events = as.integer(NA),
                               denom  = as.integer(NA)))


    ## Combination
    ans_unadj$strata       <- NA
    ans_ePsMatch$strata    <- NA
    ans_eDrsSMatch$strata  <- NA
    ans_unadj$method       <- "unadj"
    ans_ePsMatch$method    <- "ePsMatch"
    ans_eDrsSMatch$method  <- "eDrsSMatch"
    ans_ePsStrata$method   <- "ePsStrata"
    ans_eDrsSStrata$method <- "eDrsSStrata"
    ans_combo <- rbind(ans_unadj,
                       ans_ePsMatch,
                       ans_eDrsSMatch,
                       ans_ePsStrata,
                       ans_eDrsSStrata)[c("method","strata","A","events","denom")]
    ## Check
    expect_equal(SiteSummarySurv(df_prep1),
                 ans_combo)

})


test_that("site-specific summaries (binary and survival) are constructed correctly", {

    ## Combined df
    binary   <- SiteSummaryBin(df_prep1)
    survival <- SiteSummarySurv(df_prep1)
    binary$outcome   <- "binary"
    survival$outcome <- "survival"
    ans <- dplyr::bind_rows(binary,
                            survival)[c("outcome","method","strata","A","events","denom")]

    expect_equal(SiteSummary(df_prep1),
                 ans)


    ## Check handling of empty matched data and NA strata indicator
    ## Binary
    binaryNa   <- SiteSummaryBin(df_prep1_NA)
    methods_binaryNa <- unique(binary$method)[-1]
    ## Unadjusted part must match up in numbers
    expect_equal(binaryNa[binaryNa$method == "unadj",],
                 binary[binary$method == "unadj",
                        c("method","strata","A","events","denom")])
    ## The rest is all NA's
    binaryNaRest <- binaryNa[binaryNa$method %in% methods_binaryNa,
                             c("strata","A","events","denom")]
    rownames(binaryNaRest) <- NULL
    expect_equal(binaryNaRest,
                 data.frame(strata = as.integer(rep(NA, length(methods_binaryNa))),
                            A      = as.integer(rep(NA, length(methods_binaryNa))),
                            events = as.integer(rep(NA, length(methods_binaryNa))),
                            denom  = as.integer(rep(NA, length(methods_binaryNa)))))

    ## Survival
    survivalNa   <- SiteSummarySurv(df_prep1_NA)
    methods_survivalNa <- unique(survival$method)[-1]
    ## Unadjusted part must match up in numbers
    expect_equal(survivalNa[survivalNa$method == "unadj",],
                 survival[survival$method == "unadj",
                          c("method","strata","A","events","denom")])
    ## The rest is all NA's
    survivalNaRest <- survivalNa[survivalNa$method %in% methods_survivalNa,
                                 c("strata","A","events","denom")]
    rownames(survivalNaRest) <- NULL
    expect_equal(survivalNaRest,
                 data.frame(strata = as.integer(rep(NA, length(methods_survivalNa))),
                            A      = as.integer(rep(NA, length(methods_survivalNa))),
                            events = as.integer(rep(NA, length(methods_survivalNa))),
                            denom  = as.integer(rep(NA, length(methods_survivalNa)))))

})


###
### Functions for within-site risk set dataset construction
################################################################################

test_that("Site-specific risk set helper works correctly", {

    ## Construct
    df0 <- data.frame(time = c(1,2,3,4, rep(10,26)),
                      event = c(1,1,1,1, rep(0,26)),
                      A = c(1,0,1,0, rep(0,8), rep(1, 18)))
    out0 <- with(df0,
                 SiteRisksetsHelper(time      = time,
                                    event     = event,
                                    A         = A,
                                    eval_time = time[event == 1]))

    ## Expected riskset dataset
    riskset0 <- data.frame(eval_time = c(0,1,2,3,4),
                           events_A0   = c(0,0,1,0,1),
                           events_A1   = c(0,1,0,1,0),
                           riskset_A0 = c(10,10,10,9,9),
                           riskset_A1 = c(20,20,19,19,18))
    ## Check
    expect_equal(out0, riskset0)

    ## Evaluate at arbitrary time points
    out0b <- with(df0, SiteRisksetsHelper(time = time, event = event, A = A, eval_time = c(0,1,2,3,4,10)))
    ## Just two time points
    expect_equal(nrow(out0b), 6)
    expect_equal(out0b$eval_time,  c(0,1,2,3,4,10))
    expect_equal(out0b$events_A0,  c(0,0,1,0,1,0))
    expect_equal(out0b$events_A1,  c(0,1,0,1,0,0))
    expect_equal(out0b$riskset_A0, c(10,10,10,9,9,8))
    expect_equal(out0b$riskset_A1, c(20,20,19,19,18,18))


    ## Generate a risk set data retaining eval_time
    out1 <- with(df_prep1, SiteRisksetsHelper(time = time, event = event, A = A, eval_time = time[event == 1]))
    print(head(out1))

    ## It must be a new data frame
    expect_true("data.frame" %in% class(out1))
    ## event times must be sorted correctly
    expect_equal(out1$eval_time, sort(out1$eval_time))
    ## risksets must monotonically decrease
    expect_true(all(diff(out1$riskset_A0) <= 0))
    expect_true(all(diff(out1$riskset_A1) <= 0))
    ## There are as many rows as there are unique event times are in the original df + 1.
    ## +1 is for time point 0
    expect_true(nrow(out1) == length(unique(df_prep1$time[df_prep1$event == 1])) + 1)


    ## Generate a risk set data from binary outcome (no time variable)
    ## Here in this dummy dataset event also serves as Y.
    out2 <- with(df0, SiteRisksetsHelper(event = event, A = A))
    print(out2)

    ## Only one time point
    expect_true(nrow(out2) == 1)
    ## Event counts should match up
    expect_equal(out2$events_A0, with(df0, sum(event[A == 0])))
    expect_equal(out2$events_A1, with(df0, sum(event[A == 1])))
    ## Observation coutns should match up
    expect_equal(out2$riskset_A0, with(df0, length(event[A == 0])))
    expect_equal(out2$riskset_A1, with(df0, length(event[A == 1])))


    ## Generate a risk set data from binary outcome (no time variable)
    out2 <- with(df_prep1, SiteRisksetsHelper(event = Y, A = A))
    print(out2)

    ## Only one time point and it's 1
    expect_equal(nrow(out2), 1)
    expect_equal(out2$eval_time, 0)
    ## Event counts should match up
    expect_equal(out2$events_A0, with(df_prep1, sum(Y[A == 0])))
    expect_equal(out2$events_A1, with(df_prep1, sum(Y[A == 1])))
    ## Observation coutns should match up
    expect_equal(out2$riskset_A0, with(df_prep1, length(Y[A == 0])))
    expect_equal(out2$riskset_A1, with(df_prep1, length(Y[A == 1])))
    ## Further tests in the PS matched cohort
    expect_equal(as.numeric(with(subset(df_prep1, ePsMatch == 1),
                                 SiteRisksetsHelper(event = Y, A = A))[,c("events_A0", "events_A1")]),
                 as.numeric(with(subset(df_prep1, ePsMatch == 1), tapply(Y, A, sum))))
    ## Further tests in the DRS matched cohort
    expect_equal(as.numeric(with(subset(df_prep1, eDrsBMatch == 1),
                                 SiteRisksetsHelper(event = Y, A = A))[,c("events_A0", "events_A1")]),
                 as.numeric(with(subset(df_prep1, eDrsBMatch == 1), tapply(Y, A, sum))))

})


test_that("Site-specific risk set helper works correctly for weighted data", {

    ## Construct weight-of-0.5's dataset for sanity checking
    df0 <- data.frame(time = c(1,2,3,4, rep(10,26)),
                      event = c(1,1,1,1, rep(0,26)),
                      A = c(1,0,1,0, rep(0,8), rep(1, 18)),
                      W = rep(0.5, 30))
    out0 <- with(df0,
                 SiteRisksetsHelper(time      = time,
                                    event     = event,
                                    A         = A,
                                    W         = W,
                                    eval_time = time[event == 1]))

    ## Expected riskset dataset
    riskset0 <- data.frame(eval_time    = c(0,1,2,3,4),
                           events_A0    = c(0,0,1,0,1),
                           events_A1    = c(0,1,0,1,0),
                           riskset_A0   = c(10,10,10,9,9),
                           riskset_A1   = c(20,20,19,19,18),
                           ## Weighted results
                           w_events_A0  = 0.5 * c(0,0,1,0,1),
                           w_events_A1  = 0.5 * c(0,1,0,1,0),
                           w_riskset_A0 = 0.5 * c(10,10,10,9,9),
                           w_riskset_A1 = 0.5 * c(20,20,19,19,18),
                           ## Variance of weights
                           varw_events_A0  = 0,
                           varw_events_A1  = 0,
                           varw_nonevents_A0 = 0,
                           varw_nonevents_A1 = 0)
    ## Check
    expect_equal(out0, riskset0)

    ## Evaluate at arbitrary time points (This does not work for variance of weights)
    out0b <- with(df0,
                  SiteRisksetsHelper(time      = time,
                                     event     = event,
                                     A         = A,
                                     W         = W,
                                     eval_time = c(0,1,2,3,4,10)))
    ## Just two time points
    expect_equal(nrow(out0b), 6)
    expect_equal(out0b$eval_time,  c(0,1,2,3,4,10))
    ## Unweighted part of results
    expect_equal(out0b$events_A0,  c(0,0,1,0,1,0))
    expect_equal(out0b$events_A1,  c(0,1,0,1,0,0))
    expect_equal(out0b$riskset_A0, c(10,10,10,9,9,8))
    expect_equal(out0b$riskset_A1, c(20,20,19,19,18,18))
    ## Weighted part of results
    expect_equal(out0b$w_events_A0,  0.5 * c(0,0,1,0,1,0))
    expect_equal(out0b$w_events_A1,  0.5 * c(0,1,0,1,0,0))
    expect_equal(out0b$w_riskset_A0, 0.5 * c(10,10,10,9,9,8))
    expect_equal(out0b$w_riskset_A1, 0.5 * c(20,20,19,19,18,18))


    ## Construct a varying weight dataset
    W1 <- seq(0.1, 0.9, length.out = 30)
    df1 <- data.frame(time = c(1,2,3,4, rep(10,26)),
                      event = c(1,1,1,1, rep(0,26)),
                      A = c(1,0,1,0, rep(0,8), rep(1, 18)),
                      W = W1)
    out1 <- with(df1,
                 SiteRisksetsHelper(time      = time,
                                    event     = event,
                                    A         = A,
                                    W         = W,
                                    eval_time = time[event == 1]))
    ## Expected riskset dataset
    riskset1 <- data.frame(
        eval_time    = c(0,1,2,3,4),
        events_A0    = c(0,0,1,0,1),
        events_A1    = c(0,1,0,1,0),
        riskset_A0   = c(10,10,10,9,9),
        riskset_A1   = c(20,20,19,19,18),
        ## Weighted results
        w_events_A0  = c(0,0,1*W1[2],0,1*W1[4]),
        w_events_A1  = c(0,1*W1[1],0,1*W1[3],0),
        w_riskset_A0 = c(sum(subset(df1, A == 0 & time >= 0)$W),
                         sum(subset(df1, A == 0 & time >= 1)$W),
                         sum(subset(df1, A == 0 & time >= 2)$W),
                         sum(subset(df1, A == 0 & time >= 3)$W),
                         sum(subset(df1, A == 0 & time >= 4)$W)),
        w_riskset_A1 = c(sum(subset(df1, A == 1 & time >= 0)$W),
                         sum(subset(df1, A == 1 & time >= 1)$W),
                         sum(subset(df1, A == 1 & time >= 2)$W),
                         sum(subset(df1, A == 1 & time >= 3)$W),
                         sum(subset(df1, A == 1 & time >= 4)$W)),
        ## Variance of weights
        ## Event part is 0 if no ties
        varw_events_A0  = 0,
        varw_events_A1  = 0,
        ## Variance of weights among nonevents in risk sets
        ## Choose those who are in risk set t, but not having event at t.
        varw_nonevents_A0 = c(var(subset(df1, A == 0 & time >= 0 & !(time == 0 & event == 1))$W),
                              var(subset(df1, A == 0 & time >= 1 & !(time == 1 & event == 1))$W),
                              var(subset(df1, A == 0 & time >= 2 & !(time == 2 & event == 1))$W),
                              var(subset(df1, A == 0 & time >= 3 & !(time == 3 & event == 1))$W),
                              var(subset(df1, A == 0 & time >= 4 & !(time == 4 & event == 1))$W)),
        varw_nonevents_A1 = c(var(subset(df1, A == 1 & time >= 0 & !(time == 0 & event == 1))$W),
                              var(subset(df1, A == 1 & time >= 1 & !(time == 1 & event == 1))$W),
                              var(subset(df1, A == 1 & time >= 2 & !(time == 2 & event == 1))$W),
                              var(subset(df1, A == 1 & time >= 3 & !(time == 3 & event == 1))$W),
                              var(subset(df1, A == 1 & time >= 4 & !(time == 4 & event == 1))$W)))
    ## Check
    expect_equal(out1, riskset1)


    ## Generate a risk set data retaining eval_time (real IPTW survival data)
    out1 <- with(df_prep1, SiteRisksetsHelper(time = time,
                                              event = event,
                                              A = A,
                                              W = ePsSIptw,
                                              eval_time = time[event == 1]))
    print(head(out1))
    print(tail(out1))
    ##
    ## It must be a new data frame
    expect_true("data.frame" %in% class(out1))
    ## event times must be sorted correctly
    expect_equal(out1$eval_time, sort(out1$eval_time))
    ## risksets must monotonically decrease
    expect_true(all(diff(out1$riskset_A0) <= 0))
    expect_true(all(diff(out1$riskset_A1) <= 0))
    expect_true(all(diff(out1$w_riskset_A0) <= 0))
    expect_true(all(diff(out1$w_riskset_A1) <= 0))
    ## There are as many rows as there are unique event times are in the original df + 1.
    ## +1 is for time point 0
    expect_true(nrow(out1) == length(unique(df_prep1$time[df_prep1$event == 1])) + 1)
    ## Variance of weights must be 0 for risk sets without any events
    expect_equal(subset(out1, events_A1 == 0)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 0)$varw_events_A1)))
    expect_equal(subset(out1, events_A1 == 0)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 0)$varw_events_A1)))
    ## Variance of weights must be 0 for risk sets without only one event
    expect_equal(subset(out1, events_A1 == 1)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 1)$varw_events_A1)))
    expect_equal(subset(out1, events_A1 == 1)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 1)$varw_events_A1)))
    ## The event count on last day should be small
    expect_true(tail(out1, 1)$events_A0 < 10)
    expect_true(tail(out1, 1)$events_A1 < 10)
    expect_true(tail(out1, 1)$w_events_A0 < 10)
    expect_true(tail(out1, 1)$w_events_A1 < 10)
    ## Weighted and unweighted are close in this uncounfounded setting
    expect_true(all(out1$w_events_A0 - out1$events_A0 <= 1))
    expect_true(all(out1$w_events_A1 - out1$events_A1 <= 1))
    expect_true(all(out1$w_riskset_A0 - out1$riskset_A0 <= 4))
    expect_true(all(out1$w_riskset_A1 - out1$riskset_A1 <= 1))


    ## Generate a risk set data retaining eval_time (real MW survival data)
    out1 <- with(df_prep1, SiteRisksetsHelper(time = time,
                                              event = event,
                                              A = A,
                                              W = ePsMw,
                                              eval_time = time[event == 1]))
    print(head(out1))
    print(tail(out1))
    ##
    ## It must be a new data frame
    expect_true("data.frame" %in% class(out1))
    ## event times must be sorted correctly
    expect_equal(out1$eval_time, sort(out1$eval_time))
    ## risksets must monotonically decrease
    expect_true(all(diff(out1$riskset_A0) <= 0))
    expect_true(all(diff(out1$riskset_A1) <= 0))
    expect_true(all(diff(out1$w_riskset_A0) <= 0))
    expect_true(all(diff(out1$w_riskset_A1) <= 0))
    ## There are as many rows as there are unique event times are in the original df + 1.
    ## +1 is for time point 0
    expect_true(nrow(out1) == length(unique(df_prep1$time[df_prep1$event == 1])) + 1)
    ## Variance of weights must be 0 for risk sets without any events
    expect_equal(subset(out1, events_A1 == 0)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 0)$varw_events_A1)))
    expect_equal(subset(out1, events_A1 == 0)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 0)$varw_events_A1)))
    ## Variance of weights must be 0 for risk sets without only one event
    expect_equal(subset(out1, events_A1 == 1)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 1)$varw_events_A1)))
    expect_equal(subset(out1, events_A1 == 1)$varw_events_A1,
                 rep(0, length(subset(out1, events_A1 == 1)$varw_events_A1)))
    ## The event count on last day should be small
    expect_true(tail(out1, 1)$events_A0 < 10)
    expect_true(tail(out1, 1)$events_A1 < 10)
    expect_true(tail(out1, 1)$w_events_A0 < 10)
    expect_true(tail(out1, 1)$w_events_A1 < 10)
    ## Weighted and unweighted are close in this uncounfounded setting
    expect_true(all(out1$w_events_A0 - out1$events_A0 <= 1))
    expect_true(all(out1$w_events_A1 - out1$events_A1 <= 1))
    expect_true(all(out1$w_riskset_A0 - out1$riskset_A0 <= 1))
    expect_true(all(out1$w_riskset_A1 - out1$riskset_A1 <= 1))


    ## Generate a risk set data from binary outcome (no time variable)
    ## Here in this dummy dataset event also serves as Y.
    out2 <- with(df0, SiteRisksetsHelper(event = event, A = A, W = W))
    print(out2)

    ## Only one time point
    expect_true(nrow(out2) == 1)
    ## Event counts (or weighted counts) should match up
    expect_equal(out2$events_A0, with(df0, sum(event[A == 0])))
    expect_equal(out2$events_A1, with(df0, sum(event[A == 1])))
    expect_equal(out2$w_events_A0, 0.5 * with(df0, sum(event[A == 0])))
    expect_equal(out2$w_events_A1, 0.5 * with(df0, sum(event[A == 1])))
    ## Observation counts (or weighted counts) should match up
    expect_equal(out2$riskset_A0, with(df0, length(event[A == 0])))
    expect_equal(out2$riskset_A1, with(df0, length(event[A == 1])))
    expect_equal(out2$w_riskset_A0, 0.5 * with(df0, length(event[A == 0])))
    expect_equal(out2$w_riskset_A1, 0.5 * with(df0, length(event[A == 1])))


    ## Generate a risk set data from binary outcome (no time variable)
    out2 <- with(df_prep1, SiteRisksetsHelper(event = Y, A = A, W = ePsSIptw))
    print(out2)

    ## Only one time point and it's 1
    expect_equal(nrow(out2), 1)
    expect_equal(out2$eval_time, 0)
    ## Event counts should match up
    expect_equal(out2$w_events_A0, with(df_prep1, sum(Y[A == 0] * ePsSIptw[A == 0])))
    expect_equal(out2$w_events_A1, with(df_prep1, sum(Y[A == 1] * ePsSIptw[A == 1])))
    ## Observation coutns should match up
    expect_equal(out2$riskset_A0, with(df_prep1, length(Y[A == 0] * ePsSIptw[A == 0])))
    expect_equal(out2$riskset_A1, with(df_prep1, length(Y[A == 1] * ePsSIptw[A == 1])))
    ## Further tests in the PS matched cohort
    expect_equal(as.numeric(
        with(subset(df_prep1, ePsMatch == 1),
             SiteRisksetsHelper(event = Y, A = A, W = ePsSIptw))[,c("w_events_A0", "w_events_A1")]),
        as.numeric(with(subset(df_prep1, ePsMatch == 1),
                        tapply(Y * ePsSIptw, A, sum))))
    ## Further tests in the DRS matched cohort
    expect_equal(as.numeric(
        with(subset(df_prep1, eDrsBMatch == 1),
             SiteRisksetsHelper(event = Y, A = A, W = ePsSIptw))[,c("w_events_A0", "w_events_A1")]),
        as.numeric(with(subset(df_prep1, eDrsBMatch == 1),
                        tapply(Y * ePsSIptw, A, sum))))

})


test_that("site-specific risk sets (binary) are correctly created with or without strata", {

    vars_relevant <- c("strata","eval_time","events_A0","events_A1","riskset_A0","riskset_A1")

    ## Binary split by PS strata
    ans_ePsStrata <- split(x = df_prep1, f = df_prep1$ePsStrata) %>%
        lapply(FUN = function(df) {
            ## work on stratum by stratum
            with(df, SiteRisksetsHelper(event = Y, A = A))
        }) %>% do.call(what = rbind)
    ans_ePsStrata$strata <- seq_len(nrow(ans_ePsStrata))
    ans_ePsStrata <- ans_ePsStrata[vars_relevant]
    rownames(ans_ePsStrata) <- NULL
    ##
    expect_equal(with(df_prep1,
                      SiteRisksetsByStrata(event = Y, A = A, strata = ePsStrata)),
                 ans_ePsStrata)


    ## Binary split by DRS strata
    ans_eDrsBStrata <- split(x = df_prep1, f = df_prep1$eDrsBStrata) %>%
        lapply(FUN = function(df) {
            ## work on stratum by stratum
            with(df, SiteRisksetsHelper(event = Y, A = A))
        }) %>% do.call(what = rbind)
    ans_eDrsBStrata$strata <- seq_len(nrow(ans_eDrsBStrata))
    ans_eDrsBStrata <- ans_eDrsBStrata[vars_relevant]
    rownames(ans_eDrsBStrata) <- NULL
    ##
    expect_equal(with(df_prep1,
                      SiteRisksetsByStrata(event = Y, A = A, strata = eDrsBStrata)),
                 ans_eDrsBStrata)


    ## Test combined results
    res  <- SiteRisksetsBin(df_prep1)

    res_ePsStrata <- subset(res, method == "ePsStrata")[vars_relevant]
    rownames(res_ePsStrata) <- NULL
    expect_equal(res_ePsStrata,
                 ans_ePsStrata)

    res_eDrsBStrata <- subset(res, method == "eDrsBStrata")[vars_relevant]
    rownames(res_eDrsBStrata) <- NULL
    expect_equal(res_eDrsBStrata,
                 ans_eDrsBStrata)


    ## Without strata
    ## Unadjusted
    ans_unadj <- with(df_prep1, SiteRisksetsHelper(event = Y, A = A))
    ans_unadj$strata <- as.numeric(NA)
    ans_unadj <- ans_unadj[vars_relevant]
    rownames(ans_unadj) <- NULL
    res_unadj <- subset(res, method == "unadj")[vars_relevant]
    rownames(res_unadj) <- NULL
    expect_equal(res_unadj,
                 ans_unadj)
    ## ePsMatch
    ans_ePsMatch <- with(subset(df_prep1, ePsMatch == 1), SiteRisksetsHelper(event = Y, A = A))
    ans_ePsMatch$strata <- as.numeric(NA)
    ans_ePsMatch <- ans_ePsMatch[vars_relevant]
    rownames(ans_ePsMatch) <- NULL
    res_ePsMatch <- subset(res, method == "ePsMatch")[vars_relevant]
    rownames(res_ePsMatch) <- NULL
    expect_equal(res_ePsMatch,
                 ans_ePsMatch)
    ## eDrsBMatch
    ans_eDrsBMatch <- with(subset(df_prep1, eDrsBMatch == 1), SiteRisksetsHelper(event = Y, A = A))
    ans_eDrsBMatch$strata <- as.numeric(NA)
    ans_eDrsBMatch <- ans_eDrsBMatch[vars_relevant]
    rownames(ans_eDrsBMatch) <- NULL
    res_eDrsBMatch <- subset(res, method == "eDrsBMatch")[vars_relevant]
    rownames(res_eDrsBMatch) <- NULL
    expect_equal(res_eDrsBMatch,
                 ans_eDrsBMatch)
    ## ePsSIptw
    ans_ePsSIptw <- with(df_prep1, SiteRisksetsHelper(event = Y, A = A, W = ePsSIptw))
    ans_ePsSIptw$strata <- as.numeric(NA)
    ans_ePsSIptw <- ans_ePsSIptw[vars_relevant]
    rownames(ans_ePsSIptw) <- NULL
    res_ePsSIptw <- subset(res, method == "ePsSIptw")[vars_relevant]
    rownames(res_ePsSIptw) <- NULL
    expect_equal(res_ePsSIptw,
                 ans_ePsSIptw)
    ## ePsMw
    ans_ePsMw <- with(df_prep1, SiteRisksetsHelper(event = Y, A = A, W = ePsMw))
    ans_ePsMw$strata <- as.numeric(NA)
    ans_ePsMw <- ans_ePsMw[vars_relevant]
    rownames(ans_ePsMw) <- NULL
    res_ePsMw <- subset(res, method == "ePsMw")[vars_relevant]
    rownames(res_ePsMw) <- NULL
    expect_equal(res_ePsMw,
                 ans_ePsMw)


    ## Length zero vector handling
    ## Empty matched binary data
    expect_equal(SiteRisksetsByStrata(event = as.integer(), A = as.integer()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer()),
                        class),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1)),
                        class))
    ## Check names
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer()),
                        names),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1)),
                        names))

    ## Empty stratified binary data
    expect_equal(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), strata = as.integer()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), strata = as.integer()),
                        class),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1), strata = as.integer(1)),
                        class))
    ## Check names
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), strata = as.integer()),
                        names),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1), strata = as.integer(1)),
                        names))

    ## Empty weighted binary data
    expect_equal(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), W = as.double()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1",
                                       "w_events_A0", "w_events_A1", "w_riskset_A0", "w_riskset_A1",
                                       "varw_events_A0", "varw_events_A1", "varw_nonevents_A0", "varw_nonevents_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer,
                                         as.double, as.double, as.double, as.double,
                                         as.double, as.double, as.double, as.double),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), W = as.double()),
                        class),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1), W = as.double(1)),
                        class))
    ## Check names
    expect_equal(sapply(SiteRisksetsByStrata(event = as.integer(), A = as.integer(), W = as.double()),
                        names),
                 sapply(SiteRisksetsByStrata(event = as.integer(1), A = as.integer(1), W = as.double(1)),
                        names))

})


test_that("site-specific risk sets (survival) are correctly created with or without strata", {

    vars_relevant <- c("strata","eval_time","events_A0","events_A1","riskset_A0","riskset_A1")

    ## Survival split by PS strata
    ans_ePsStrata <- split(x = df_prep1, f = df_prep1$ePsStrata) %>%
        lapply(FUN = function(df) {
            ## work on stratum by stratum
            with(df,
                 SiteRisksetsHelper(time = time, event = event, A = A, eval_time = 0:366))
        }) %>% do.call(what = rbind)
    ans_ePsStrata$strata <- rep(1:10, each = length(0:366))
    ans_ePsStrata <- ans_ePsStrata[vars_relevant]
    ## Drop no event rows except for time 0
    ans_ePsStrata <- subset(ans_ePsStrata, events_A0 + events_A1 > 0 | eval_time == 0)
    rownames(ans_ePsStrata) <- NULL
    ##
    expect_equal(with(df_prep1,
                      SiteRisksetsByStrata(time = time, event = event, A = A,
                                           strata = ePsStrata)),
                 ans_ePsStrata)


    ## Survival split by DRS strata
    ans_eDrsSStrata <- split(x = df_prep1, f = df_prep1$eDrsSStrata) %>%
        lapply(FUN = function(df) {
            ## work on stratum by stratum
            with(df,
                 SiteRisksetsHelper(time = time, event = event, A = A, eval_time = 0:366))
        }) %>% do.call(what = rbind)
    ans_eDrsSStrata$strata <- rep(1:10, each = length(0:366))
    ans_eDrsSStrata <- ans_eDrsSStrata[vars_relevant]
    rownames(ans_eDrsSStrata) <- NULL
    ## Drop no event rows except for time 0
    ans_eDrsSStrata <- subset(ans_eDrsSStrata, events_A0 + events_A1 > 0 | eval_time == 0)
    rownames(ans_eDrsSStrata) <- NULL
    ##
    expect_equal(with(df_prep1,
                      SiteRisksetsByStrata(time = time, event = event, A = A,
                                           strata = eDrsSStrata)),
                 ans_eDrsSStrata)


    ## Test combined results
    res <- SiteRisksetsSurv(df_prep1)

    res_ePsStrata <- subset(res, method == "ePsStrata")[vars_relevant]
    rownames(res_ePsStrata) <- NULL
    expect_equal(res_ePsStrata,
                 ans_ePsStrata)

    res_eDrsSStrata <- subset(res, method == "eDrsSStrata")[vars_relevant]
    rownames(res_eDrsSStrata) <- NULL
    expect_equal(res_eDrsSStrata,
                 ans_eDrsSStrata)

    ## Check if there are unusually many events on the administrative censoring day
    expect_equal(sum(subset(res, eval_time == 366)[,"events_A0"] > 10), 0)
    expect_equal(sum(subset(res, eval_time == 366)[,"events_A1"] > 10), 0)
    expect_equal(sum(subset(res, eval_time == 366 & method %in% c("ePsSIptw","ePsMw"))[,"w_events_A0"] > 10), 0)


    ## Without strata
    ## Unadjusted
    ans_unadj <- with(df_prep1,
                      SiteRisksetsHelper(time = time, event = event, A = A, eval_time = 0:366))
    ans_unadj$strata <- as.numeric(NA)
    ans_unadj <- ans_unadj[vars_relevant]
    rownames(ans_unadj) <- NULL
    res_unadj <- subset(res, method == "unadj")[vars_relevant]
    rownames(res_unadj) <- NULL
    expect_equal(res_unadj,
                 ans_unadj)
    ## ePsMatch
    ans_ePsMatch <- with(subset(df_prep1, ePsMatch == 1),
                         SiteRisksetsHelper(time = time, event = event, A = A, eval_time = 0:366))
    ans_ePsMatch$strata <- as.numeric(NA)
    ans_ePsMatch <- ans_ePsMatch[vars_relevant]
    rownames(ans_ePsMatch) <- NULL
    res_ePsMatch <- subset(res, method == "ePsMatch")[vars_relevant]
    rownames(res_ePsMatch) <- NULL
    expect_equal(res_ePsMatch,
                 ans_ePsMatch)
    ## eDrsSMatch
    ans_eDrsSMatch <- with(subset(df_prep1, eDrsSMatch == 1),
                           SiteRisksetsHelper(time = time, event = event, A = A, eval_time = 0:366))
    ans_eDrsSMatch$strata <- as.numeric(NA)
    ans_eDrsSMatch <- ans_eDrsSMatch[vars_relevant]
    rownames(ans_eDrsSMatch) <- NULL
    res_eDrsSMatch <- subset(res, method == "eDrsSMatch")[vars_relevant]
    rownames(res_eDrsSMatch) <- NULL
    expect_equal(res_eDrsSMatch,
                 ans_eDrsSMatch)
    ## ePsSIptw
    ans_ePsSIptw <- with(df_prep1, SiteRisksetsHelper(time = time, event = event, A = A,
                                                      W = ePsSIptw, eval_time = 0:366))
    ans_ePsSIptw$strata <- as.numeric(NA)
    ans_ePsSIptw <- ans_ePsSIptw[vars_relevant]
    rownames(ans_ePsSIptw) <- NULL
    res_ePsSIptw <- subset(res, method == "ePsSIptw")[vars_relevant]
    rownames(res_ePsSIptw) <- NULL
    expect_equal(res_ePsSIptw,
                 ans_ePsSIptw)
    ## ePsMw
    ans_ePsMw <- with(df_prep1, SiteRisksetsHelper(time = time, event = event, A = A,
                                                   W = ePsMw, eval_time = 0:366))
    ans_ePsMw$strata <- as.numeric(NA)
    ans_ePsMw <- ans_ePsMw[vars_relevant]
    rownames(ans_ePsMw) <- NULL
    res_ePsMw <- subset(res, method == "ePsMw")[vars_relevant]
    rownames(res_ePsMw) <- NULL
    expect_equal(res_ePsMw,
                 ans_ePsMw)


    ## Length zero vector handling
    ## Empty matched survival data
    expect_equal(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer()),
                        class),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1)),
                        class))
    ## Check names
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer()),
                        names),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1)),
                        names))

    ## Empty stratified survival data
    expect_equal(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(),
                                      strata = as.integer()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(),
                                             strata = as.integer()),
                        class),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1),
                                             strata = as.integer(1)),
                        class))
    ## Check Names
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(),
                                             strata = as.integer()),
                        names),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1),
                                             strata = as.integer(1)),
                        names))

    ## Empty weighted survival data
    expect_equal(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(), W = as.double()),
                 CreateNaDf(labels = c("strata","eval_time",
                                       "events_A0","events_A1","riskset_A0","riskset_A1",
                                       "w_events_A0", "w_events_A1", "w_riskset_A0", "w_riskset_A1",
                                       "varw_events_A0", "varw_events_A1", "varw_nonevents_A0", "varw_nonevents_A1"),
                            as_funcs = c(as.integer, as.double,
                                         as.integer, as.integer, as.integer, as.integer,
                                         as.double, as.double, as.double, as.double,
                                         as.double, as.double, as.double, as.double),
                            nrow = 1))
    ## Check classes
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(), W = as.double()),
                        class),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1), W = as.double(1)),
                        class))
    ## Check names
    expect_equal(sapply(SiteRisksetsByStrata(time = as.double(), event = as.integer(), A = as.integer(), W = as.double()),
                        names),
                 sapply(SiteRisksetsByStrata(time = as.double(1), event = as.integer(1), A = as.integer(1), W = as.double(1)),
                        names))

})


test_that("site-specific risk sets (binary and survival) are constructed correctly", {

    ## Binary
    binary <- SiteRisksetsBin(df_prep1)
    ## All method are found
    methods_bin <- c("unadj","ePsMatch", "eDrsBMatch", "ePsStrata", "eDrsBStrata",
                     "ePsSIptw", "ePsMw")
    expect_true(all(binary$method %in% methods_bin))
    expect_true(all(methods_bin %in% binary$method))
    ## For binary evaluation time is all zero
    expect_equal(unique(binary$eval_time), 0)
    ## Strata should be NA for non-stratification methods
    expect_true(all(is.na(binary[!grepl("Strata", binary$method), "strata"])))
    ## Weighted parts are NA for non-weighting methods
    expect_true(all(is.na(binary[!grepl("Iptw|Mw", binary$method),
                                 c("w_events_A0", "w_events_A1",
                                   "w_riskset_A0", "w_riskset_A1")])))

    ## Survival
    survival <- SiteRisksetsSurv(df_prep1)
    ## All method are found
    methods_surv <- c("unadj","ePsMatch", "eDrsSMatch", "ePsStrata", "eDrsSStrata",
                      "ePsSIptw", "ePsMw")
    expect_true(all(survival$method %in% methods_surv))
    expect_true(all(methods_surv %in% survival$method))
    ## Strata should be NA for non-stratification methods
    expect_true(all(is.na(survival[!grepl("Strata", survival$method), "strata"])))
    ## Weighted parts are NA for non-weighting methods
    expect_true(all(is.na(binary[!grepl("Iptw|Mw", binary$method),
                                 c("w_events_A0", "w_events_A1",
                                   "w_riskset_A0", "w_riskset_A1")])))

    ## Combined df
    binary$outcome   <- "binary"
    survival$outcome <- "survival"
    ans <- dplyr::bind_rows(binary,
                            survival)[c("outcome","method","strata",
                                        "eval_time",
                                        "events_A0", "events_A1",
                                        "riskset_A0", "riskset_A1",
                                        ## Weighted part
                                        "w_events_A0", "w_events_A1",
                                        "w_riskset_A0", "w_riskset_A1",
                                        ## Weighted part
                                        "varw_events_A0", "varw_events_A1",
                                        "varw_nonevents_A0", "varw_nonevents_A1")]

    res <- SiteRisksets(df_prep1)
    expect_equal(res, ans)


    ## Examine NA handling
    ## binary
    ans_binaryNA <- binary[!duplicated(binary$method),
                           !(names(binary) %in% "outcome")]
    ## Put NA's except for unadj row and method column
    ans_binaryNA[!(ans_binaryNA$method %in% "unadj"),
                 !(names(ans_binaryNA) %in% "method")] <- NA
    rownames(ans_binaryNA) <- NULL
    ## Test
    expect_equal(SiteRisksetsBin(df_prep1_NA),
                 ans_binaryNA)

    ## survival (keep all unadjusted rows and non-duplicated methods)
    ans_survivalNA <- survival[(survival$method %in% "unadj" | !duplicated(survival$method)),
                               !(names(survival) %in% "outcome")]
    ## Put NA's except for unadj row and method column
    ans_survivalNA[!(ans_survivalNA$method %in% "unadj"),
                   !(names(ans_survivalNA) %in% "method")] <- NA
    rownames(ans_survivalNA) <- NULL
    ## Test
    expect_equal(SiteRisksetsSurv(df_prep1_NA),
                 ans_survivalNA)

    ## Both outcomes (keep all unadjusted rows and non-duplicated methods)
    ansNA <- res[(res$method %in% "unadj" | !duplicated(res[c("outcome","method")])), ]
    ansNA[!(ansNA$method %in% "unadj"),
          !(names(ansNA) %in% c("outcome","method"))] <- NA
    rownames(ansNA) <- NULL
    expect_equal(SiteRisksets(df_prep1_NA),
                 ansNA)

})


###
### Functions for individual-level data construction
################################################################################

test_that("individual-level datasets are extracted correctly", {

    ## Extract observed and derived variables except for covariates
    ans <- df_prep1[,c("A", "Y", "time", "event",
                       "ePs", "eDrsB", "eDrsS",
                       "ePsSIptw", "ePsMw",
                       "ePsStrata", "eDrsBStrata", "eDrsSStrata",
                       "ePsMatch", "eDrsBMatch", "eDrsSMatch")]

    expect_equal(SiteDataset(df_prep1),
                 ans)
})


test_that("individual-level counterfactual datasets are extracted correctly", {

    expect_equal(SiteTruth(df_prep1),
                 ## Returned as is
                 df_prep1)
})
