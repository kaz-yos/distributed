################################################################################
### Test data generator functions
##
## Created on: 2016-05-05
## Author: Kazuki Yoshida
################################################################################

library(testthat)

### Context (1 for each file)
context("### Test 01.GenerateData.R")


set.seed(20160505)

mean1 <- 0
sd1   <- 1
pp    <- seq(0.1, 0.6, 0.1)
dfX   <- AssignCovariates(n = 10^4, mean1 = mean1, sd1 = sd1, pp = pp)
cat("\n")
print(head(dfX, 20))


###
### Variable generation
################################################################################

test_that("Covariates are generated correctly", {

    expect_equal(class(dfX), "data.frame")
    expect_equal(dim(dfX), c(10^4, 8))
    expect_equal(names(dfX), c("id", paste0("X", 1:7)))
    expect_equal(round(sd(dfX$X1)), sd1)
    expect_equal(round(mean(dfX$X1)), mean1)
    expect_equal(as.numeric(round(colMeans(dfX[,paste0("X", 2:7)]),1)),
                 pp)
})


test_that("Treatment is assigned correctly", {

    ## No covariate effect scenario (1:1 randomization)
    alpha0 <- 0
    alphaX <- rep(0,7)
    dfXA <- AssignTreatment(dfX = dfX, alpha0 = alpha0, alphaX = alphaX)
    cat("\n")
    print(head(dfXA, 20))

    expect_true(all(dfXA$lpA == 0))
    expect_true(all(dfXA$pA == 0.5))
    expect_true(round(mean(dfXA$A), 1) == 0.5)


    ## Strong positive effect of all binary covariates
    alpha0 <- 10
    alphaX <- c(0, rep(10,6))
    dfXA <- AssignTreatment(dfX = dfX, alpha0 = alpha0, alphaX = alphaX)
    cat("\n")
    print(head(dfXA, 20))
    ## Mostly 1
    expect_true(round(mean(dfXA$A), 1) >= 0.95)


    ## Strong negative effect of all binary covariates
    alpha0 <- -10
    alphaX <- c(0, rep(-10,6))
    dfXA <- AssignTreatment(dfX = dfX, alpha0 = alpha0, alphaX = alphaX)
    cat("\n")
    print(head(dfXA, 20))
    ## Mostly 1
    expect_true(round(mean(dfXA$A), 1) <= 0.05)

})


test_that("Disease is assigned correctly", {

    ## Randomize treatment
    alpha0 <- 0
    alphaX <- rep(0,7)
    dfXA <- AssignTreatment(dfX = dfX, alpha0 = alpha0, alphaX = alphaX)
    cat("\n")
    print(head(dfXA, 20))


    ## Randomize outcome
    beta0 <- 0
    betaX <- rep(0,7)
    betaA <- 0
    betaXA <- rep(0,7)
    dfXAY <- AssignOutcomeBin(dfXA = dfXA, beta0 = beta0, betaX = betaX, betaA = betaA, betaXA = betaXA)
    cat("\n")
    print(head(dfXAY, 20))

    expect_true(all(dfXAY$lpY == 0))
    expect_true(all(dfXAY$trueDrs == 0.5))
    expect_true(round(mean(dfXAY$Y), 1) == 0.5)

    ## Strong positive effect of all binary covariates
    beta0 <- 10
    betaX <- c(0, rep(10,6))
    dfXAY <- AssignOutcomeBin(dfXA = dfXA, beta0 = beta0, betaX = betaX, betaA = betaA, betaXA = betaXA)
    cat("\n")
    print(head(dfXAY, 20))
    ## Mostly 1
    expect_true(mean(dfXAY$Y) >= 0.95)


    ## Strong negative effect of all binary covariates
    beta0 <- -10
    betaX <- c(0, rep(-10,6))
    dfXAY <- AssignOutcomeBin(dfXA = dfXA, beta0 = beta0, betaX = betaX, betaA = betaA, betaXA = betaXA)
    cat("\n")
    print(head(dfXAY, 20))
    ## Mostly 1
    expect_true(mean(dfXAY$Y) <= 0.05)

})


test_that("Time to event is assigned correctly", {

    ## Randomize treatment
    alpha0 <- 0
    alphaX <- rep(0,7)
    dfXA <- AssignTreatment(dfX = dfX, alpha = alpha0, alphaX = alphaX)
    cat("\n")
    print(head(dfXA, 20))


    ## Randomize outcome (no covariate or treatment effect)
    betaX <- rep(0,7)
    betaA <- 0
    betaXA <- rep(0,7)
    dfXAY <- AssignOutcomeSurv(dfXA = dfXA, betaX = betaX, betaA = betaA, betaXA = betaXA,
                               lambda = 0.1, lambda_c = 0.05, Tmax = 10^5)
    cat("\n")
    print(head(dfXAY, 20))
    expect_true(length(unique(dfXAY$C)) > 1)
    expect_true(all(dfXAY$lpLH == 0))
    ## Under no treatment effect both counterfactuals are the same
    expect_equal(dfXAY$rate, dfXAY$rate0)
    expect_equal(dfXAY$rate, dfXAY$rate1)
    ## Under no covariate effect at all all rates are the same
    expect_true(all(dfXAY$rate == 0.1))
    ## Mean true event time
    expect_true(round(mean(dfXAY$T)) == 10)
    ## Mean true censoring time
    expect_true(round(mean(dfXAY$C)) == 20)
    ## Under no treatment effect: censoring prob = lambda_c / (lambda + lambda_c)
    expect_true(abs(mean(dfXAY$event == 0) - 1/3) < 0.03)


    ## Administrative censoring
    dfXAY <- AssignOutcomeSurv(dfXA = dfXA, betaX = betaX, betaA = betaA, betaXA = betaXA,
                               lambda = 0.1, lambda_c = 0.05, Tmax = 1)
    ## All observed times should be shorter than 1 year (366 by ceiling).
    expect_true(all(dfXAY$time <= 366))


    ## Strong positive effect of all binary covariates
    beta0 <- 10
    betaX <- c(0, rep(10,6))
    dfXAY <- AssignOutcomeSurv(dfXA = dfXA, betaX = betaX, betaA = betaA, betaXA = betaXA,
                               lambda = 0.1, lambda_c = 0.005, Tmax = 10^5)
    cat("\n")
    print(head(dfXAY, 20))
    ## Under no treatment effect both counterfactuals are the same
    expect_equal(dfXAY$rate, dfXAY$rate0)
    expect_equal(dfXAY$rate, dfXAY$rate1)
    ## Events happen quickly
    expect_true(mean(dfXAY$T < 1) >= 0.90)
    ## Censoring probability
    expect_true(abs(mean(dfXAY$event == 0) - 0.005/(0.1+0.005)) < 0.03)


    ## Strong negative effect of all binary covariates
    beta0 <- -10
    betaX <- c(0, rep(-10,6))
    dfXAY <- AssignOutcomeSurv(dfXA = dfXA, betaX = betaX, betaA = betaA, betaXA = betaXA,
                               lambda = 0.1, lambda_c = 0.005, Tmax = 10^22)
    cat("\n")
    print(head(dfXAY, 20))
    ## Under no treatment effect both counterfactuals are the same
    expect_equal(dfXAY$rate, dfXAY$rate0)
    expect_equal(dfXAY$rate, dfXAY$rate1)
    ## Events happen slowly
    expect_true(mean(dfXAY$T < 1) <= 0.10)
    ## Censoring probability
    expect_true(abs(mean(dfXAY$event == 0) - 0.005/(0.1+0.005)) < 0.03)

})

###
### Single site generation
################################################################################

test_that("A single center is generated correctly (also table generation)", {

    ## Data generation
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(0,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    cat("\n")
    print.data.frame(head(df, 20))
    ## print method
    print(df)
    expect_output(print(df),
                  "Size: 100000")
    expect_output(print(df),
                  "not ready")
    ## data structure
    expect_true("ResSite" %in% class(df))
    expect_true("data.frame" %in% class(df))
    expect_equal(attr(df, "ScenarioResSite"),
                 list(n = 10^5,
                      alphas = c(alpha0 = 0,
                                 alphaX = rep(0,7)),
                      betas = c(beta0 = 0,
                                betaX = rep(0,7),
                                betaA = 0,
                                betaXA = rep(0,7)),
                      survParams = c(10, 10, 1)))
    expect_true(nrow(df) == 10^5)

    ## Different covariate counts
    p <- 16
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(0,p)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,p),
                                      betaA = 0,
                                      betaXA = rep(0,p)),
                            survParams = c(10, 10, 1))
    expect_equal(p, length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(df))))

    ## summary method
    ## params empty if truth = FALSE is specified
    expect_true(is.null(summary(df)$params))
    ## params populated if truth = TRUE is specified
    expect_true(!is.null(summary(df, truth = TRUE)$params))

    ## Data generation and tabling
    tab <- GenerateOneCenterTable(n = 10^5,
                                  alphas = c(alpha0 = 0,
                                             alphaX = rep(0,7)),
                                  betas = c(beta0 = 0,
                                            betaX = rep(0,7),
                                            betaA = 0,
                                            betaXA = rep(0,7)),
                                  survParams = c(10, 10, 1))
    cat("\n")
    print(tab)
    expect_equal(class(tab), "ParamsTableOne")
    expect_equal(class(tab$TableOneOverall), "TableOne")
    expect_equal(class(tab$TableOneByA), "TableOne")
    expect_equal(ncol(print(tab, print = FALSE)), 4)

})


test_that("Check data generation under various scenarios", {

    ## Randomized treatment and outcome
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(0,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    ## Names of cvariates and counter factuals
    nameCovCounter <- Filter(f = function(elt) {grepl("^X|^pY|^rate", elt)}, x = names(df))
    ## Treatment assignment is roughtly 50:50
    expect_equal(as.numeric(round(diff(table(df$A) / 10^3))), 0)
    ## Covariates and counterfactuals are balanced
    expect_equal(round(colMeans(subset(df[,nameCovCounter], df$A == 0)), 1),
                 round(colMeans(subset(df[,nameCovCounter], df$A == 1)), 1))


    ## All covariates are instruments (positive association)
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    ## Names of covariates and counter factuals
    nameCov <- Filter(f = function(elt) {grepl("^X", elt)}, x = names(df))
    ## All covariates are larger in the treated
    expect_true(all(colMeans(subset(df[,nameCov], df$A == 0)) <
                    colMeans(subset(df[,nameCov], df$A == 1))))

    ## All covariates are instruments (negative association)
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(-1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    ## All covariates are larger in the treated
    expect_true(all(colMeans(subset(df[,nameCov], df$A == 0)) >
                    colMeans(subset(df[,nameCov], df$A == 1))))

    ## All covariates are predictors (positive association)
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(0,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(1,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    ## All covariates are still balanced
    expect_equal(round(colMeans(subset(df[,nameCov], df$A == 0)), 1),
                 round(colMeans(subset(df[,nameCov], df$A == 1)), 1))

    ## All covariates are predictors (positive association)
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(0,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(1,7),
                                      betaA = 0,
                                      betaXA = rep(0,7)),
                            survParams = c(10, 10, 1))
    ## All covariates are still balanced
    expect_equal(round(colMeans(subset(df[,nameCov], df$A == 0)), 1),
                 round(colMeans(subset(df[,nameCov], df$A == 1)), 1))

    ## All covariates are confounders (positive association); treatment protective
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(1,7),
                                      betaA = -3,
                                      betaXA = rep(0,7)),
                            survParams = c(1, 1, 1))
    ## Mean treated counterfactuals are better than mean untreated counterfactuals
    with(df, {
        expect_true(mean(pY0) > mean(pY1))
        expect_true(mean(rate0) > mean(rate1))
    })
    ## Untreated counterfactual matches observed in the untreated
    with(subset(df, A == 0), {
        expect_equal(pY, pY0)
        expect_equal(rate, rate0)
    })
    ## Treated counterfactual matches observed in the treated
    with(subset(df, A == 1), {
        expect_equal(pY, pY1)
        expect_equal(rate, rate1)
    })


    ## All covariates are confounders (positive association); treatment harmful
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(1,7),
                                      betaA = +3,
                                      betaXA = rep(0,7)),
                            survParams = c(1, 1, 1))
    ## Mean treated counterfactuals are worse than mean untreated counterfactuals
    with(df, {
        expect_true(mean(pY0) < mean(pY1))
        expect_true(mean(rate0) < mean(rate1))
    })
    ## Untreated counterfactual matches observed in the untreated
    with(subset(df, A == 0), {
        expect_equal(pY, pY0)
        expect_equal(rate, rate0)
    })
    ## Treated counterfactual matches observed in the treated
    with(subset(df, A == 1), {
        expect_equal(pY, pY1)
        expect_equal(rate, rate1)
    })

    ## All covariates are instruments (positive assoc.); More protective if X5 = 1
    ##
    ## Covariates were made instruments to avoid differential baseline probabilities
    ## of disease (by X5) that affects assessment of effect heterogeneity on the
    ## bounded natural scale ([0,1] for probability; [0,Inf] for rate).
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = -3,
                                      betaXA = c(0,0,0,0,-1,0,0)),
                            survParams = c(1, 1, 1))
    ## Mean counterfactual effect is larger in X5 == 1
    expect_true(with(subset(df, X5 == 1), abs(mean(pY1) - mean(pY0))) >
                with(subset(df, X5 == 0), abs(mean(pY1) - mean(pY0))))
    expect_true(with(subset(df, X5 == 1), abs(mean(rate1) - mean(rate0))) >
                with(subset(df, X5 == 0), abs(mean(rate1) - mean(rate0))))

    ## All covariates are instruments (positive assoc.); Less protective if X5 = 1
    ##
    ## Covariates were made instruments to avoid differential baseline probabilities
    ## of disease (by X5) that affects assessment of effect heterogeneity on the
    ## bounded natural scale ([0,1] for probability; [0,Inf] for rate).
    df <- GenerateOneCenter(n = 10^5,
                            alphas = c(alpha0 = 0,
                                       alphaX = rep(1,7)),
                            betas = c(beta0 = 0,
                                      betaX = rep(0,7),
                                      betaA = -3,
                                      betaXA = c(0,0,0,0,+1,0,0)),
                            survParams = c(1, 1, 1))
    ## Mean counterfactual effect is smaller in X5 == 1
    expect_true(with(subset(df, X5 == 1), abs(mean(pY1) - mean(pY0))) <
                with(subset(df, X5 == 0), abs(mean(pY1) - mean(pY0))))
    expect_true(with(subset(df, X5 == 1), abs(mean(rate1) - mean(rate0))) <
                with(subset(df, X5 == 0), abs(mean(rate1) - mean(rate0))))

})


###
### Multiple site generation
################################################################################

test_that("Multiple centers are generated correctly", {

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

    lstDf <- GenerateDistResNet(lstN          = rep(list(n), K),
                                lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                                lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
                                lstSurvParams = rep(list(survParams), K))

    expect_true(length(lstDf) == K)
    expect_true(all(sapply(lstDf, function(site) {
        "ResSite" %in% class(site)
    })))
    expect_true(all(sapply(lstDf, function(site) {
        "data.frame" %in% class(site)
    })))
    ## Scenario parameters are kept as an attribute
    expect_equal(attr(lstDf, "ScenarioDistResNet"),
                 GenerateScenarioDistResNet(
                     lstN          = rep(list(n), K),
                     lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                     lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
                     lstSurvParams = rep(list(survParams), K)))

    ## print method
    print(lstDf)
    expect_output(print(lstDf), "Distributed research network with 3 sites")
    expect_output(print(lstDf), "Name: site3")
    expect_output(print(lstDf), "Size: 1000")

    ## summary method
    lstParamsTableOne <- summary(lstDf)
    expect_true(all("ParamsTableOne" %in% sapply(lstParamsTableOne, class)))
    expect_true(is.null(lstParamsTableOne[[1]]$params))
    ##  with parameters
    lstParamsTableOne2 <- summary(lstDf, truth = TRUE)
    expect_true(!is.null(lstParamsTableOne2[[1]]$params))
    expect_true(all(sapply(lapply(lstParamsTableOne2, "[[", "params"), length) == 4))

    ## If not all arguments have the same length, an error arises
    expect_error(GenerateDistResNet(lstN          = rep(list(n), K),
                                    lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                                    lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K-1),
                                    lstSurvParams = rep(list(survParams), K)))
    ## Inconsistent parameter length gives an error
    expect_error(GenerateDistResNet(lstN          = rep(list(n), K),
                                    lstAlphas     = rep(list(c(alpha0, alphaX[-1])), K),
                                    lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
                                    lstSurvParams = rep(list(survParams), K)))
    expect_error(GenerateDistResNet(lstN          = rep(list(n), K),
                                    lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                                    lstBetas      = rep(list(c(beta0, betaX[-1], betaA, betaXA)), K),
                                    lstSurvParams = rep(list(survParams), K)))
    expect_error(GenerateDistResNet(lstN          = rep(list(n), K),
                                    lstAlphas     = rep(list(c(alpha0, alphaX)), K),
                                    lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA[-1])), K),
                                    lstSurvParams = rep(list(survParams), K)))
})


###
### Generate scenarios (parameter sets)
################################################################################

test_that("a single scenario dataset is generated correctly", {

    ScenarioDistResNet <-
        GenerateScenarioDistResNet(list(1, 2, 3, 4),
                                   list(c(0.1, -seq_len(7)),
                                        c(0.2, -seq_len(7)),
                                        c(0.3, -seq_len(7)),
                                        c(0.4, -seq_len(7))),
                                   list(c(0.1, seq_len(7), 100, 10 + seq_len(7)),
                                        c(0.2, seq_len(7), 100, 10 + seq_len(7)),
                                        c(0.3, seq_len(7), 100, 10 + seq_len(7)),
                                        c(0.4, seq_len(7), 100, 10 + seq_len(7))),
                                   list(c(1, 1, 1),
                                        c(1, 1, 1),
                                        c(1, 1, 1),
                                        c(1, 1, 1)))
    expect_true("ScenarioDistResNet" %in% class(ScenarioDistResNet))

    ## Test print method
    expect_output(print(ScenarioDistResNet),
                  "0.1 ; -1 -2 -3 -4 -5 -6 -7")
    expect_output(print(ScenarioDistResNet),
                  "0.1 ; 1 2 3 4 5 6 7 ; 100 ; 11 12 13 14 15 16 17")


    ## Error detections
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

    ## If not all arguments have the same length, an error arises
    expect_error(GenerateScenarioDistResNet(
        lstN          = rep(list(n), K),
        lstAlphas     = rep(list(c(alpha0, alphaX)), K),
        lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K-1),
        lstSurvParams = rep(list(survParams), K)))
    ## Inconsistent parameter length gives an error
    expect_error(GenerateScenarioDistResNet(
        lstN          = rep(list(n), K),
        lstAlphas     = rep(list(c(alpha0, alphaX[-1])), K),
        lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA)), K),
        lstSurvParams = rep(list(survParams), K)))
    expect_error(GenerateScenarioDistResNet(
        lstN          = rep(list(n), K),
        lstAlphas     = rep(list(c(alpha0, alphaX)), K),
        lstBetas      = rep(list(c(beta0, betaX[-1], betaA, betaXA)), K),
        lstSurvParams = rep(list(survParams), K)))
    expect_error(GenerateScenarioDistResNet(
        lstN          = rep(list(n), K),
        lstAlphas     = rep(list(c(alpha0, alphaX)), K),
        lstBetas      = rep(list(c(beta0, betaX, betaA, betaXA[-1])), K),
        lstSurvParams = rep(list(survParams), K)))

})

## Create size set
## List of K scalars indicating site sample sizes
lstNBase <- list(10^5,2*10^4,2*10^4,5*10^3)

## Create treatment model coefficient sets
## List of K vectors
lstAlphasBase <- list(c(0.1, rep(0,7)),
                      c(0.2, rep(0,7)),
                      c(0.3, rep(0,7)),
                      c(0.4, rep(0,7)))
lstAlphasAlt1 <- list(c(0.2, rep(0,7)),
                      c(0.3, rep(0,7)),
                      c(0.4, rep(0,7)),
                      c(0.5, rep(0,7)))

## Create outcome model coefficient sets
## Intercept, covariates, treatment, interaction
## List of K vectors
lstBetasBase <- list(c(0.1, rep(0,7), 0, rep(0,7)),
                     c(0.2, rep(0,7), 0, rep(0,7)),
                     c(0.3, rep(0,7), 0, rep(0,7)),
                     c(0.4, rep(0,7), 0, rep(0,7)))
lstBetasAlt1 <- list(c(0.2, rep(0,7), 0, rep(0,7)),
                     c(0.3, rep(0,7), 0, rep(0,7)),
                     c(0.4, rep(0,7), 0, rep(0,7)),
                     c(0.5, rep(0,7), 0, rep(0,7)))

## Create survival outcome model parameter sets
## List of K vectors
lstSurvParamsBase <- list(c(1, 1, 1),
                          c(1, 1, 1),
                          c(1, 1, 1),
                          c(1, 1, 1))
lstSurvParamsAlt1 <- list(c(1/2, 1/2, 1),
                          c(1/2, 1/2, 1),
                          c(1/2, 1/2, 1),
                          c(1/2, 1/2, 1))

## Actually generate scenarios (mixing)
Scenarios <- GenerateScenarios(lstLstN          = list(lstNBase),
                               lstLstAlphas     = list(lstAlphasBase,
                                                       lstAlphasAlt1),
                               lstLstBetas      = list(lstBetasBase,
                                                       lstBetasAlt1),
                               lstLstSurvParams = list(lstSurvParamsBase,
                                                       lstSurvParamsAlt1),
                               mix = TRUE)
cat(length(Scenarios), "scenarios\n")
print(Scenarios)

## Actually generate scenarios (no mixing)
ScenariosNoMix <- GenerateScenarios(lstLstN          = list(lstNBase,
                                                            lstNBase),
                                    lstLstAlphas     = list(lstAlphasBase,
                                                            lstAlphasAlt1),
                                    lstLstBetas      = list(lstBetasBase,
                                                            lstBetasAlt1),
                                    lstLstSurvParams = list(lstSurvParamsBase,
                                                            lstSurvParamsAlt1),
                                    mix = FALSE)


test_that("Scenarios are generated correctly with all possible combinations", {

    ## Parameter set mixing scenarios
    ## Correct class as a parameter set for a DRN
    expect_true("Scenarios" %in% class(Scenarios))
    ## Each scenario is a scenario object for a DistResNet
    for (Scenario in Scenarios) {
        expect_true("ScenarioDistResNet" %in% class(Scenario))
    }

    ## Number of scenarios should be
    expect_true(length(Scenarios) == 1*2*2*2)

    ## 4 list elements (N, Alphas, Betas, SurvParams) for each scenario
    expect_true(unique(sapply(Scenarios, length)) == 4)

    ## All scenarios here have four centers
    expect_true(all(sapply(Scenarios, sapply,
                           length) == 4))

})


test_that("Scenarios are generated correctly without combinations", {

    ## Parameter set non-mixing scenarios
    ## Correct class as a parameter set for a DRN
    expect_true("Scenarios" %in% class(ScenariosNoMix))
    ## Each scenario is a scenario object for a DistResNet
    for (Scenario in ScenariosNoMix) {
        expect_true("ScenarioDistResNet" %in% class(Scenario))
    }

    ## Number of scenarios should be
    expect_true(length(ScenariosNoMix) == 2)

    ## 4 list elements (N, Alphas, Betas, SurvParams) for each scenario
    expect_true(unique(sapply(ScenariosNoMix, length)) == 4)

    ## All scenarios here have four centers
    expect_true(all(sapply(ScenariosNoMix, sapply,
                           length) == 4))

    ## This equality should hold
    scenario1 <- GenerateScenarioDistResNet(lstNBase,
                                            lstAlphasBase,
                                            lstBetasBase,
                                            lstSurvParamsBase)
    scenario2 <- GenerateScenarioDistResNet(lstNBase,
                                            lstAlphasAlt1,
                                            lstBetasAlt1,
                                            lstSurvParamsAlt1)
    ansScenariosNoMix <- list(scenario1, scenario2)
    class(ansScenariosNoMix) <- c("Scenarios", class(ansScenariosNoMix))
    expect_equal(ScenariosNoMix,
                 ansScenariosNoMix)
})


## Function to test data generation for one scenario
TestOneScenarioGen <- function(scenarioParams, lstNBase) {
    ## Realize data for the first scenario
    DistResNet <- GenerateDistResNet(lstN          = scenarioParams[[1]],
                                     lstAlphas     = scenarioParams[[2]],
                                     lstBetas      = scenarioParams[[3]],
                                     lstSurvParams = scenarioParams[[4]])
    ## DistResNet containing ResSite
    expect_true("DistResNet" %in% class(DistResNet))
    expect_true(all("ResSite" %in% sapply(DistResNet, class)))
    ## Really a list of data frames
    expect_true("list" %in% class(DistResNet))
    expect_true(all("data.frame" %in% sapply(DistResNet, class)))

    ## There are four centers.
    expect_true(length(DistResNet) == 4)

    ## Check sample sizes match the prespecified sizes
    expect_true(all(mapply(FUN = function(a,b) {
        a == nrow(b)
    },
    lstNBase,
    DistResNet)))
}


test_that("Data are generated correctly for multiple scenarios", {

    ## Loop over all scenarios, and check each scenario data
    lapply(Scenarios,
           TestOneScenarioGen, lstNBase = lstNBase)

})
