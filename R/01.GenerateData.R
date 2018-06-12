################################################################################
### Simulation data generation functions
##
## Created on: 2016-04-28
## Author: Kazuki Yoshida
################################################################################


###
### Covariate creation
################################################################################
## Independent covariates.
## First one is a continuous and modeled after an age variable.
## The others emulate baseline comorbidities of varying frequencies.

##' Create indepedent covariates given mean, sd, and prevalence
##'
##' The first covariate is a Normal(mean1, sd1) random variable, whereas the other
##' covariates are Bernoulli(p) random variables where p's are specified by the
##' elements of the pp vector.
##'
##' @param n number of individuals to create
##' @param mean1 mean for first continuous normal covariate
##' @param sd1 sd for first continuous normal covariate
##' @param pp vector of prevalences for other Bernoulli covariates
##'
##' @return n x (m+1) data frame including ID and m covariates
##'
##' @author Kazuki Yoshida
##'
##' @export
AssignCovariatesNormBin <- function(n, mean1, sd1, pp) {
    ## ID
    id <- seq_len(n)
    ## One continuous covariate
    X1 <- rnorm(n = n, mean = mean1, sd = sd1)
    out <- data.frame(id = id,
                      X1 = X1)
    ## binary covariates (start indexing at 2)
    Xs <- lapply(pp, function(p) {
        rbinom(n = n, size = 1, prob = p)
    })
    names(Xs) <- paste0("X", seq_along(pp) + 1)
    ## Return as a data frame
    cbind(out, do.call(cbind, Xs))
}

##' Create independent covariates by default method
##'
##' The first covariate is a Normal(0,1) random variable, whereas the other covariates are Bernoulli(p) random variables where p's are specified by the elements of the pp vector.
##'
##' @param n number of individuals to create
##' @param p number of covariates
##'
##' @return n x (p+1) data frame including ID and p covariates
##'
##' @author Kazuki Yoshida
##'
##' @export
AssignCovariatesNormBinDefault <- function(n, p) {
    ## Probabilities for binary covariates (1st is continuous, thus, 1 fewer)
    ## Range probabilities from 0.05 to 0.50
    mean1 <- 0
    sd1 <- 1
    pp <- seq(from = 0.05, to = 0.50, length.out = (p - 1))
    AssignCovariatesNormBin(n, mean1 = mean1, sd1 = sd1, pp = pp)
}

##' Construct independent covariate generator
##'
##' Construct an independent covariate generator that generates the first covariate as Normal(mean1,sd1) random variable, whereas the other covariates are Bernoulli(p) random variables where p's are specified by equal interval sequence of (p2, ..., pp)
##'
##' @param mean1 mean for first continuous normal covariate
##' @param sd1 sd for first continuous normal covariate
##' @param p2 prevalence of binary X2
##' @param pLast prevalence of binary Xlast
##'
##' @return covariate generator that takes n and p as arguments
##'
##' @author Kazuki Yoshida
##'
##' @export
ConstructAssignCovariatesNormBin <- function(mean1, sd1, p2, pLast) {
    function(n, p) {
        ## Probabilities for binary covariates (1st is continuous, thus, 1 fewer)
        pVector <- seq(from = p2, to = pLast, length.out = (p - 1))
        AssignCovariatesNormBin(n, mean1 = mean1, sd1 = sd1, pp = pVector)
    }
}


###
### Treatment assignment
################################################################################
## Assign treatment using a logistic regression model.
## True propensity scores are recorded.

##' Assign treatment based on true propensity score model
##'
##' A binary treatment is assigned as a Bernoulli(p_i) random variable where
##' p_i is an individual-specific treatment probability (true propensity score)
##' calculated from covariate values and coefficients.
##'
##' @param dfX data frame including covariates
##' @param alpha0 a scalar value coeffcient for the intercept
##' @param alphaX a vector of length 7 specifiying true coefficients for 7 covariates
##'
##' @return data frame with additional true propensity score pA and treatment A.
##'
##' @author Kazuki Yoshida
##'
##' @export
AssignTreatment <- function(dfX, alpha0, alphaX) {
    ## Sanity check
    assert_that(class(dfX) == "data.frame")
    assert_that(length(alpha0) == 1)
    assert_that(length(alphaX) == length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(dfX))))
    assert_that(names(dfX)[1] == "id")

    ## Log odds of treatment
    lpA <- as.numeric(alpha0 + as.matrix(dfX[,-1]) %*% alphaX)
    ## True PS
    dfX$pA <- expit(lpA)
    ## Treatment assignment
    dfX$A <- rbinom(n = length(dfX$pA), size = 1, prob = dfX$pA)
    ## return with additional variables
    dfX
}


###
### Outcome assignment
################################################################################

###  Binary outcome assignment
##' Assign binary outcome based on covariates and treatment
##'
##' A binary outcome is assigned as a Bernoulli(p_i) random variable where
##' p_i is an individual-specific outcome probability calculated from covariate
##' and treatment values and coefficients.
##'
##' @param dfXA data frame including covariates and treatment assingment A
##' @param beta0 scalar value coeffcient for the intercept
##' @param betaX vector specifying true coefficients for covariates
##' @param betaA scalar value for treatment effect
##' @param betaXA vector specifying true coefficients for covariates-treatment interaction. The vector must be the same length as \code{betaX}.
##'
##' @return data frame with additional counterfactual and observed variables related to the outcome.
##'
##' @author Kazuki Yoshida
##'
##' @export
AssignOutcomeBin <- function(dfXA, beta0, betaX, betaA, betaXA) {
    ## Sanity check
    assert_that(class(dfXA) == "data.frame")
    assert_that(length(beta0) == 1)
    assert_that(length(betaX) == length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(dfXA))))
    assert_that(length(betaA) == 1)
    assert_that(length(betaXA) == length(betaX))
    assert_that(names(dfXA)[1] == "id")

    ## Covariate names
    covs <- paste0("X", seq_len(length(betaX)))

    ## Log odds of disease
    lpY  <- as.numeric(beta0 + (as.matrix(dfXA[,covs]) %*% betaX) + (betaA * dfXA$A) + ((as.matrix(dfXA[,covs]) %*% betaXA) * matrix(dfXA$A)))
    ## Counterfactuals
    lpY0 <- as.numeric(beta0 + (as.matrix(dfXA[,covs]) %*% betaX) + (betaA * 0)      + ((as.matrix(dfXA[,covs]) %*% betaXA) * 0))
    lpY1 <- as.numeric(beta0 + (as.matrix(dfXA[,covs]) %*% betaX) + (betaA * 1)      + ((as.matrix(dfXA[,covs]) %*% betaXA) * 1))

    ## Probability of disease
    dfXA$pY <- expit(lpY)
    ## Counterfactuals
    dfXA$pY0 <- expit(lpY0)
    dfXA$pY1 <- expit(lpY1)

    ## Disease assignment
    dfXA$Y <- rbinom(n = length(dfXA$pY), size = 1, prob = dfXA$pY)
    ## return with additional variables
    dfXA
}


###  Survival outcome assignment
##' Assign survival time outcome
##'
##' A time to event outcome is assigned as an Exponential(lambda_i) random
##' variable where lambda_i is an individual-specific time-constant hazard
##' of an event calculated from covariate and treatment values and coefficients.
##' Censoring times are generated from an Exponential(lambda_c_i) distribution
##' where lambda_c_i is an individual-specific time-constant hazard of censoring,
##' depending on the covariates only. The censoring time distribution does not
##' depend on the treatment to avoid collider stratification bias. An additional
##' administrative censoring time is defined by Tmax, which marks the end of the
##' follow up for those who have follow up times longer than Tmax.
##' The latent true event time and censoring time variables are coded in years.
##' The final observed time variable is coded in days (\code{ceiling} is applied
##' to avoid fractional numbers) to emulate the granularity of typical claims data.
##'
##' @param dfXA data frame including covariates and treatment assingment A
##' @param betaX vector specifying true coefficients for covariates
##' @param betaA scalar value for treatment effect
##' @param betaXA vector specifying true coefficients for covariates-treatment interaction. The vector must be the same length as \code{betaX}.
##' @param lambda scalar value defining the constant baseline hazard for the event outcome. The corresponding time scale is in years.
##' @param lambda_c scalar value defining the constant baseline hazard for censoring. The corresponding time scale is in years.
##' @param Tmax scalar value defining the end of follow up. Use to introduce a constant value administrative censoring time. The corresponding time scale is in years. To avoid administrative censoring, set this to a very high value.
##'
##' @return data frame with added individual-specific rate (\code{rate}), counterfactual rates (\code{rate0} and \code{rate1}), event time (\code{T}; years), censoring time (\code{C}; years), observed time (\code{time}; integer days), and event indicator (\code{event}).
##'
##' @author Kazuki Yoshida
##'
##' @export
AssignOutcomeSurv <- function(dfXA, betaX, betaA, betaXA, lambda, lambda_c, Tmax) {
    assert_that(class(dfXA) == "data.frame")
    assert_that(length(betaX) == length(Filter(f = function(elt) {grepl("^X", elt)}, x = names(dfXA))))
    assert_that(length(betaA) == 1)
    assert_that(length(betaXA) == length(betaX))
    assert_that(length(lambda) == 1)
    assert_that(length(lambda_c) == 1)
    assert_that(length(Tmax) == 1)
    assert_that(names(dfXA)[1] == "id")

    ## Number of observations to generate
    n <- nrow(dfXA)

    ## Covariate names
    covs <- paste0("X", seq_len(length(betaX)))

    ## Event
    ## Linear predictor for individual-specific log hazard ratios
    ## This corresponds to predict(coxph_model, type = "lp").
    lp  <- as.numeric((as.matrix(dfXA[,covs]) %*% betaX) + (betaA * dfXA$A) + ((as.matrix(dfXA[,covs]) %*% betaXA) * matrix(dfXA$A)))
    ## Counterfactuals
    lp0 <- as.numeric((as.matrix(dfXA[,covs]) %*% betaX) + (betaA * 0)      + ((as.matrix(dfXA[,covs]) %*% betaXA) * 0))
    lp1 <- as.numeric((as.matrix(dfXA[,covs]) %*% betaX) + (betaA * 1)      + ((as.matrix(dfXA[,covs]) %*% betaXA) * 1))
    ## Individual-specific rate: lambda_i = lambda * exp(lp)
    dfXA$rate  <- lambda * exp(lp)
    ## Counterfactuals
    dfXA$rate0 <- lambda * exp(lp0)
    dfXA$rate1 <- lambda * exp(lp1)
    ## Assign event time from exponential distribution
    dfXA$T <- rexp(n = n, rate = dfXA$rate)

    ## Censoring
    ## Censoring time based on the individual-specific censoring rates
    ## These are still non-informative conditional on the covariates.
    ## There is no dependency on the treatment.
    ## Censoring probability = lambda_c / (lambda + lambda_c) if not treated.
    dfXA$C <- rexp(n = n,
                   rate = as.numeric(lambda_c * exp(as.matrix(dfXA[,covs]) %*% betaX)))

    ## Observed time and censoring indicator
    dfXA$time  <- pmin(dfXA$C, dfXA$T, Tmax)
    ## Event happened if the latent event time is the observed time.
    dfXA$event <- as.numeric(dfXA$time == dfXA$T)

    ## Make discrete day level data
    dfXA$time <- ceiling(dfXA$time * 365.25)

    ## return
    dfXA
}


###
### Data center creation
################################################################################

## Helper function for GenerateOneCenter
## Should only be used within GenerateOneCenter
DissectParams <- function(alphas, betas, survParams) {
    ## Number of covariates
    p <- length(alphas) - 1
    ## Extract coefficients for readability
    ## First element is regarded as intercept
    alpha0  <- alphas[1]
    alphaX  <- alphas[-1]
    ## First element is the intercept
    beta0   <- betas[1]
    ## 2nd through (p + 1)-th elements are the main effects of X
    betaX   <- betas[seq(from = 2, to = (p + 1))]
    ## (p + 2)-th element is the main effect of treatment
    betaA   <- betas[(p + 2)]
    ## (p + 3)-th through (2 * p + 2)-th elements are the interaction effects
    betaXA  <- betas[seq(from = (p + 3), to = (2 * p + 2))]
    ## First is baseline rate, second is maximum observation time
    lambda   <- survParams[1]
    lambda_c <- survParams[2]
    Tmax     <- survParams[3]
    ## Return as a list
    list(alpha0   = alpha0,
         alphaX   = alphaX,
         beta0    = beta0,
         betaX    = betaX,
         betaA    = betaA,
         betaXA   = betaXA,
         lambda   = lambda,
         lambda_c = lambda_c,
         Tmax     = Tmax)
}


##' Generate one data center
##'
##' Generate one data center with information of n patients. The parameters are given as four raw vectors. This function is intentionally designed to be a low level function working with raw numeric vector parameters for maximum flexibility. In actual simulation runs, this function should be called by \code{\link{GenerateDistResNet}}, not directly.
##'
##' @param n study site-specific sample size
##' @param AssignCovariates covariate generation functions that takes n and p as the only arguments.
##' @param alphas parameter vector for treatment model including c(alpha0, alphaX)
##' @param betas parameter vector for outcome model shared among binary and survival outcome models including \code{c(beta0, betaX, betaA, betaXA)}.
##' @param survParams vector of two. The first element is the baseline hazard of events in the exponential event time outcome model (\code{lambda}). The second element is the baseline hazard of censoring in the exponential censoring time model (\code{lambda_c}).
##'
##' @return data frame with id, covariates, treatment (A) related variables, binary outcome (Y) related variables, and survival outcome (time, event) related variables.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateOneCenter <- function(n, AssignCovariates, alphas, betas, survParams) {
    ## All arguments are numeric
    assert_that(is.numeric(n) & is.numeric(alphas) & is.numeric(betas) & is.numeric(survParams))
    ## Sample size is a scalar
    assert_that(length(n) == 1)
    ## Sample size is positive
    assert_that(n > 0)
    ## Sample size is not a fractional number
    assert_that(round(n) == n)
    ## Number of covariates in the X matrix (excluding intercept)
    p <- length(alphas) - 1
    ## intercept and covariates. There must be at least one continuous and binary.
    assert_that(p >= 2)
    ## intercept, covariates, treatment, covariate*treatment coefficients
    assert_that(length(betas) == (2 * p + 2))
    ## lambda and lambda_c
    assert_that(length(survParams) == 3)

    ## Extract coefficients for readability
    params <- DissectParams(alphas = alphas, betas = betas, survParams = survParams)
    ## Sanity check for the number of covariates
    assert_that(length(params$alphaX) == length(params$betaX))

    ## Generate dataset
    ## Covariates
    out <- AssignCovariates(n = n, p = p) %>%
        ## Treatment A
        AssignTreatment(dfX        = .,
                        alpha0     = params$alpha0,
                        alphaX     = params$alphaX) %>%
        ## Binary outcome
        AssignOutcomeBin(dfXA      = .,
                         beta0     = params$beta0,
                         betaX     = params$betaX,
                         betaA     = params$betaA,
                         betaXA    = params$betaXA) %>%
        ## Survival outcome
        AssignOutcomeSurv(dfXA     = .,
                          betaX    = params$betaX,
                          betaA    = params$betaA,
                          betaXA   = params$betaXA,
                          lambda   = params$lambda,
                          lambda_c = params$lambda_c,
                          Tmax     = params$Tmax)

    ## Check and return
    assert_that(class(out) == "data.frame")
    assert_that(all(sapply(out, class) %in% c("integer","numeric")))
    assert_that(nrow(out) == n)
    assert_that(names(out)[1] == "id")
    ## Keep parameters as attributes for peeking the truth
    attr(out, which = "ScenarioResSite") <-
        list(n = n,
             AssignCovariates = AssignCovariates,
             alphas = alphas,
             betas = betas,
             survParams = survParams)
    ## Add class Research Site
    class(out) <- c("ResSite", class(out))
    out
}


##' print method for a ResSite class object
##'
##' Print the sample size and analysis readiness of a research site in a distributed research network.
##'
##' @param x ResSite object to print
##' @param ... for compatibility with generic functions.
##'
##' @return NULL
##'
##' @author Kazuki Yoshida
##'
##' @export
print.ResSite <- function(x, ...) {
    cat("### Size:", nrow(x), "\n")
    if ("ResSiteReady" %in% class(x)) {
        cat("### Site ready for analysis\n")
    } else {
        cat("### Site not ready for analysis\n")
    }
}


##' summary method for a ResSite class object
##'
##' Summarizes variables in a research site in a distributed research network.
##'
##' @inheritParams print.ResSite
##' @param truth whether to include the unobserved latent variables including counterfactual probabilities of the binary outcome assignment, counterfactual rates of the survival events, as well as true underlying 1-year survival \code{S1yr}.
##' @param ... for compatibility with generic functions.
##'
##' @return \code{ParamsTableOne} object containing parameter values (if \code{truth = TRUE}) and TableOne object defined in \code{tableone}. If \code{truth} is TRUE, then it will include parameters, latent variables, and the true one-year survival probability.
##'
##' @author Kazuki Yoshida
##'
##' @export
summary.ResSite <- function(x, truth = FALSE, ...) {

    ## Calculate true 1 year survival
    x$S1yr <- as.numeric(x$T >= 1)

    ## Pick variables to include in table depending on whether truth is requested.
    if (truth) {
        ## Show all variables except for A (stratifying variable)
        vars <- names(x)[!(names(x) %in% c("id"))]
    } else {
        ## Just show observed variables
        vars <- Filter(f = function(elt) {grepl("^X|^Y|^A|^time|^event", elt)}, x = names(x))
    }

    ## Create TableOne object avoid testing to save time
    tab    <- tableone::CreateTableOne(vars = vars, data = x, test = FALSE)
    tabByA <- tableone::CreateTableOne(vars = vars, strata = "A", data = x, test = FALSE)
    ## Combine in a list
    out <- list(TableOneOverall = tab,
                TableOneByA     = tabByA)
    ## Do not include parameters if latent variables are not included
    if (truth) {
        out$params <- attr(x, which = "ScenarioResSite")
    } else {
        out$params <- NULL
    }

    class(out) <- "ParamsTableOne"
    out
}


##' Generate one data center and its table
##'
##' Generate one data center with information of n patients, and then generate a baseline characteristics table using the \code{tableone} package.
##'
##' @inheritParams GenerateOneCenter
##' @param ... parameters for tableone::CreateTableOne function
##'
##' @return ParamsTableOne object containing parameter values and TableOne object defined in \code{tableone}. The true one-year survival variable is added to the variables already in the generated dataset.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateOneCenterTable <- function(n, AssignCovariates, alphas, betas, survParams, ...) {
    ## Generate one center
    ResSite <- GenerateOneCenter(n = n, AssignCovariates = AssignCovariates, alphas = alphas, betas = betas, survParams = survParams)

    ## Just call summary method for a ParamsTableOne object
    summary(ResSite, truth = TRUE)
}


##' print method for the ParamsTableOne object
##'
##' Print the table summarizing variables. If included, parameter values are also printed. Invisibly return a matrix containing the table.
##'
##' @param x ParamsTableOne object generated by \code{\link{summary.ResSite}} or \code{\link{GenerateOneCenterTable}}.
##' @param print whether to print results
##' @param ... additional parameters passed to \code{print.TableOne}
##'
##' @return return a matrix object. Print results.
##'
##' @author Kazuki Yoshida
##'
##' @export
print.ParamsTableOne <- function(x, print = TRUE, ...) {

    ## Generate overall and table stratified by A
    matOverall <- print(x$TableOneOverall, printToggle = FALSE, ...)
    matByA     <- print(x$TableOneByA, printToggle = FALSE, smd = TRUE, ...)

    ## Combine
    out <- cbind(matOverall,
                 matByA)

    if (print) {
        ## If parameters are included show them
        if (!is.null(x$params)) {
            ## Number of covariates in the X matrix (excluding intercept)
            params <- DissectParams(alphas     = x$params$alphas,
                                    betas      = x$params$betas,
                                    survParams = x$params$survParams)
            ## Show parameter values
            cat("\n### Parameter values\n")
            cat("### n        :", x$params$n, "\n")
            cat("### alpha0   :", params$alpha0, "\n")
            cat("### alphaX   :", params$alphaX, "\n")
            cat("### beta0    :", params$beta0, "\n")
            cat("### betaX    :", params$betaX, "\n")
            cat("### betaA    :", params$betaA, "\n")
            cat("### betaXA   :", params$betaXA, "\n")
            cat("### lambda   :", params$lambda, "\n")
            cat("### lambda_c :", params$lambda_c, "\n")
            cat("### Tmax     :", params$Tmax, "\n")
        }
        ## Show the table
        cat("\n### Baseline table\n")
        print(out, quote = FALSE)
    }
    ## Always return the matrix object
    invisible(out)
}


###
### Distributed research network creation
################################################################################

##' Generate a distributed research network.
##'
##' Generate a distributed research network of K sites.
##'
##' @param lstN list of sizes of centers. Number of centers (K) is also determined by its length.
##' @param lstAssignCovariates list of covariate generation functions that takes n and p as the only arguments.
##' @param lstAlphas list of alpha vectors potentially varying among centers. See \code{\link{GenerateOneCenter}} for construction of an alpha vector.
##' @param lstBetas list of beta vectors potentially varying among centers. See \code{\link{GenerateOneCenter}} for construction of a beta vector.
##' @param lstSurvParams list of survival outcome parameters potentially varying among centers.
##'
##' @return list of K data frames, each of which contains a complete dataset from a single data center
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateDistResNet <- function(lstN, lstAssignCovariates, lstAlphas, lstBetas, lstSurvParams) {

    ## Check for same lengths
    assert_that(all(length(lstN) == sapply(list(lstAssignCovariates, lstAlphas, lstBetas, lstSurvParams), length)))

    out <- mapply(FUN = GenerateOneCenter,
                  ## vector or list of parameters
                  ## i-th elements from all of them are fed to FUN
                  lstN,
                  lstAssignCovariates,
                  lstAlphas,
                  lstBetas,
                  lstSurvParams,
                  ## Return as a list
                  SIMPLIFY = FALSE)

    ## Keep parameters as an attribute
    attr(out, which = "ScenarioDistResNet") <-
        GenerateScenarioDistResNet(lstN = lstN,
                                   lstAssignCovariates = lstAssignCovariates,
                                   lstAlphas = lstAlphas,
                                   lstBetas = lstBetas,
                                   lstSurvParams = lstSurvParams)
    ## Give site names
    names(out) <- paste0("site", seq_along(out))
    ## Add class Distributed Research Network
    class(out) <- c("DistResNet", class(out))
    out
}


##' print method for a DistResNet class object
##'
##' Run the print method for each center (ResSite) and print results.
##'
##' @param x DistResNet object
##' @param ... for compatibility with generic functions.
##'
##' @return NULL. Print results
##'
##' @author Kazuki Yoshida
##'
##' @export
print.DistResNet <- function(x, ...) {
    cat("### Distributed research network with", length(x), "sites \n")
    cat("\n")
    ## Index over sites
    ii <- seq_along(x)
    for (i in ii) {
        cat("### Name:", names(x)[[i]], " \n")
        ## invoke print.ResSite
        print(x[[i]])
        cat("\n")
    }
}


##' summary method for a DistResNet class object
##'
##' Summarizes variables for each research site in a distributed research network.
##'
##' @param x DistResNet object
##' @inheritParams summary.ResSite
##' @param ... for compatibility with generic functions.
##'
##' @return a list of \code{ParamsTableOne} objects.
##'
##' @author Kazuki Yoshida
##'
##' @export
summary.DistResNet <- function(x, truth = FALSE, ...) {
    ## Show parameters if truth is asked for.
    if (truth)  {
        print(attr(x, "ScenarioDistResNet"))
    }
    ## Just loop over sites
    lapply(x, summary, truth = truth)
}


###
### Scenario creation
################################################################################

##' Generate ScenarioDistResNet object containing parameters for a DRN.
##'
##' Generate a ScenarioDistResNet object, a set of parameters specifying a distributed research network data generation. The number of sites are implied by the length of parameters (K).
##'
##' @param lstN list of K sample sizes
##' @param lstAssignCovariates list of covariate generation functions that takes n and p as the only arguments.
##' @param lstAlphas list of K vectors, each of which specifies alpha parameters for the treatment assignment model.
##' @param lstBetas list of K vectors, each of which specifies beta parameters for the outcome assignment model.
##' @param lstSurvParams list of K vectors, each of which specifies survival parameters.
##'
##' @return ScenarioDistResNet object.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateScenarioDistResNet <- function(lstN, lstAssignCovariates, lstAlphas, lstBetas, lstSurvParams) {
    ## All arguments must be lists.
    assert_that(class(lstN) == "list")
    assert_that(class(lstAssignCovariates) == "list")
    assert_that(class(lstAlphas) == "list")
    assert_that(class(lstBetas) == "list")
    assert_that(class(lstSurvParams) == "list")

    ## They all have to have the same length
    site_count_K <- unique(sapply(list(lstN, lstAssignCovariates, lstAlphas, lstBetas, lstSurvParams), length))
    assert_that(length(site_count_K) == 1)

    ## In each site, betas are twice as long as alphas
    for (i in seq_along(lstN)) {
        assert_that(length(lstBetas[[i]]) == 2 * length(lstAlphas[[i]]))
    }

    ## Combine them in a list
    out <- list(lstN = lstN,
                lstAssignCovariates = lstAssignCovariates,
                lstAlphas = lstAlphas,
                lstBetas = lstBetas,
                lstSurvParams = lstSurvParams)

    ## Add correct class
    class(out) <- c("ScenarioDistResNet", class(out))
    out
}

##' Generate Scenarios object containing multiple DRN scenarios.
##'
##' Generate a Scenario object, a list of ScenarioDistResNet objects, each of which corresponds to one simulation scenario by combining the possible values given.
##'
##' @param lstLstN list of lists of sample sizes
##' @param lstLstAlphas list of lists of covariate generating functions
##' @param lstLstAlphas list of lists of alpha parameters
##' @param lstLstBetas list of lists of beta parameters
##' @param lstLstSurvParams list of lists of survival parameters
##' @param mix whether to mix parameter sets. If FALSE, the first parameter sets for all domains are combined as the first scenario, etc. The parameter set has to be of the same length in this case.
##'
##' @return Scenarios object, which is a list of ScenarioDistResNet objects. Each ScenarioDistResNet a parameter set for a single scenario.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateScenarios <- function(lstLstN, lstLstAssignCovariates, lstLstAlphas, lstLstBetas, lstLstSurvParams, mix) {

    ## All of arguments must be lists.
    assert_that(class(lstLstN) == "list")
    assert_that(class(lstLstAssignCovariates) == "list")
    assert_that(class(lstLstAlphas) == "list")
    assert_that(class(lstLstBetas) == "list")
    assert_that(class(lstLstSurvParams) == "list")

    ## All list elements are also lists.
    assert_that(unique(unlist(lapply(lstLstN, class))) == "list")
    assert_that(unique(unlist(lapply(lstLstAssignCovariates, class))) == "list")
    assert_that(unique(unlist(lapply(lstLstAlphas, class))) == "list")
    assert_that(unique(unlist(lapply(lstLstBetas, class))) == "list")
    assert_that(unique(unlist(lapply(lstLstSurvParams, class))) == "list")

    ## All list elements are of the same length
    lengthLstN <- unique(unlist(lapply(lstLstN, length)))
    lengthLstAlphas <- unique(unlist(lapply(lstLstAlphas, length)))
    lengthLstAssignCovariates <- unique(unlist(lapply(lstLstAssignCovariates, length)))
    lengthLstBetas <- unique(unlist(lapply(lstLstBetas, length)))
    lengthLstSurvParams <- unique(unlist(lapply(lstLstSurvParams, length)))
    assert_that(identical(lengthLstN, lengthLstAssignCovariates))
    assert_that(identical(lengthLstN, lengthLstAlphas))
    assert_that(identical(lengthLstN, lengthLstBetas))
    assert_that(identical(lengthLstN, lengthLstSurvParams))

    ## If not mixing, all parameter set arguments must be of the same lengths.
    if (mix == FALSE) {
        len_args <- sapply(list(lstLstN, lstLstAssignCovariates, lstLstAlphas, lstLstBetas, lstLstSurvParams),
                           length)
        assert_that(length(unique(len_args)) == 1)
    }

    if (mix) {
        ## If mixing, use the expand.grid technique to obtain all possible combinations.
        ## Make sure this ordring is the same as GenData arguments
        outGrid <- expand.grid(lstLstN = lstLstN,
                               lstLstAssignCovariates = lstLstAssignCovariates,
                               lstLstAlphas = lstLstAlphas,
                               lstLstBetas = lstLstBetas,
                               lstLstSurvParams = lstLstSurvParams)

        ## Create a list of ScenarioDistResNet objects
        lstScenarios <- lapply(seq_len(nrow(outGrid)), function(i) {
            ## Take each low out as a list
            lstOut <- lapply(outGrid[i,], "[[", 1)

            ## Construct a scenario for a single distributed research network
            out <- GenerateScenarioDistResNet(lstN = lstOut[[1]],
                                              lstAssignCovariates = lstOut[[2]],
                                              lstAlphas = lstOut[[3]],
                                              lstBetas = lstOut[[4]],
                                              lstSurvParams = lstOut[[5]])
            out
        })

    } else {
        ## If not mixing just take one by one
        lstScenarios <- mapply(FUN = GenerateScenarioDistResNet,
                               lstLstN,
                               lstLstAssignCovariates,
                               lstLstAlphas,
                               lstLstBetas,
                               lstLstSurvParams,
                               SIMPLIFY = FALSE)
    }

    ## Scenarios class
    class(lstScenarios) <- c("Scenarios", class(lstScenarios))
    lstScenarios
}


##' print method for a ScenarioDistResNet class object
##'
##' Print parameters for a scenario for a distributed research network.
##'
##' @param x ScenarioDistResNet object
##' @param ... for compatibility with generic functions.
##'
##' @return NULL. Print results
##'
##' @author Kazuki Yoshida
##'
##' @export
print.ScenarioDistResNet <- function(x, ...) {
    cat("### Scenario for a distributed research network.\n")
    cat("###", length(x$lstN), "sites are included.\n")

    cat("### Sample sizes\n")
    for (elt in x$lstN) {
        cat(elt, "\n")
    }

    ## Number of covariates
    nCovs <- length(x$lstAlphas[[1]]) - 1

    cat("### Alpha parameters (treatment model)\n")
    cat("### Intercept ; Covariates\n")
    for (elt in x$lstAlphas) {
        cat(elt[1],
            ";",
            elt[-1],
            "\n")
    }

    cat("### Beta parameters (outcome model)\n")
    cat("### Intercept ; Covariates ; Treatment ; Interactions\n")
    for (elt in x$lstBetas) {
        cat(elt[1],
            ";",
            elt[1 + seq_len(nCovs)],
            ";",
            elt[1 + 1 + nCovs],
            ";",
            elt[1 + 1 + nCovs + seq_len(nCovs)],
            "\n")
    }

    cat("### Survival parameters (survival outcome only)\n")
    cat("### Event rate ; Censoring rate ; Administrative censoring\n")
    for (elt in x$lstSurvParams) {
        cat(elt, "\n")
    }
    cat("\n")
}

##' print method for a Scenarios class object
##'
##' Print parameters for each scenario included in a set of scenarios (Scenarios object).
##'
##' @param x Scenarios object
##' @param ... for compatibility with generic functions.
##'
##' @return NULL. Print results
##'
##' @author Kazuki Yoshida
##'
##' @export
print.Scenarios <- function(x, ...) {
    cat("### Simulation scenarios\n")
    cat("###", length(x), "scenarios are contained.\n")
    cat("#####################################################\n\n")

    if (is.null(names(x))) {
        ## If no scenarios names are available
        ## Print each scenario
        for (i in seq_along(x)) {
            cat("### Scenario", i, "\n")
            print(x[[i]])
        }
    } else {
        ## If scenarios names are available
        ## Print each scenario
        for (i in seq_along(x)) {
            cat("### Scenario", i,
                paste0("(", names(x)[i], ")"),
                "\n")
            print(x[[i]])
        }
    }
}


###
### Data generation and saving
################################################################################

##' Generate R distributed research network iterations for a single scenario
##'
##' Given a parameter set for a single scenario, generate R iterations of distributed research network. Each distributed research network is generated by \code{\link{GenerateDistResNet}} and within-site data preparation is also conducted by \code{\link{RequestSiteDataPreparation}}.
##'
##' @param ScenarioDistResNet A ScenarioDistResNet object contained in a Scenarios object generated by \code{\link{GenerateScenarios}}.
##' @param R scalar value. Data are generated R times.
##' @param scenarioCount scalar value. Indicate the scenario count. Used for data file name.
##' @param partCount scalar value. Indicates which subpart of
##'
##' @return Use for its side effect. No return value. Save a data file in the working directory. The data file contains \code{ScenarioDistResNet}, \code{scenarioCount}, \code{R}, and \code{lstIter}. \code{lstIter} is a list of length \code{R}. Each element is a list of \code{K} centers.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateDistResNetRtimesAndSave <- function(ScenarioDistResNet, R, scenarioCount, partCount) {
    ## Sanity check
    assert_that("ScenarioDistResNet" %in% class(ScenarioDistResNet))
    assert_that(is.numeric(R))
    assert_that(length(R) == 1)
    assert_that(is.numeric(scenarioCount))
    assert_that(length(scenarioCount) == 1)
    assert_that(is.numeric(partCount))
    assert_that(length(partCount) == 1)

    ## Generate R iterations for a single scenario
    ## This is not parallelized. Parallelize the outer loop over scenarios.
    ## List of R iterations. Each iteration is a list of K centers.
    lstIter <- lapply(seq_len(R), function(i) {

        ## Generate a distributed research network object
        DistResNet <- GenerateDistResNet(lstN = ScenarioDistResNet[[1]],
                                         lstAssignCovariates = ScenarioDistResNet[[2]],
                                         lstAlphas = ScenarioDistResNet[[3]],
                                         lstBetas = ScenarioDistResNet[[4]],
                                         lstSurvParams = ScenarioDistResNet[[5]])
        DistResNet
    })

    ## Generate file name
    fileName <- sprintf("ScenarioRaw%03d_part%03d_R%d.RData",
                        scenarioCount,
                        partCount,
                        R)

    ## Save
    save(ScenarioDistResNet, R, scenarioCount, partCount, lstIter,
         file = fileName)

    ## No return value
    NULL
}


##' Generate R iterations for each scenario
##'
##' Given an iteration count and a list of scenario parameters for multiple scenarios (\code{Scenarios} object), generate R datasets for each scenario. If a scenario within a \code{Scenarios} object is corrupt (not correctly a \code{ScenarioDistResNet} object), it is skipped by \code{try}.
##'
##' @param Scenarios a Scenarios object. Each list element is a specification for a given scenario.
##' @param parts how many subfiles to create for each scenario. Use this to ease parallelization across cluster nodes.
##' @param R iteration count for each subfile. The total iteration count for a given scenario is \code{parts * R}.
##'
##' @return Use for its side effect. There is no return value. Data files are created in the same folder.
##'
##' @author Kazuki Yoshida
##'
##' @export
GenerateDataForAllScenarios <- function(Scenarios, parts, R) {
    ## Sanity check
    assert_that("Scenarios" %in% class(Scenarios))
    assert_that(is.numeric(parts))
    assert_that(length(parts) == 1)
    assert_that(is.numeric(R))
    assert_that(length(R) == 1)

    ## Message for information
    cat("### Generating data files.\n")
    cat("### Generating data files for", length(Scenarios), "scenarios.\n")
    cat("###", parts,
        "subfiles with", R,
        "iterations each for each scenario.\n")
    cat("### Total iteration count:", parts*R, "for each scenario.\n")

    ## Generate all possible (scenarioCount, partCount) tuples
    ## and put in a list
    lstTuples <-
        expand.grid(scenarioCount = seq_along(Scenarios),
                    partCount     = seq_len(parts)) %>%
        t %>%
        as.data.frame %>%
        as.list

    ## Reproducible parallel execution
    ## Parallelization is at the file level, i.e., (scenario, part) tuple level
    foreach::foreach(tuple = lstTuples) %dorng% {

        ## Assign values for readability
        scenarioCount <- tuple[1]
        partCount     <- tuple[2]
        cat("### Generating Scenario", scenarioCount,
            "Part", partCount, "\n")

        ## Generate the partCount-th subfile for the scenarioCount-th scenario.
        try(GenerateDistResNetRtimesAndSave(
            ScenarioDistResNet = Scenarios[[scenarioCount]],
            R                  = R,
            scenarioCount      = scenarioCount,
            partCount          = partCount))

        ## No return value
        NULL
    }

    ## This is a void function that returns NULL.
    NULL
}
