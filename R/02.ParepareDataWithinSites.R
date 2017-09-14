################################################################################
### Functions for data preparation within each data site
##
## Created on: 2016-04-28
## Author: Kazuki Yoshida
################################################################################


###
### Within-site data processing
################################################################################

###  Summary score estimation

##' Estimate summary confounder scores (binary logistic model)
##'
##' Estimate propensity score or disease risk score depending on the model specification. This is for binary outcomes.
##'
##' @param formula formula to specify the model used to estimate the summary score. For propensity score, specify the treatment model. For disease risk score, specify the outcome model in the absence of treatment.
##' @param data data frame of data for model building. Specify the entire dataset for propensity score estimation. Specify the untreated subset for disease risk score estimation.
##' @param newdata data frame of data for prediction. Specify the entire dataset.
##'
##' @return vector of summary scores predicted for \code{newdata}.
##'
##' @author Kazuki Yoshida
##'
##' @export
EstimateScoreBin <- function(formula, data, newdata) {
    ## Sanity check
    assert_that(class(formula) == "formula")
    assert_that(is.data.frame(data))
    assert_that(is.data.frame(newdata))

    ## Try to fit binary logistic model
    model1 <- try(glm(formula = formula,
                      family  = binomial(link = "logit"),
                      data    = data))

    ## Condition on model fit success
    if (is.error(model1)) {
        ## If model fitting fails, return NA's
        out <- rep(NA, nrow(newdata))

    } else {
        ## Otherwise, predict
        ## Probability of the predicted "outcome"
        out <- predict(object = model1, newdata = newdata, type = "response")
    }

    ## Make sure prediction for all observations
    assert_that(length(out) == nrow(newdata))

    as.numeric(out)
}


##' Estimate disease risk score (Cox model)
##'
##' Estimate disease risk score using Cox proportional hazards model. This is for survival outcomes. The disease risk score is defined as the linear predictor from the \code{coxph} model object. This is the log hazard ratio (range [-Inf, Inf]) of the event compared to the mean covariate individual.
##'
##' @param formula formula to specify the model used to estimate the summary score. For disease risk score, specify the outcome model in the absence of treatment.
##' @param data data frame of data for model building. Specify the untreated subset for disease risk score estimation.
##' @param newdata data frame of data for prediction. Specify the entire dataset.
##'
##' @return vector of linear predictors (log hazard ratios) for \code{newdata}.
##'
##' @author Kazuki Yoshida
##'
##' @export
EstimateScoreSurv <- function(formula, data, newdata) {
    ## Sanity check
    assert_that(class(formula) == "formula")
    assert_that(is.data.frame(data))
    assert_that(is.data.frame(newdata))

    ## Fit binary logistic model
    model1 <- try(coxph(formula = formula,
                        data    = data))

    ## Condition on model fit success
    if (is.error(model1)) {
        ## If model fitting fails, return NA's
        out <- rep(NA, nrow(newdata))

    } else {
        ## Otherwise, predict
        ## Estimated linear predictor (log hazard ratio) for each individual
        out <- predict(object = model1, newdata = newdata, type = "lp")
    }

    ## Make sure prediction for all observations
    assert_that(length(out) == nrow(newdata))

    as.numeric(out)
}


###  IPW creation

##' Create IPW
##'
##' Create inverse probability weights given probability of status assignment and observed status.
##'
##' @param score vector of predicted probability of being assigned \code{status == 1}.
##' @param status vector of observed status (binary 0, 1).
##' @param type character value containing either one of \code{naive} (non-stabilized naive IPW), \code{stabilized} (IPW stabilized by marginal probabilities of \code{status}), or \code{mw} (matching weights; stabilized by \code{pmin(score, 1 - score)}).
##'
##' @return vector of IPW.
##'
##' @author Kazuki Yoshida
##'
##' @export
CreateIpw <- function(score, status, type) {
    assert_that(length(score) == length(status))
    ## Only allow 0,1 status
    assert_that(all(status %in% c(0,1)))
    ## Check for valid type
    assert_that(class(type) == "character")
    assert_that(all(type %in% c("naive","stabilized","mw")))

    ## Naive IPW
    ipw <- ((1 / score) * status) + ((1 / (1 - score)) * (1 - status))

    ## Choose depending on specification
    if (type == "stabilized") {
        ## Marginal probabilities of status
        pStatus0 <- mean(status == 0)
        pStatus1 <- mean(status == 1)
        ## Stabilize
        ipw * ((status * pStatus1) + ((1 - status) * pStatus0))

    } else if (type == "mw") {
        ## Matching weights
        ipw * pmin(score, 1 - score)

    } else if (type == "naive") {
        ipw
    }
}


###  Strata creation

##' Create strata
##'
##' Create strata by categorizing summary scores into \code{nStrata} quantiles. If there are fewer discrete scores than \code{nStrata}, the rank will be returned. Sometimes a single value may have many data points causing quantile threshold to be identical values. The duplicated threshold values are dropped in this case, giving levels fewer than nStrata.
##'
##' @param score summary score to be categorized.
##' @param nStrata number of strata.
##'
##' @return vector strata as numeric
##'
##' @author Kazuki Yoshida
##'
##' @export
CreateStrata <- function(score, nStrata) {
    assert_that(is.numeric(score))
    assert_that(nStrata > 0)

    ## Create a factor version for discreteness assessment
    scoreFactor <- factor(score)

    if (all(is.na(score))) {
        ## If all scores are NA, just give NA
        scoreCat <- rep(as.integer(NA), length(score))

    } else if (nlevels(scoreFactor) <= nStrata) {
        ## If there are fewer valid discrete values than nStrata, just use ranks
        scoreCat <- as.integer(scoreFactor)

    } else {
        ## Otherwise peform cutting

        ## nStrata + 1 separators are required for nStrata intervals
        probs <- seq(from = 0, to = 1, length.out = (nStrata + 1))
        ## Quantile thresholds
        thresholds <- quantile(x = score, probs = probs)
        ## Do not allow duplication
        thresholds_nodup <- unique(thresholds)
        ## Message if duplication were dropped
        if (length(thresholds) != length(thresholds_nodup)) {
            message("CreateStrata: nStrata = ", nStrata, ", actual strata: ", length(thresholds_nodup) - 1)
        }

        ## Categorize score using non-duplicated thresholds
        scoreCat <- as.integer(cut(x = score,
                                   ## Separators must be number of categories + 1
                                   breaks = thresholds_nodup,
                                   ## This gives integers
                                   labels = FALSE,
                                   ## [a,b) FALSE means NOT closed on the right
                                   right = FALSE,
                                   ## The lowest value is included in the lowest category.
                                   include.lowest = TRUE))
    }

    scoreCat
}


###  Match creation

## MatchIt package supports nearest neighbor matching.
## Additionally caliper use needs some care.
## http://r.iq.harvard.edu/docs/matchit/2.4-14/Additional_Arguments_f3.html

## Matching package appears to do nearest neighbor matching, but it is not clearly stated.
##
## https://www.jstatsoft.org/article/view/v042i07

##' Match across status
##'
##'
##' @param score summary score vector used as matching metric.
##' @param status binary status variable vector to match across.
##' @param logit whether to perform matching on the logit scale. Use if using scores spanning [0,1] only.
##' @param caliper caliper multipliation factor. The actual caliper is \code{caliper * sd(score)}.
##'
##' @return vector of matched pair IDs containing NA for unmatched.
##'
##' @author Kazuki Yoshida
##'
##' @export
CreateMatches <- function(score, status, logit, caliper) {

    ## Sanity check
    assert_that(length(score) == length(status))
    assert_that(is.logical(logit))

    ## Logit transform score if requested
    if (logit) {
        score <- logit(score)
    }

    ## Condition on data characteristic to avoid errors in Match()
    if (all(is.na(score))) {
        ## If all scores are NA (score model failure), give NA's
        out <- rep(0, length(score))

    } else if (length(unique(status)) == 1) {
        ## If everybody is in one treatment group (non-existent treated or untreated), give NA's
        out <- rep(0, length(score))

    } else {
        ## Otherwise, proceed with matching
        m_out <- try(Matching::Match(Tr = status, X = score,
                                     ## 1:1 matching without replacement
                                     M = 1, replace = FALSE, ties = FALSE,
                                     ## Caliper definition
                                     caliper = caliper,
                                     ## This is just to avoid unnecessary warnings.
                                     estimand = ifelse(test = (mean(status) <= 0.5),
                                                       ## Fewer treated
                                                       yes = "ATT",
                                                       ## Fewer untreated
                                                       no = "ATC")))

        if (is.error(m_out)) {
            ## If matching procedure gives an error no matches
            out <- rep(0, length(score))

        } else {
            ## If a valid matching output is obtained, proceed with matching
            ## index of matched people
            index_matched <- c(m_out$index.treated, m_out$index.control)

            ## Create a 0 vector
            out <- as.numeric(rep(0, length(status)))

            ## Put 1 for observations matched
            out[index_matched] <- 1
        }
    }

    ## Return 0,1 vector indicating matched subjects
    out
}


###  Prepare all helper variables

##' Prepare PS- and DRS-related helper variables
##'
##' Add estimated PS and DRS to the dataset as well as a match indicator, weights, and strata based on these summary scores.
##'
##' @param x ResSite object containing a dataset from a single center given by \code{\link{GenerateOneCenter}}. All covariates (variable names starting with X) are used in PS and DRS estimation.
##'
##' @return ResSiteReady object, a data frame with additional PS, DRS, strata, and match indicator.
##'
##' @author Kazuki Yoshida
##'
##' @export
SitePrepareHelperVariables <- function(x) {
    ## Must be an unready ResSite object
    assert_that("ResSite" %in% class(x))
    assert_that(!("ResSiteReady" %in% class(x)))

    ## Extract covariate names (All variables starting with X are used)
    covNames <- Filter(f = function(elt) {grepl("^X", elt)}, x = names(x))
    ## Right-handed side of the formulas (covariates only)
    rhs <- paste(covNames, collapse = " + ")
    ## Generate formulas
    formulaA    <- paste0("A ~ ", rhs)
    formulaY    <- paste0("Y ~ ", rhs)
    formulaSurv <- paste0("Surv(time, event) ~ ", rhs)

    ## Summary score estimation ##
    ## Propensity score estimation
    x$ePs   <- EstimateScoreBin(formula  = as.formula(formulaA),
                                data     = x,
                                newdata  = x)
    ## Binary disease risk score estimation. range [0,1]
    x$eDrsB <- EstimateScoreBin(formula  = as.formula(formulaY),
                                ## Use untreated subset
                                data     = x[x$A == 0, ],
                                newdata  = x)
    ## Survival disease risk score estimation. range [-Inf,Inf]
    x$eDrsS <- EstimateScoreSurv(formula = as.formula(formulaSurv),
                                 ## Use untreated subset
                                 data    = x[x$A == 0, ],
                                 newdata = x)

    ## Weights creation ##
    ## IPTW
    x$ePsSIptw <- CreateIpw(score  = x$ePs,
                            status = x$A,
                            type   = "stabilized")
    ## MW
    x$ePsMw    <- CreateIpw(score  = x$ePs,
                            status = x$A,
                            type   = "mw")

    ## Stratification variable addition ##
    ## PS strata
    x$ePsStrata   <- CreateStrata(score = x$ePs, nStrata = 10)
    ## Binary DRS strata
    x$eDrsBStrata <- CreateStrata(score = x$eDrsB, nStrata = 10)
    ## Survival DRS strata
    x$eDrsSStrata <- CreateStrata(score = x$eDrsS, nStrata = 10)

    ## Matched cohort variable addition ##
    ## caliper is for 0.2 * SD(score)
    ## PS matched cohort indicator
    x$ePsMatch   <- CreateMatches(score = x$ePs,   status = x$A, logit = TRUE,  caliper = 0.2)
    ## Binary DRS matched cohort indicator (logit transform)
    x$eDrsBMatch <- CreateMatches(score = x$eDrsB, status = x$A, logit = TRUE,  caliper = 0.2)
    ## Survival DRS matched cohort indicator (use as is)
    x$eDrsSMatch <- CreateMatches(score = x$eDrsS, status = x$A, logit = FALSE, caliper = 0.2)

    ## Return data frame with additional class attribute for sanity checking
    class(x) <- c("ResSiteReady", class(x))
    x
}


###
### Site-specific regression
################################################################################

###  Binary outcome (Logistic regression)

## Helper function for SiteSpecificLogistic
## Model fitting only
FitSiteLogistic <- function(df) {
    ## Unadjusted
    unadj        <- try(glm(formula = Y ~ A,
                            family  = binomial(link = "logit"),
                            data    = ValidateDf(df = df, expo_var = "A", bin_var = "Y")))
    ## PS-matched cohort logistic regression
    ePsMatch     <- try(glm(formula = Y ~ A,
                            family  = binomial(link = "logit"),
                            data    = ValidateDf(df[df$ePsMatch == 1, ], expo_var = "A", bin_var = "Y")))
    ## DRS-matched cohort logistic regression
    eDrsBMatch   <- try(glm(formula = Y ~ A,
                            family  = binomial(link = "logit"),
                            data    = ValidateDf(df[df$eDrsBMatch == 1, ], expo_var = "A", bin_var = "Y")))

    ## PS-stratified conditional logistic regression
    ePsStrataE   <- try(clogit(formula = Y ~ A + strata(ePsStrata),
                               data    = ValidateDf(df = df, expo_var = "A", bin_var = "Y"),
                               method  = c("exact", "approximate", "efron", "breslow")[3]))
    ## DRS-stratified conditional logistic regression
    eDrsBStrataE <- try(clogit(formula = Y ~ A + strata(eDrsBStrata),
                               data    = ValidateDf(df = df, expo_var = "A", bin_var = "Y"),
                               method  = c("exact", "approximate", "efron", "breslow")[3]))

    ## PS-weighted cohort logistic regression
    ePsSIptw     <- try(svyglm(formula = Y ~ A,
                               family  = quasibinomial(link = "logit"),
                               design  = svydesign(ids = ~ 1,
                                                   weights = ~ ePsSIptw,
                                                   data = ValidateDf(df = df, expo_var = "A", bin_var = "Y"))))
    ## MW cohort logistic regression
    ePsMw        <- try(svyglm(formula = Y ~ A,
                               family  = quasibinomial(link = "logit"),
                               design  = svydesign(ids = ~ 1,
                                                   weights = ~ ePsMw,
                                                   data = ValidateDf(df = df, expo_var = "A", bin_var = "Y"))))

    ## Breslow
    ## PS-stratified conditional logistic regression
    ePsStrataB   <- try(clogit(formula = Y ~ A + strata(ePsStrata),
                               data    = ValidateDf(df, expo_var = "A", event_var = "Y"),
                               method  = c("exact", "approximate", "efron", "breslow")[4]))
    ## DRS-stratified conditional logistic regression
    eDrsBStrataB <- try(clogit(formula = Y ~ A + strata(eDrsBStrata),
                               data    = ValidateDf(df, expo_var = "A", event_var = "Y"),
                               method  = c("exact", "approximate", "efron", "breslow")[4]))

    ## Return model objects
    ## They all have appropriate coef and vcov methods.
    list(unadj        = unadj,
         ePsMatch     = ePsMatch,
         eDrsBMatch   = eDrsBMatch,
         ePsStrataE   = ePsStrataE,
         eDrsBStrataE = eDrsBStrataE,
         ePsSIptw     = ePsSIptw,
         ePsMw        = ePsMw,
         ePsStrataB   = ePsStrataB,
         eDrsBStrataB = eDrsBStrataB)
}

##' Run site-specific logistic regression
##'
##' Run site-specific logistic regression using both PS and DRS. Matched analyses, stratified analyses (conditional logistic regression via \code{survival::clogit}), and weighted analyses (sandwich variance estimator via the \code{survey} package) are conducted.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing point estimates and variance estimates.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRegressionBin <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Fit
    lstFits <- FitSiteLogistic(df)

    ## Extract coef/var
    lst <- lapply(lstFits,
                  CoefVar, variable = "A")

    ## Return as a data frame with method name
    out <- data.frame(method = names(lst),
                      coef   = sapply(lst, "[", 1),
                      var    = sapply(lst, "[", 2),
                      stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}


###  Survival outcome (Cox regression)
## Helper function for SiteSpecificCox
## Model fitting only
FitSiteCox <- function(df) {
    ## Unadjusted Cox regression
    unadjE       <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[3]))
    ## PS-matched cohort Cox regression
    ePsMatchE    <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df[df$ePsMatch == 1, ], expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[3]))
    ## DRS-matched cohort Cox regression
    eDrsSMatchE  <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df[df$eDrsSMatch == 1, ], expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[3]))

    ## PS-stratified Cox regression
    ePsStrataE   <- try(coxph(formula = Surv(time, event) ~ A + strata(ePsStrata),
                              data = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[3]))
    ## DRS-stratified Cox regression
    eDrsSStrataE <- try(coxph(formula = Surv(time, event) ~ A + strata(eDrsSStrata),
                              data = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[3]))

    ## PS-weighted cohort Cox regression
    ePsSIptwE    <- try(svycoxph(formula = Surv(time, event) ~ A,
                                 design  = svydesign(ids = ~ 1,
                                                     weights = ~ ePsSIptw,
                                                     data = ValidateDf(df, expo_var = "A", event_var = "event")),
                                 method = c("exact", "approximate", "efron", "breslow")[3]))
    ## MW cohort Cox regression
    ePsMwE       <- try(svycoxph(formula = Surv(time, event) ~ A,
                                 design  = svydesign(ids = ~ 1,
                                                     weights = ~ ePsMw,
                                                     data = ValidateDf(df, expo_var = "A", event_var = "event")),
                                 method = c("exact", "approximate", "efron", "breslow")[3]))

    ## Breslow
    ## Unadjusted Cox regression
    unadjB       <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[4]))
    ## PS-matched cohort Cox regression
    ePsMatchB    <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df[df$ePsMatch == 1, ], expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[4]))
    ## DRS-matched cohort Cox regression
    eDrsSMatchB  <- try(coxph(formula = Surv(time, event) ~ A,
                              data    = ValidateDf(df[df$eDrsSMatch == 1, ], expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[4]))

    ## PS-stratified Cox regression
    ePsStrataB   <- try(coxph(formula = Surv(time, event) ~ A + strata(ePsStrata),
                              data = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[4]))
    ## DRS-stratified Cox regression
    eDrsSStrataB <- try(coxph(formula = Surv(time, event) ~ A + strata(eDrsSStrata),
                              data = ValidateDf(df, expo_var = "A", event_var = "event"),
                              method = c("exact", "approximate", "efron", "breslow")[4]))

    ## PS-weighted cohort Cox regression
    ePsSIptwB    <- try(svycoxph(formula = Surv(time, event) ~ A,
                                 design  = svydesign(ids = ~ 1,
                                                     weights = ~ ePsSIptw,
                                                     data = ValidateDf(df, expo_var = "A", event_var = "event")),
                                 method = c("exact", "approximate", "efron", "breslow")[4]))
    ## MW cohort Cox regression
    ePsMwB       <- try(svycoxph(formula = Surv(time, event) ~ A,
                                 design  = svydesign(ids = ~ 1,
                                                     weights = ~ ePsMw,
                                                     data = ValidateDf(df, expo_var = "A", event_var = "event")),
                                 method = c("exact", "approximate", "efron", "breslow")[4]))

    ## Return model objects
    ## Do not save these. They include the original dataset and are large.
    list(unadjE        = unadjE,
         ePsMatchE     = ePsMatchE,
         eDrsSMatchE   = eDrsSMatchE,
         ePsStrataE    = ePsStrataE,
         eDrsSStrataE  = eDrsSStrataE,
         ePsSIptwE     = ePsSIptwE,
         ePsMwE        = ePsMwE,
         ## Breslow
         unadjB        = unadjB,
         ePsMatchB     = ePsMatchB,
         eDrsSMatchB   = eDrsSMatchB,
         ePsStrataB    = ePsStrataB,
         eDrsSStrataB  = eDrsSStrataB,
         ePsSIptwB     = ePsSIptwB,
         ePsMwB        = ePsMwB)

}

##' Run site-specific Cox regression
##'
##' Run site-specific Cox regression using both PS and DRS. Matched analyses, stratified analyses, and weighted analyses are conducted.
##'
##' @inheritParams SiteRegression
##'
##' @return list of pairs of a point estimate and a variance estimate for the treatment effect.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRegressionSurv <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Fit
    lstFits <- FitSiteCox(df)

    ## Extract coef/var
    lst <- lapply(lstFits,
                  CoefVar, variable = "A")

    ## Return as a data frame with method name
    out <- data.frame(method = names(lst),
                      coef   = sapply(lst, "[", 1),
                      var    = sapply(lst, "[", 2),
                      stringsAsFactors = FALSE)
    rownames(out) <- NULL
    out
}


###  Both outcomes

## Undocumented helper function
## Given a list of binary outcome df and survival outcome df,
## stack them up in a df with an additional column outcome.
NameOutcomesAndStackUp <- function(lstBinSurv) {
    binary   <- lstBinSurv[["binary"]]
    survival <- lstBinSurv[["survival"]]
    ## cbind outcome column and stack up
    dplyr::bind_rows(cbind(outcome = rep("binary", nrow(binary)),
                           binary,
                           stringsAsFactors = FALSE),
                     cbind(outcome = rep("survival", nrow(survival)),
                           survival,
                           stringsAsFactors = FALSE))
}

##' Run site-specific regression analyses
##'
##' Run site-specific regression analyses and return a data frame.
##'
##' @param df data frame for a single data center preapred by \code{\link{SitePrepareHelperVariables}}. It must have class \code{ResSiteReady}.
##'
##' @return list containing one list for the binary analyses (see \code{\link{SiteRegressionBin}}) and another for the survival analyses (see \code{\link{SiteRegressionSurv}}).
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRegression <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Call respective functions
    list(binary   = SiteRegressionBin(df),
         survival = SiteRegressionSurv(df)) %>%
        NameOutcomesAndStackUp
}


###
### Site-specific summary aggregation
################################################################################

###  Binary outcome

##' Summarize event counts and trial counts by groups.
##'
##' Summarize event counts and trial counts of binary outcome (coded 0 or 1) classified
##' by a grouping variable. Intended for summarizing these counts by the binary
##' treatment groups. When there are empty cells, the event counts and the total
##' counts are both forced to be zeros.
##'
##' @param var variable name of a binary variable (coded 0 or 1) in the data frame.
##' @param by variable name of a categorical variable used as a grouping variable. This should usually be a factor if called in a stratified analysis to avoid differential lengths.
##' @param data data frame containing these variables.
##'
##' @return data frame containing the \code{events} (event count) and the \code{denom} (trial count) information for each treatment group. When the data frame supplied is empty (0 rows), a single-row data frame filled with NA is returned.
##'
##' @author Kazuki Yoshida
##'
##' @export
SummarizeEventsTotalsBy <- function(var, by, data) {
    assert_that(length(var) == 1)
    assert_that(length(by) == 1)
    assert_that("data.frame" %in% class(data))
    ## Binary
    assert_that(all(data[,var] %in% c(0,1)))

    if (nrow(data) == 0) {
        ## If no row, return a tibble with 1 row filled with NA's
        out <- CreateNaDf(labels = c("A", "events", "denom"),
                          as_funcs = c(as.integer, as.integer, as.integer),
                          nrow = 1) %>% (dplyr::as_data_frame)

    } else {
        ## Otherwise construct summaries
        out <- data %>%
            ## Using by is ok. It's an argument.
            dplyr::group_by_(.dots = c(by)) %>%
            dplyr::summarize_(events = as.formula(paste0("~ sum(", var, ")")),
                              denom = as.formula(paste0("~ length(", var, ")")))
    }
}

##' Summarize stratum-specific event counts and trial counts by groups.
##'
##' Summarize stratum-specific event counts and trial counts by groups and return
##' a list of lists containing events and totals.
##'
##' @inheritParams SummarizeEventsTotalsBy
##' @param strata variable name of a stratifying variable.
##'
##' @return data frame containing one stratum and treatment-specific row of \code{events} (event count) and \code{denom} (trial count) for each stratum. When the data frame supplied is empty (0 rows), a single-row data frame filled with NA is returned.
##'
##' @author Kazuki Yoshida
##'
##' @export
SummarizeStratumEventsTotalsBy <- function(var, by, strata, data) {
    assert_that(length(var) == 1)
    assert_that(length(by) == 1)
    assert_that(length(strata) == 1)
    assert_that("data.frame" %in% class(data))
    ## Binary
    assert_that(all(data[,var] %in% c(0,1)))

    if (nrow(data) == 0) {
        ## If no row, return a tibble with 1 row filled with NA's
        out <- CreateNaDf(labels = c("strata", "A", "events", "denom"),
                          as_funcs = c(as.integer, as.integer, as.integer, as.integer),
                          nrow = 1) %>% (dplyr::as_data_frame)

    } else {

        out <- data %>%
            ## Using strata, by is ok. They are arguments.
            dplyr::group_by_(.dots = c(strata, by)) %>%
            dplyr::summarize_(events = as.formula(paste0("~ sum(", var, ")")),
                              denom = as.formula(paste0("~ length(", var, ")")))
        names(out)[1] <- "strata"
    }

    out
}


##' Obtain site-specific summaries for binary outcome
##'
##' Obtain site-specific summaries for binary outcome using PS and DRS-matched cohorts as well as PS and DRS-stratified cohorts.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing results for PS/DRS-matched cohorts and PS/DRS-stratified cohorts.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteSummaryBin <- function(df) {

    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Unadjusted summaries
    unadj       <- SummarizeEventsTotalsBy(var = "Y", by = "A", data = df)

    ## Matched summaries are calculated in appropriate subsets
    ePsMatch    <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                           data = df[df$ePsMatch %in% 1, ])
    eDrsBMatch  <- SummarizeEventsTotalsBy(var = "Y", by = "A",
                                           data = df[df$eDrsBMatch %in% 1, ])

    ## Stratified summaries
    ePsStrata   <- SummarizeStratumEventsTotalsBy(var = "Y", by = "A", strata = "ePsStrata",
                                                  data = df[!is.na(df$ePsStrata),])
    eDrsBStrata <- SummarizeStratumEventsTotalsBy(var = "Y", by = "A", strata = "eDrsBStrata",
                                                  data = df[!is.na(df$eDrsBStrata),])

    ## Combination
    ## Add NA strata for unstratified results
    unadj$strata       <- NA
    ePsMatch$strata    <- NA
    eDrsBMatch$strata  <- NA
    ## Add method names
    unadj$method       <- "unadj"
    ePsMatch$method    <- "ePsMatch"
    eDrsBMatch$method  <- "eDrsBMatch"
    ePsStrata$method   <- "ePsStrata"
    eDrsBStrata$method <- "eDrsBStrata"
    ## row combine
    dplyr::bind_rows(unadj,
                     ePsMatch,
                     eDrsBMatch,
                     ePsStrata,
                     eDrsBStrata)[c("method","strata","A","events","denom")] %>%
        as.data.frame
}


###  Survival outcome

##' Summarize event counts and total person-time by groups.
##'
##' Summarize event counts and total person-time of survival outcome (coded 0 or 1) classified
##' by a grouping variable. Intended for summarizing these counts by the binary
##' treatment groups. When there are empty cells, the event counts and the total
##' person-time are both forced to be zeros.
##'
##' @param event variable name of a binary event indicator variable (coded 0 or 1) in the data frame.
##' @param time variable name of a follow-up time variable in the data frame.
##' @inheritParams SummarizeEventsTotalsBy
##'
##' @return data frame containing the \code{events} (event count) element the \code{denom} (person-time) information for each treatment group. When the data frame supplied is empty (0 rows), a single-row data frame filled with NA is returned.
##'
##' @author Kazuki Yoshida
##'
##' @export
SummarizeEventsPtBy <- function(event, time, by, data) {
    assert_that(length(event) == 1)
    assert_that(length(time) == 1)
    assert_that(length(by) == 1)
    assert_that("data.frame" %in% class(data))
    ## Binary
    assert_that(all(data[,event] %in% c(0,1)))
    ## Valid time
    assert_that(all(data[,time] >= 0))

    if (nrow(data) == 0) {
        ## If no row, return a tibble with 1 row filled with NA's
        out <- CreateNaDf(labels = c("A", "events", "denom"),
                          as_funcs = c(as.integer, as.integer, as.integer),
                          nrow = 1) %>% (dplyr::as_data_frame)

    } else {

        out <- data %>%
            ## Using by is ok. It's an argument.
            dplyr::group_by_(.dots = c(by)) %>%
            dplyr::summarize_(events = as.formula(paste0("~ sum(", event, ")")),
                              denom  = as.formula(paste0("~ sum(", time, ")")))
    }

    out
}

##' Summarize stratum-specific event counts and total person-time by groups.
##'
##' Summarize stratum-specific event counts and total person-time by groups and return
##' a list of lists containing events and ptime.
##'
##' @inheritParams SummarizeEventsPtBy
##' @param strata variable name of a stratifying variable.
##'
##' @return data frame containing one stratum and treatment-specific row of \code{events} (event count) and \code{denom} (person-time) for each stratum. When the data frame supplied is empty (0 rows), a single-row data frame filled with NA is returned.
##'
##' @author Kazuki Yoshida
##'
##' @export
SummarizeStratumEventsPtBy <- function(event, time, by, strata, data) {
    assert_that(length(event) == 1)
    assert_that(length(time) == 1)
    assert_that(length(by) == 1)
    assert_that("data.frame" %in% class(data))
    ## Binary
    assert_that(all(data[,event] %in% c(0,1)))
    ## Valid time
    assert_that(all(data[,time] >= 0))

    if (nrow(data) == 0) {
        ## If no row, return a tibble with 1 row filled with NA's
        out <- CreateNaDf(labels = c("strata", "A", "events", "denom"),
                          as_funcs = c(as.integer, as.integer, as.integer, as.integer),
                          nrow = 1) %>% (dplyr::as_data_frame)

    } else {

        out <- data %>%
            ## Using strata, by is ok. They are arguments.
            dplyr::group_by_(.dots = c(strata, by)) %>%
            dplyr::summarize_(events = as.formula(paste0("~ sum(", event, ")")),
                              denom  = as.formula(paste0("~ sum(", time, ")")))
        names(out)[1] <- "strata"
    }

    out
}


##' Obtain site-specific summaries for survival outcome
##'
##' Obtain site-specific summaries for survival outcome using PS and DRS-matched cohorts as well as PS and DRS-stratified cohorts.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing results for PS/DRS-matched cohorts and PS/DRS-stratified cohorts.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteSummarySurv <- function(df) {

    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Unadjusted summaries
    unadj       <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                       data = df)

    ## Matched summaries are calculated in appropriate subsets
    ePsMatch    <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                       data = df[df$ePsMatch %in% 1, ])
    eDrsSMatch  <- SummarizeEventsPtBy(event = "event", time = "time", by = "A",
                                       data = df[df$eDrsSMatch %in% 1, ])

    ## Stratified summaries
    ePsStrata   <- SummarizeStratumEventsPtBy(event = "event", time = "time", by = "A",
                                              strata = "ePsStrata",
                                              data = df[!is.na(df$ePsStrata),])
    eDrsSStrata <- SummarizeStratumEventsPtBy(event = "event", time = "time", by = "A",
                                              strata = "eDrsSStrata",
                                              data = df[!is.na(df$eDrsSStrata),])

    ## Combination
    ## Add NA strata for unstratified results
    unadj$strata       <- NA
    ePsMatch$strata    <- NA
    eDrsSMatch$strata  <- NA
    ## Add method names
    unadj$method       <- "unadj"
    ePsMatch$method    <- "ePsMatch"
    eDrsSMatch$method  <- "eDrsSMatch"
    ePsStrata$method   <- "ePsStrata"
    eDrsSStrata$method <- "eDrsSStrata"
    ## row combine
    dplyr::bind_rows(unadj,
                     ePsMatch,
                     eDrsSMatch,
                     ePsStrata,
                     eDrsSStrata)[c("method","strata","A","events","denom")] %>%
        as.data.frame
}


###  Both outcomes

##' Run site-specific summary data creation
##'
##' Run site-specific summary data creation and return a data frame.
##'
##'
##' @inheritParams SiteRegression
##'
##' @return list containing one list for the binary analyses (see \code{\link{SiteRegressionBin}}) and another for the survival analyses (see \code{\link{SiteRegressionSurv}}).
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteSummary <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Call respective functions
    list(binary   = SiteSummaryBin(df),
         survival = SiteSummarySurv(df)) %>%
        NameOutcomesAndStackUp
}


###
### Site-specific riskset construction
################################################################################

###  Helper functions that should be called within strata (or matched data)

SiteRisksetsHelperUnweighted <- function(time = NULL, event, A, eval_times = NULL) {
    ## event and treatment must be of the same length
    assert_that(length(event) == length(A))
    ## event indicator has to be binary
    assert_that(all(event %in% c(0,1)))
    ## exposure indicator has to be binary
    assert_that(all(A %in% c(0,1)))

    ## time must be positive if existing
    if (is.null(time)) {
        ## If it does not exist, assume a binary outcome.
        ## impute time 0 for all as a dummy time variable.
        time <- rep(0, length(event))
        ## Evaluate at 0 only
        eval_times <- 0

    } else {
        ## If it exists, check the length
        assert_that(length(time) == length(event))
        ## All positive
        assert_that(all(time >= 0))
        assert_that(all(eval_times >= 0))
    }

    ## Construct time points to evaluate risk sets and events
    eval_times <- sort(unique(c(0, eval_times)))

    ## FIXME: This is likely very inefficient.
    ## loop over event times
    out <- lapply(seq_along(eval_times), function(i) {

        ## Current evaluation time point
        eval_time <- eval_times[i]
        ## Previous evaluation time point
        if (i == 1) {
            ## If evaluation at first point (should be time 0),
            ## give some smaller number
            prev_time <- eval_time - 1
        } else {
            ## Otherwise give the previous time point
            prev_time <- eval_times[i - 1]
        }

        ## One-row df
        data.frame(eval_time = eval_time,
                   ## Cumulative events
                   ##  Number of events happening at or before eval_time
                   events_A0 = sum(event[A == 0 & time > prev_time & time <= eval_time]),
                   events_A1 = sum(event[A == 1 & time > prev_time & time <= eval_time]),
                   ## Current risk sets
                   ##  Number of individuals remaining at eval_time
                   riskset_A0    = sum(A == 0 & time >= eval_time),
                   riskset_A1    = sum(A == 1 & time >= eval_time))

    }) %>% do.call(rbind, .)

    ## Return df without intermediate variable
    out
}


## Requires weights
SiteRisksetsHelperWeighted <- function(time = NULL, event, A, W, eval_times = NULL) {
    ## event and treatment must be of the same length
    assert_that(length(event) == length(A))
    assert_that(length(event) == length(W))
    ## event indicator has to be binary
    assert_that(all(event %in% c(0,1)))
    ## exposure indicator has to be binary
    assert_that(all(A %in% c(0,1)))
    ## weights must be nonnegative
    assert_that(all(W >= 0))

    ## time must be positive if existing
    if (is.null(time)) {
        ## If it does not exist, assume a binary outcome.
        ## impute time 0 for all as a dummy time variable.
        time <- rep(0, length(event))
        ## Evaluate at 0 only
        eval_times <- 0

    } else {
        ## If it exists, check the length
        assert_that(length(time) == length(event))
        ## All positive
        assert_that(all(time >= 0))
        assert_that(all(eval_times >= 0))
    }

    ## Construct time points to evaluate risk sets and events
    ## Always include 0
    eval_times <- sort(unique(c(0, eval_times)))

    ## FIXME: This is likely very inefficient.
    ## loop over event times
    out <- lapply(seq_along(eval_times), function(i) {

        ## Current evaluation time point
        eval_time <- eval_times[i]
        ## Previous evaluation time point
        if (i == 1) {
            ## If evaluation at first point (should be time 0),
            ## give some smaller number
            prev_time <- eval_time - 1
        } else {
            ## Otherwise give the previous time point
            prev_time <- eval_times[i - 1]
        }

        ## One-row df
        data.frame(
            eval_time = eval_time,
            ## Cumulative events
            ##  Number of weighted events happening at or before eval_time after prev_time
            w_events_A0 = sum(W[event == 1 & A == 0 &
                                time > prev_time & time <= eval_time]),
            w_events_A1 = sum(W[event == 1 & A == 1 &
                                time > prev_time & time <= eval_time]),
            ## Current risk sets
            ##  Number of weighted individuals remaining at eval_time
            w_riskset_A0    = sum(W[A == 0 & time >= eval_time]),
            w_riskset_A1    = sum(W[A == 1 & time >= eval_time]),
            ## Variances
            ## Variance of weights for events happening at eval_time
            varw_events_A0  = safe_var(W[event == 1 & A == 0 &
                                         time > prev_time & time <= eval_time]),
            varw_events_A1  = safe_var(W[event == 1 & A == 1 &
                                         time > prev_time & time <= eval_time]),
            ## Variance of weights for nonevent individuals remaining at eval_time
            ## Those who are in the risk set at eval_time, but NOT having event at eval_time
            ## Not safe if not all event times are in eval_times.
            ## Do not eval at coarser intervals than events.
            varw_nonevents_A0 = safe_var(W[A == 0 & time >= eval_time &
                                           ## Not having event exactly at event_time
                                           !(time == eval_time & event == 1)]),
            varw_nonevents_A1 = safe_var(W[A == 1 & time >= eval_time &
                                           ## Not having event exactly at event_time
                                           !(time == eval_time & event == 1)]))

    }) %>% do.call(rbind, .)

    ## Return df without intermediate variable
    out
}

##' Construct risk set data
##'
##' Construct risk set data classified by exposure status given
##' the time variable, event status variable, exposure status
##' variable, as well as time points to evaluate events and risk set
##' sizes at. The time points for evaluation must include the observed
##' event times. Otherwise, an interval may contain events and the
##' riskset size definition at the end of the interval become ambiguous.
##' As the observed time variable in the simulated dataset is coded in
##' integer days, specifying all days in the intended follow up duration
##' will assure inclusion of all event times.
##'
##' @param time vector of the observed time variable. If omitted, taken to be all 0. Omit when \code{event} is a binary outcome variable without an accompanying time variable.
##' @param event vector of the event status binary variable (must be 0, 1). A binary outcome variable can also be used if the outcome of interest is such a variable.
##' @param A vector of the exposure status variable (must be 0, 1).
##' @param W vector of the weights. Omit if not weighting.
##' @param eval_times vector of time points to evaluate risk sets and event counts. Time 0 is always included. Required if \code{time} is given. When \code{time} is not given, it is taken to be 0 only.
##'
##' @return reduced risk set data frame having one row for each unique evaluation time.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRisksetsHelper <- function(time = NULL, event, A, W = NULL, eval_times = NULL) {

    ## Always conduct unweighted riskset construction
    out <- SiteRisksetsHelperUnweighted(time = time, event = event, A = A,
                                        eval_times = eval_times)
    nrow_out <- nrow(out)

    ## If weights are given, additionally conduct weighted risksets construction
    if (!is.null(W)) {
        out_w <- SiteRisksetsHelperWeighted(time = time, event = event, A = A,
                                            W = W,
                                            eval_times = eval_times)
        ## Merge these two
        out <- merge(x = out,
                     y = out_w,
                     by = "eval_time",
                     ## outer joint
                     all.x = TRUE,  all.y = TRUE)

        ## Check the row count
        assert_that(nrow_out == nrow(out))
    }

    ## Return
    out
}

###  Helper functions that are strata-aware

## Actual worker for non-empty data
SiteRisksetsByStrataHelperValid <- function(time = NULL, event, A, W = NULL, strata = NULL) {
    ## strata and A must be of the same length if strata exists
    if (!is.null(strata)) {
        assert_that(length(strata) == length(A))
        assert_that(all(!is.na(strata)))
    }

    ## Call helper
    if (is.null(strata)) {
        ## If not stratifying, just call the helper
        out <- SiteRisksetsHelper(time = time, event = event, A = A,
                                  W = W,
                                  ## Evaluate at relevant time points when events occurred
                                  eval_times = time[event == 1])
        out$strata <- as.integer(NA)

    } else {
        ## If stratifying call for each stratum
        out <- lapply(sort(unique(strata)), function(stratum) {
            out <- SiteRisksetsHelper(time  = time[strata == stratum],
                                      event = event[strata == stratum],
                                      A     = A[strata == stratum],
                                      W     = W[strata == stratum],
                                      eval_times = (time[strata == stratum])[(event[strata == stratum]) == 1])
            out$strata <- stratum
            out
        }) %>% do.call(what = rbind)
    }

    ## Record where the strata column is
    posStrata <- which(names(out) == "strata")

    ## Bring strata column to the left
    cbind(out[c(posStrata)],
          ## Add the reminder
          out[-posStrata])
}

## Actual worker for empty data to return NA data frame
SiteRisksetsByStrataHelperEmpty <- function(time = NULL, event, A, W = NULL, strata = NULL) {
    ## time variable is irrelevant, as the returned df always has a single row.

    if (is.null(W) & is.null(strata)) {
        ## Empty matched data
        out <- CreateNaDf(labels = c("strata","eval_time",
                                     "events_A0","events_A1","riskset_A0","riskset_A1"),
                          as_funcs = c(as.integer, as.double,
                                       as.integer, as.integer, as.integer, as.integer),
                          nrow = 1)

    } else if (!is.null(strata)) {
        ## Empty stratified data (same as above)
        out <- CreateNaDf(labels = c("strata","eval_time",
                                     "events_A0","events_A1","riskset_A0","riskset_A1"),
                          as_funcs = c(as.integer, as.double,
                                       as.integer, as.integer, as.integer, as.integer),
                          nrow = 1)

    } else if (!is.null(W)) {
        ## Empty weighted data
        CreateNaDf(labels = c("strata","eval_time",
                              "events_A0","events_A1","riskset_A0","riskset_A1",
                              "w_events_A0", "w_events_A1", "w_riskset_A0", "w_riskset_A1",
                              "varw_events_A0", "varw_events_A1", "varw_nonevents_A0", "varw_nonevents_A1"),
                   as_funcs = c(as.integer, as.double,
                                as.integer, as.integer, as.integer, as.integer,
                                as.double, as.double, as.double, as.double,
                                as.double, as.double, as.double, as.double),
                   nrow = 1)

    } else {
        error("Both strata and W are supplied. Anomalous data!")
    }
}

##' Construct risk set data stratifying on a variable
##'
##' Gives a stratified risk set data. See \code{\link{SiteRisksetsHelper}} for details
##' about the risk set data within each stratum. If not \code{strata} vector is supplied,
##' one stratum with value NA is assumed.
##'
##' @inheritParams SiteRisksetsHelper
##' @param strata vector of stratifying variable
##'
##' @return reduced stratified risk set data frame having one row for each unique evaluation time for each strata.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRisksetsByStrata <- function(time = NULL, event, A, W = NULL, strata = NULL) {
    ## strata and A must be of the same length if strata exists
    if (!is.null(strata)) {
        assert_that(length(strata) == length(A))
        assert_that(all(!is.na(strata)))
    }

    if (length(event) == 0 | length(A) == 0) {
        ## If empty data are given, call NA riskset data creator
        out <- SiteRisksetsByStrataHelperEmpty(time = time, event = event, A = A, W = W, strata = strata)

    } else {
        ## Otherwise called the true worker function
        out <- SiteRisksetsByStrataHelperValid(time = time, event = event, A = A, W = W, strata = strata)
    }

    out
}


###  Binary outcome

##' Obtain site-specific risk sets for binary outcome
##'
##' Obtain site-specific risk sets for binary outcome using PS and DRS-matched cohorts as well as PS and DRS-stratified cohorts.
##'
##' @inheritParams SiteRegression
##'
##' @return list containing results for PS/DRS-matched risk sets and PS/DRS-stratified risk sets.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRisksetsBin <- function(df) {

    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## IMPORTANT:
    ## event variable is outcome Y
    ## A variable is treatment A
    ## Do not specify time variables

    ## Unadjusted
    unadj       <- SiteRisksetsByStrata(event  = df[, "Y"],
                                        A      = df[, "A"])

    ## Matched summaries are calculated in appropriate subsets
    ePsMatch    <- SiteRisksetsByStrata(event  = df[df$ePsMatch %in% 1, "Y"],
                                        A      = df[df$ePsMatch %in% 1, "A"])
    eDrsBMatch  <- SiteRisksetsByStrata(event  = df[df$eDrsBMatch %in% 1, "Y"],
                                        A      = df[df$eDrsBMatch %in% 1, "A"])

    ## Stratified summaries
    ePsStrata   <- SiteRisksetsByStrata(event  = df[!is.na(df$ePsStrata), "Y"],
                                        A      = df[!is.na(df$ePsStrata), "A"],
                                        strata = df[!is.na(df$ePsStrata), "ePsStrata"])
    eDrsBStrata <- SiteRisksetsByStrata(event  = df[!is.na(df$eDrsBStrata), "Y"],
                                        A      = df[!is.na(df$eDrsBStrata), "A"],
                                        strata = df[!is.na(df$eDrsBStrata), "eDrsBStrata"])

    ## Weighted summaries
    ePsSIptw    <- SiteRisksetsByStrata(event  = df[!is.na(df$ePsSIptw), "Y"],
                                        A      = df[!is.na(df$ePsSIptw), "A"],
                                        W      = df[!is.na(df$ePsSIptw), "ePsSIptw"])
    ePsMw       <- SiteRisksetsByStrata(event  = df[!is.na(df$ePsMw), "Y"],
                                        A      = df[!is.na(df$ePsMw), "A"],
                                        W      = df[!is.na(df$ePsMw), "ePsMw"])

    ## Name methods
    unadj$method       <- "unadj"
    ePsMatch$method    <- "ePsMatch"
    eDrsBMatch$method  <- "eDrsBMatch"
    ePsStrata$method   <- "ePsStrata"
    eDrsBStrata$method <- "eDrsBStrata"
    ePsSIptw$method    <- "ePsSIptw"
    ePsMw$method       <- "ePsMw"

    ## row combine
    out <- dplyr::bind_rows(unadj,
                            ePsMatch,
                            eDrsBMatch,
                            ePsStrata,
                            eDrsBStrata,
                            ePsSIptw,
                            ePsMw)

    ## Move method to first column
    posMethod <- which(names(out) == "method")
    cbind(out[posMethod],
          out[-posMethod])
}


###  Survival outcome

##' Obtain site-specific risk sets for survival outcome
##'
##' Obtain site-specific risk sets for survival outcome using PS and DRS-matched cohorts as well as PS and DRS-stratified cohorts.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing point estimates and variance estimates.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRisksetsSurv <- function(df) {

    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Unadjusted
    unadj       <- SiteRisksetsByStrata(time   = df[, "time"],
                                        event  = df[, "event"],
                                        A      = df[, "A"])

    ## Matched summaries are calculated in appropriate subsets
    ePsMatch    <- SiteRisksetsByStrata(time   = df[df$ePsMatch %in% 1, "time"],
                                        event  = df[df$ePsMatch %in% 1, "event"],
                                        A      = df[df$ePsMatch %in% 1, "A"])
    eDrsSMatch  <- SiteRisksetsByStrata(time   = df[df$eDrsSMatch %in% 1, "time"],
                                        event  = df[df$eDrsSMatch %in% 1, "event"],
                                        A      = df[df$eDrsSMatch %in% 1, "A"])

    ## Stratified summaries
    ePsStrata   <- SiteRisksetsByStrata(time   = df[!is.na(df$ePsStrata), "time"],
                                        event  = df[!is.na(df$ePsStrata), "event"],
                                        A      = df[!is.na(df$ePsStrata), "A"],
                                        strata = df[!is.na(df$ePsStrata), "ePsStrata"])
    eDrsSStrata <- SiteRisksetsByStrata(time   = df[!is.na(df$eDrsSStrata), "time"],
                                        event  = df[!is.na(df$eDrsSStrata), "event"],
                                        A      = df[!is.na(df$eDrsSStrata), "A"],
                                        strata = df[!is.na(df$eDrsSStrata), "eDrsSStrata"])

    ## Weighted summaries
    ePsSIptw    <- SiteRisksetsByStrata(time   = df[!is.na(df$ePsSIptw), "time"],
                                        event  = df[!is.na(df$ePsSIptw), "event"],
                                        A      = df[!is.na(df$ePsSIptw), "A"],
                                        W      = df[!is.na(df$ePsSIptw), "ePsSIptw"])
    ePsMw       <- SiteRisksetsByStrata(time   = df[!is.na(df$ePsMw), "time"],
                                        event  = df[!is.na(df$ePsMw), "event"],
                                        A      = df[!is.na(df$ePsMw), "A"],
                                        W      = df[!is.na(df$ePsMw), "ePsMw"])

    ## Name methods
    unadj$method       <- "unadj"
    ePsMatch$method    <- "ePsMatch"
    eDrsSMatch$method  <- "eDrsSMatch"
    ePsStrata$method   <- "ePsStrata"
    eDrsSStrata$method <- "eDrsSStrata"
    ePsSIptw$method    <- "ePsSIptw"
    ePsMw$method       <- "ePsMw"

    ## row combine
    out <- dplyr::bind_rows(unadj,
                            ePsMatch,
                            eDrsSMatch,
                            ePsStrata,
                            eDrsSStrata,
                            ePsSIptw,
                            ePsMw)

    ## Move method to first column
    posMethod <- which(names(out) == "method")
    cbind(out[posMethod],
          out[-posMethod])
}


###  Both outcomes

##' Run site-specific risk set data creation
##'
##' Run site-specific risk set data creation and return a data frame.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing results for the binary and survival outcomes.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteRisksets <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Call respective functions
    list(binary   = SiteRisksetsBin(df),
         survival = SiteRisksetsSurv(df)) %>%
        NameOutcomesAndStackUp
}


###
### Site-specific individual-level data extraction
################################################################################

###  Both outcomes

##' Run site-specific individual-level data extraction
##'
##' Run site-specific individual-level data extraction and return a data frame.
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing individual-level observed data.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteDataset <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Extract data for sharing
    df[c("A", "Y", "time", "event",
         "ePs", "eDrsB", "eDrsS",
         "ePsSIptw", "ePsMw",
         "ePsStrata", "eDrsBStrata", "eDrsSStrata",
         "ePsMatch", "eDrsBMatch", "eDrsSMatch")]
}


##' Run site-specific individual-level data extraction (counterfactuals)
##'
##' Run site-specific individual-level data extraction including all variables both observed and latent (e.g., counterfactuals).
##'
##' @inheritParams SiteRegression
##'
##' @return data frame containing individual-level data including counterfactuals.
##'
##' @author Kazuki Yoshida
##'
##' @export
SiteTruth <- function(df) {
    ## Sanity check
    assert_that(("ResSiteReady" %in% class(df)))

    ## Extract data for sharing
    df
}
