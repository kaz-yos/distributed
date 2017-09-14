################################################################################
### Functions for assessing methods
##
## Created on: 2016-04-23
## Author: Kazuki Yoshida
################################################################################

### Extract matching names
ExtColNames <- function(pat, data) {
    Filter(f = function(elt) {grepl(pat, elt)}, x = colnames(data))
}


### Calcualte average coefficient
AverageCoefs <- function(coefData) {
    out <- colMeans(coefData)
    names(out) <- gsub("coef", "meanCoef", x = names(out))
    out
}


### Calcualte average estimated variance
AverageVars <- function(varData) {
    out <- colMeans(varData)
    names(out) <- gsub("var", "meanVar", x = names(out))
    out
}


### Calculate true variance
TrueVars <- function(coefData) {
    out <- apply(coefData, MARGIN = 2, FUN = var)
    names(out) <- gsub("coef", "trueVar", x = names(out))
    out
}


### Calculate bias
Bias <- function(trueCoefs, meanCoefs) {
    repCount <- round(length(meanCoefs) / length(trueCoefs))
    trueCoefs <- rep(trueCoefs, repCount)
    out <- meanCoefs - trueCoefs
    names(out) <- gsub("meanCoef", "bias", names(meanCoefs))
    out
}


### Calculate MSE
Mse <- function(trueVars, biases) {
    names(biases) <- NULL
    out <- trueVars + biases^2
    names(out) <- gsub("trueVar", "mse", names(trueVars))
    out
}


### Calculate rejection rate
RejRate <- function(alpha, pvalData) {
    out <- colMeans(pvalData > alpha)
    names(out) <- gsub("pval", "rejRate", names(out))
    out
}


### Report all metrics
ReportOneScenario <- function(data, vecParams, R, scenarioCount, trueCoefs, alpha) {

    coefData   <- data[, ExtColNames("coef", data)]
    varData    <- data[, ExtColNames("var", data)]
    pvalData   <- data[, ExtColNames("pval", data)]

    meanCoefs  <- AverageCoefs(coefData)
    biases     <- Bias(trueCoefs, meanCoefs)

    meanVars   <- AverageVars(varData)
    trueVars   <- TrueVars(coefData)

    mses       <- Mse(trueVars, biases)
    rejRates   <- RejRate(alpha, pvalData)

    c(scenarioCount = scenarioCount,
      R = R,
      Q = unique(data$Q)[1],
      vecParams,
      meanNk = mean(data$meanNk),
      meanCoefs,
      biases,
      meanVars,
      trueVars,
      mses,
      rejRates)
}
