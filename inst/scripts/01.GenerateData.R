#!/usr/local/bin/Rscript

################################################################################
### Simulation data generation script
##
## Created on: 2016-04-23
## Author: Kazuki Yoshida
################################################################################


###
### Capture data filename argument
################################################################################

## Specify the core count as the first argument
nCores <- as.numeric(commandArgs(trailingOnly = TRUE)[1])

## Execution not allowed without nCores
stopifnot(!is.na(nCores))


### Prepare environment
################################################################################

## sink() if being run non-interactively
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0("./log/", basename(..scriptFileName..), ".txt"), split = TRUE)
}
options(width = 120)

## Record start time
start_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")

## Configure parallelization
## Parallel backend for foreach (also loads foreach and parallel; includes doMC)
library(doParallel)
## Reproducible parallelization
library(doRNG)
## Used by parallel::mclapply() as default
options(mc.cores = nCores)
## Used by doParallel as default
options(cores = nCores)
## Register doParallel as the parallel backend for foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = nCores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Load packages
library(distributed)
library(tidyverse)


cat("
###
### Specify parameter sets
################################################################################\n")

###  Create size set
lstNBase  <- list(10^5,2*10^4,2*10^4,5*10^3)
lstN4_20k <- list(2*10^4,2*10^4,2*10^4,2*10^4)
lstN4_10k <- list(1*10^4,1*10^4,1*10^4,1*10^4)
lstN8     <- list(10^5,2*10^4,2*10^4,5*10^3, 10^5,2*10^4,2*10^4,5*10^3)
lstN4_5k  <- list(5*10^3,5*10^3,5*10^3,5*10^3)
lstN4_1k  <- list(1*10^3,1*10^3,1*10^3,1*10^3)
lstN4_500  <- list(5*10^2,5*10^2,5*10^2,5*10^2)

###  Create covariate generator sets
lstAssignCovariatesNormBinDefault <- list(AssignCovariatesNormBinDefault,
                                          AssignCovariatesNormBinDefault,
                                          AssignCovariatesNormBinDefault,
                                          AssignCovariatesNormBinDefault)
lstAssignCovariatesNormBinSmall <-
    list(ConstructAssignCovariatesNormBin(mean1 = 0.0, sd1 = 1, p2 = 0.05, pLast = 0.50),
         ConstructAssignCovariatesNormBin(mean1 = 0.1, sd1 = 1, p2 = 0.07, pLast = 0.50),
         ConstructAssignCovariatesNormBin(mean1 = 0.2, sd1 = 1, p2 = 0.09, pLast = 0.50),
         ConstructAssignCovariatesNormBin(mean1 = 0.3, sd1 = 1, p2 = 0.11, pLast = 0.50))
lstAssignCovariatesNormBinModerate <-
    list(ConstructAssignCovariatesNormBin(mean1 = 0.0, sd1 = 1, p2 = 0.05, pLast = 0.50),
         ConstructAssignCovariatesNormBin(mean1 = 0.2, sd1 = 1, p2 = 0.07, pLast = 0.48),
         ConstructAssignCovariatesNormBin(mean1 = 0.4, sd1 = 1, p2 = 0.09, pLast = 0.46),
         ConstructAssignCovariatesNormBin(mean1 = 0.6, sd1 = 1, p2 = 0.11, pLast = 0.44))
lstAssignCovariatesNormBinLarge <-
    list(ConstructAssignCovariatesNormBin(mean1 = 0.0, sd1 = 1, p2 = 0.05, pLast = 0.50),
         ConstructAssignCovariatesNormBin(mean1 = 0.3, sd1 = 1, p2 = 0.08, pLast = 0.47),
         ConstructAssignCovariatesNormBin(mean1 = 0.6, sd1 = 1, p2 = 0.11, pLast = 0.44),
         ConstructAssignCovariatesNormBin(mean1 = 0.9, sd1 = 1, p2 = 0.14, pLast = 0.41))

###  Create treatment model coefficient sets
lstAlphasBaseNew <- list(c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)))
## Treatment prevalence changes
alpha0_25perc <- -2.4
alpha0_10perc <- -3.7
alpha0_05perc <- -4.5
lstAlphasTx25 <- list(c(alpha0 = alpha0_25perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_25perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_25perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_25perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)))
lstAlphasTx10 <- list(c(alpha0 = alpha0_10perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_10perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_10perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_10perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)))
lstAlphasTx05 <- list(c(alpha0 = alpha0_05perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_05perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_05perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)),
                      c(alpha0 = alpha0_05perc, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 7)))
## Varying number of confounders
lstAlphasVaryCon <- list(c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 5)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 10)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 20)),
                         c(alpha0 = -0.85, alphaX = 0.7 * seq(log(0.2), log(5), length.out = 40)))
## Null association
lstAlphasNull <- list(c(alpha0 = 0, alphaX = rep(0, 7)),
                      c(alpha0 = 0, alphaX = rep(0, 7)),
                      c(alpha0 = 0, alphaX = rep(0, 7)),
                      c(alpha0 = 0, alphaX = rep(0, 7)))


###  Create outcome model coefficient sets
## Intercept, covariates, treatment, interaction
lstBetasBaseNew <- list(c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                          betaA = 0,
                          betaXA = rep(0,7)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                          betaA = 0,
                          betaXA = rep(0,7)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                          betaA = 0,
                          betaXA = rep(0,7)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                          betaA = 0,
                          betaXA = rep(0,7)))
## Treatment 25% with same confounding for survival
lstBetasTx25 <- list(c(beta0 = -3.5,
                       betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)))
## Treatment 10% with same confounding for survival
lstBetasTx10 <- list(c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)))
## Treatment 5% with same confounding for survival (same as 10%)
lstBetasTx05 <- list(c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)),
                     c(beta0 = -3.5,
                       betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                       betaA = 0,
                       betaXA = rep(0, 7)))
## Rare disease around 1%
lstBetasRareDis1 <- list(c(beta0 = -5.5,
                           betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                           betaA = 0,
                           betaXA = rep(0,7)),
                         c(beta0 = -5.5,
                           betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                           betaA = 0,
                           betaXA = rep(0,7)),
                         c(beta0 = -5.5,
                           betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                           betaA = 0,
                           betaXA = rep(0,7)),
                         c(beta0 = -5.5,
                           betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                           betaA = 0,
                           betaXA = rep(0,7)))
## Rare disease around 0.1%
lstBetasRareDis.1 <- list(c(beta0 = -7.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = 0,
                            betaXA = rep(0,7)),
                          c(beta0 = -7.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = 0,
                            betaXA = rep(0,7)),
                          c(beta0 = -7.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = 0,
                            betaXA = rep(0,7)),
                          c(beta0 = -7.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = 0,
                            betaXA = rep(0,7)))
## Rare disease around 0.01%
lstBetasRareDis.01 <- list(c(beta0 = -10.0,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -10.0,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -10.0,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -10.0,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)))
## Rare disease around 1% (cHR 1.2)
lstBetasRareDis1HR1.2 <- list(c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0,7)))
## Rare disease around 0.1% (cHR 1.2)
lstBetasRareDis.1HR1.2 <- list(c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(1.2),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(1.2),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(1.2),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(1.2),
                                 betaXA = rep(0,7)))
## Rare disease around 0.01% (cHR 1.2)
lstBetasRareDis.01HR1.2 <- list(c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(1.2),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(1.2),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(1.2),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(1.2),
                                  betaXA = rep(0,7)))
## Rare disease around 1% (cHR 0.8)
lstBetasRareDis1HR0.8 <- list(c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0,7)),
                              c(beta0 = -5.5,
                                betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0,7)))
## Rare disease around 0.1% (cHR 0.8)
lstBetasRareDis.1HR0.8 <- list(c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(0.8),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(0.8),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(0.8),
                                 betaXA = rep(0,7)),
                               c(beta0 = -7.5,
                                 betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                 betaA = log(0.8),
                                 betaXA = rep(0,7)))
## Rare disease around 0.01% (cHR 0.8)
lstBetasRareDis.01HR0.8 <- list(c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(0.8),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(0.8),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(0.8),
                                  betaXA = rep(0,7)),
                                c(beta0 = -10.0,
                                  betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                                  betaA = log(0.8),
                                  betaXA = rep(0,7)))
## Vary disease incidence
lstBetasVaryDisInc <- list(c(beta0 = -3.5,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -5.5,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -7.5,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)),
                           c(beta0 = -10.0,
                             betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                             betaA = 0,
                             betaXA = rep(0,7)))
## Harmful treatment
lstBetasHarmHR1.2 <- list(c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(1.2),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(1.2),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(1.2),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(1.2),
                            betaXA = rep(0,7)))
lstBetasHarmHR2.0 <- list(c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(2.0),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(2.0),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(2.0),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(2.0),
                            betaXA = rep(0,7)))
## Harmful treatment 25% with same confounding for survival
lstBetasHarmHR1.2Tx25 <- list(c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)))
## Harmful treatment 10% with same confounding for survival
lstBetasHarmHR1.2Tx10 <- list(c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)))
## Harmful treatment 5% with same confounding for survival (same as 10%)
lstBetasHarmHR1.2Tx05 <- list(c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(1.2),
                                betaXA = rep(0, 7)))
## Protective treatment
lstBetasProtHR0.8 <- list(c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(0.8),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(0.8),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(0.8),
                            betaXA = rep(0,7)),
                          c(beta0 = -3.5,
                            betaX = 0.5 * seq(log(0.2), log(5), length.out = 7),
                            betaA = log(0.8),
                            betaXA = rep(0,7)))
## Protective treatment 25% with same confounding for survival
lstBetasProtHR0.8Tx25 <- list(c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.48 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)))
## Protective treatment 10% with same confounding for survival
lstBetasProtHR0.8Tx10 <- list(c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)))
## Protective treatment 5% with same confounding for survival (same as 10%)
lstBetasProtHR0.8Tx05 <- list(c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)),
                              c(beta0 = -3.5,
                                betaX = 0.46 * seq(log(0.2), log(5), length.out = 7),
                                betaA = log(0.8),
                                betaXA = rep(0, 7)))
## Varying number of confounders
lstBetasVaryCon <- list(c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 5),
                          betaA = 0,
                          betaXA = rep(0,5)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 10),
                          betaA = 0,
                          betaXA = rep(0,10)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 20),
                          betaA = 0,
                          betaXA = rep(0,20)),
                        c(beta0 = -3.5,
                          betaX = 0.5 * seq(log(0.2), log(5), length.out = 40),
                          betaA = 0,
                          betaXA = rep(0,40)))
## Null covariates, harmful treatment (HR2.0)
lstBetasNullHR2.0 <- list(c(beta0 = -3.5,
                            betaX = rep(0, 7),
                            betaA = log(2.0),
                            betaXA = rep(0, 7)),
                          c(beta0 = -3.5,
                            betaX = rep(0, 7),
                            betaA = log(2.0),
                            betaXA = rep(0, 7)),
                          c(beta0 = -3.5,
                            betaX = rep(0, 7),
                            betaA = log(2.0),
                            betaXA = rep(0, 7)),
                          c(beta0 = -3.5,
                            betaX = rep(0, 7),
                            betaA = log(2.0),
                            betaXA = rep(0, 7)))

###  Create survival outcome model parameter sets
##
lstSurvParamsBaseNew <- list(c(-log(0.97), -log(0.9), 1),
                             c(-log(0.97), -log(0.9), 1),
                             c(-log(0.97), -log(0.9), 1),
                             c(-log(0.97), -log(0.9), 1))
##
lstSurvParamsTx25 <- list(c(-log(0.975), -log(0.91), 1),
                          c(-log(0.975), -log(0.91), 1),
                          c(-log(0.975), -log(0.91), 1),
                          c(-log(0.975), -log(0.91), 1))
##
lstSurvParamsTx10 <- list(c(-log(0.978), -log(0.91), 1),
                          c(-log(0.978), -log(0.91), 1),
                          c(-log(0.978), -log(0.91), 1),
                          c(-log(0.978), -log(0.91), 1))
##
lstSurvParamsTx05 <- list(c(-log(0.98), -log(0.915), 1),
                          c(-log(0.98), -log(0.915), 1),
                          c(-log(0.98), -log(0.915), 1),
                          c(-log(0.98), -log(0.915), 1))
##
lstSurvParamsRareDis1 <- list(c(-log(0.995), -log(0.9), 1),
                              c(-log(0.995), -log(0.9), 1),
                              c(-log(0.995), -log(0.9), 1),
                              c(-log(0.995), -log(0.9), 1))
##
lstSurvParamsRareDis.1 <- list(c(-log(0.9995), -log(0.9), 1),
                               c(-log(0.9995), -log(0.9), 1),
                               c(-log(0.9995), -log(0.9), 1),
                               c(-log(0.9995), -log(0.9), 1))
##
lstSurvParamsRareDis.01 <- list(c(-log(0.99995), -log(0.9), 1),
                                c(-log(0.99995), -log(0.9), 1),
                                c(-log(0.99995), -log(0.9), 1),
                                c(-log(0.99995), -log(0.9), 1))
##
lstSurvParamsVaryDisInc <- list(c(-log(0.97), -log(0.9), 1),
                                c(-log(0.995), -log(0.9), 1),
                                c(-log(0.9995), -log(0.9), 1),
                                c(-log(0.99995), -log(0.9), 1))

cat("
###
### Create scenarios
################################################################################\n")

## Actually generate scenarios
Scenarios <- GenerateScenarios(
    lstLstN          = list("base" = lstNBase,
                            "tx25" = lstNBase,
                            "tx10" = lstNBase,
                            "tx05" = lstNBase,
                            "dis1perc" = lstNBase,
                            "dis0.1perc" = lstNBase,
                            "dis0.01perc" = lstNBase,
                            "dis_var_inc" = lstN4_20k,
                            "harm_hr1.2" = lstNBase,
                            "harm_hr2.0" = lstNBase,
                            "2_5k_sites" = lstN4_5k[rep(1,2)],
                            "4_5k_sites" = lstN4_5k[rep(1,4)],
                            "8_5k_sites" = lstN4_5k[rep(1,8)],
                            "vary_conf_count" = lstN4_20k,
                            "noise_cov_hr2.0" = lstNBase,
                            "instruments_hr2.0" = lstNBase,
                            "rct_hr2.0" = lstNBase,
                            "new_tx25" = lstNBase,
                            "new_tx10" = lstNBase,
                            "new_tx05" = lstNBase,
                            "16_5k_sites" = lstN4_5k[rep(1,16)],
                            "new_harm_hr1.2_tx25" = lstNBase,
                            "new_harm_hr1.2_tx10" = lstNBase,
                            "new_harm_hr1.2_tx05" = lstNBase,
                            "new_prot_hr0.8" = lstNBase,
                            "new_prot_hr0.8_tx25" = lstNBase,
                            "new_prot_hr0.8_tx10" = lstNBase,
                            "new_prot_hr0.8_tx05" = lstNBase,
                            "new_dis1perc_hr1.2" = lstNBase,
                            "new_dis0.1perc_hr1.2" = lstNBase,
                            "new_dis0.01perc_hr1.2" = lstNBase,
                            "new_dis1perc_hr0.8" = lstNBase,
                            "new_dis0.1perc_hr0.8" = lstNBase,
                            "new_dis0.01perc_hr0.8" = lstNBase,
                            "2_10k_sites" = lstN4_10k[rep(1,2)],
                            "4_10k_sites" = lstN4_10k[rep(1,4)],
                            "8_10k_sites" = lstN4_10k[rep(1,8)],
                            "16_10k_sites" = lstN4_10k[rep(1,16)],
                            "2_20k_sites" = lstN4_20k[rep(1,2)],
                            "4_20k_sites" = lstN4_20k[rep(1,4)],
                            "8_20k_sites" = lstN4_20k[rep(1,8)],
                            "16_20k_sites" = lstN4_20k[rep(1,16)],
                            "2_1k_sites" = lstN4_1k[rep(1,2)],
                            "4_1k_sites" = lstN4_1k[rep(1,4)],
                            "8_1k_sites" = lstN4_1k[rep(1,8)],
                            "16_1k_sites" = lstN4_1k[rep(1,16)],
                            "2_500_sites" = lstN4_500[rep(1,2)],
                            "4_500_sites" = lstN4_500[rep(1,4)],
                            "8_500_sites" = lstN4_500[rep(1,8)],
                            "16_500_sites" = lstN4_500[rep(1,16)],
                            "covariates_small_hr0.8" = lstNBase,
                            "covariates_small_hr1.0" = lstNBase,
                            "covariates_small_hr1.2" = lstNBase,
                            "covariates_moderate_hr0.8" = lstNBase,
                            "covariates_moderate_hr1.0" = lstNBase,
                            "covariates_moderate_hr1.2" = lstNBase,
                            "covariates_large_hr0.8" = lstNBase,
                            "covariates_large_hr1.0" = lstNBase,
                            "covariates_large_hr1.2" = lstNBase),
    ##
    lstLstAssignCovariates = list("base" = lstAssignCovariatesNormBinDefault,
                                  "tx25" = lstAssignCovariatesNormBinDefault,
                                  "tx10" = lstAssignCovariatesNormBinDefault,
                                  "tx05" = lstAssignCovariatesNormBinDefault,
                                  "dis1perc" = lstAssignCovariatesNormBinDefault,
                                  "dis0.1perc" = lstAssignCovariatesNormBinDefault,
                                  "dis0.01perc" = lstAssignCovariatesNormBinDefault,
                                  "dis_var_inc" = lstAssignCovariatesNormBinDefault,
                                  "harm_hr1.2" = lstAssignCovariatesNormBinDefault,
                                  "harm_hr2.0" = lstAssignCovariatesNormBinDefault,
                                  "2_5k_sites" = lstAssignCovariatesNormBinDefault[rep(1,2)],
                                  "4_5k_sites" = lstAssignCovariatesNormBinDefault[rep(1,4)],
                                  "8_5k_sites" = lstAssignCovariatesNormBinDefault[rep(1,8)],
                                  "vary_conf_count" = lstAssignCovariatesNormBinDefault,
                                  "noise_cov_hr2.0" = lstAssignCovariatesNormBinDefault,
                                  "instruments_hr2.0" = lstAssignCovariatesNormBinDefault,
                                  "rct_hr2.0" = lstAssignCovariatesNormBinDefault,
                                  "new_tx25" = lstAssignCovariatesNormBinDefault,
                                  "new_tx10" = lstAssignCovariatesNormBinDefault,
                                  "new_tx05" = lstAssignCovariatesNormBinDefault,
                                  "16_5k_sites" = lstAssignCovariatesNormBinDefault[rep(1,16)],
                                  "new_harm_hr1.2_tx25" = lstAssignCovariatesNormBinDefault,
                                  "new_harm_hr1.2_tx10" = lstAssignCovariatesNormBinDefault,
                                  "new_harm_hr1.2_tx05" = lstAssignCovariatesNormBinDefault,
                                  "new_prot_hr0.8" = lstAssignCovariatesNormBinDefault,
                                  "new_prot_hr0.8_tx25" = lstAssignCovariatesNormBinDefault,
                                  "new_prot_hr0.8_tx10" = lstAssignCovariatesNormBinDefault,
                                  "new_prot_hr0.8_tx05" = lstAssignCovariatesNormBinDefault,
                                  "new_dis1perc_hr1.2" = lstAssignCovariatesNormBinDefault,
                                  "new_dis0.1perc_hr1.2" = lstAssignCovariatesNormBinDefault,
                                  "new_dis0.01perc_hr1.2" = lstAssignCovariatesNormBinDefault,
                                  "new_dis1perc_hr0.8" = lstAssignCovariatesNormBinDefault,
                                  "new_dis0.1perc_hr0.8" = lstAssignCovariatesNormBinDefault,
                                  "new_dis0.01perc_hr0.8" = lstAssignCovariatesNormBinDefault,
                                  "2_10k_sites" = lstAssignCovariatesNormBinDefault[rep(1,2)],
                                  "4_10k_sites" = lstAssignCovariatesNormBinDefault[rep(1,4)],
                                  "8_10k_sites" = lstAssignCovariatesNormBinDefault[rep(1,8)],
                                  "16_10k_sites" = lstAssignCovariatesNormBinDefault[rep(1,16)],
                                  "2_20k_sites" = lstAssignCovariatesNormBinDefault[rep(1,2)],
                                  "4_20k_sites" = lstAssignCovariatesNormBinDefault[rep(1,4)],
                                  "8_20k_sites" = lstAssignCovariatesNormBinDefault[rep(1,8)],
                                  "16_20k_sites" = lstAssignCovariatesNormBinDefault[rep(1,16)],
                                  "2_1k_sites" = lstAssignCovariatesNormBinDefault[rep(1,2)],
                                  "4_1k_sites" = lstAssignCovariatesNormBinDefault[rep(1,4)],
                                  "8_1k_sites" = lstAssignCovariatesNormBinDefault[rep(1,8)],
                                  "16_1k_sites" = lstAssignCovariatesNormBinDefault[rep(1,16)],
                                  "2_500_sites" = lstAssignCovariatesNormBinDefault[rep(1,2)],
                                  "4_500_sites" = lstAssignCovariatesNormBinDefault[rep(1,4)],
                                  "8_500_sites" = lstAssignCovariatesNormBinDefault[rep(1,8)],
                                  "16_500_sites" = lstAssignCovariatesNormBinDefault[rep(1,16)],
                                  "covariates_small_hr0.8" = lstAssignCovariatesNormBinSmall,
                                  "covariates_small_hr1.0" = lstAssignCovariatesNormBinSmall,
                                  "covariates_small_hr1.2" = lstAssignCovariatesNormBinSmall,
                                  "covariates_moderate_hr0.8" = lstAssignCovariatesNormBinModerate,
                                  "covariates_moderate_hr1.0" = lstAssignCovariatesNormBinModerate,
                                  "covariates_moderate_hr1.2" = lstAssignCovariatesNormBinModerate,
                                  "covariates_large_hr0.8" = lstAssignCovariatesNormBinLarge,
                                  "covariates_large_hr1.0" = lstAssignCovariatesNormBinLarge,
                                  "covariates_large_hr1.2" = lstAssignCovariatesNormBinLarge),
    ##
    lstLstAlphas     = list("base" = lstAlphasBaseNew,
                            "tx25" = lstAlphasTx25,
                            "tx10" = lstAlphasTx10,
                            "tx05" = lstAlphasTx05,
                            "dis1perc" = lstAlphasBaseNew,
                            "dis0.1perc" = lstAlphasBaseNew,
                            "dis0.01perc" = lstAlphasBaseNew,
                            "dis_var_inc" = lstAlphasBaseNew,
                            "harm_hr1.2" = lstAlphasBaseNew,
                            "harm_hr2.0" = lstAlphasBaseNew,
                            "2_5k_sites" = lstAlphasBaseNew[rep(1,2)],
                            "4_5k_sites" = lstAlphasBaseNew[rep(1,4)],
                            "8_5k_sites" = lstAlphasBaseNew[rep(1,8)],
                            "vary_conf_count" = lstAlphasVaryCon,
                            "noise_cov_hr2.0" = lstAlphasNull,
                            "instruments_hr2.0" = lstAlphasBaseNew,
                            "rct_hr2.0" = lstAlphasNull,
                            "new_tx25" = lstAlphasTx25,
                            "new_tx10" = lstAlphasTx10,
                            "new_tx05" = lstAlphasTx05,
                            "16_5k_sites" = lstAlphasBaseNew[rep(1,16)],
                            "new_harm_hr1.2_tx25" = lstAlphasTx25,
                            "new_harm_hr1.2_tx10" = lstAlphasTx10,
                            "new_harm_hr1.2_tx05" = lstAlphasTx05,
                            "new_prot_hr0.8" = lstAlphasBaseNew,
                            "new_prot_hr0.8_tx25" = lstAlphasTx25,
                            "new_prot_hr0.8_tx10" = lstAlphasTx10,
                            "new_prot_hr0.8_tx05" = lstAlphasTx05,
                            "new_dis1perc_hr1.2" = lstAlphasBaseNew,
                            "new_dis0.1perc_hr1.2" = lstAlphasBaseNew,
                            "new_dis0.01perc_hr1.2" = lstAlphasBaseNew,
                            "new_dis1perc_hr0.8" = lstAlphasBaseNew,
                            "new_dis0.1perc_hr0.8" = lstAlphasBaseNew,
                            "new_dis0.01perc_hr0.8" = lstAlphasBaseNew,
                            "2_10k_sites" = lstAlphasBaseNew[rep(1,2)],
                            "4_10k_sites" = lstAlphasBaseNew[rep(1,4)],
                            "8_10k_sites" = lstAlphasBaseNew[rep(1,8)],
                            "16_10k_sites" = lstAlphasBaseNew[rep(1,16)],
                            "2_20k_sites" = lstAlphasBaseNew[rep(1,2)],
                            "4_20k_sites" = lstAlphasBaseNew[rep(1,4)],
                            "8_20k_sites" = lstAlphasBaseNew[rep(1,8)],
                            "16_20k_sites" = lstAlphasBaseNew[rep(1,16)],
                            "2_1k_sites" = lstAlphasBaseNew[rep(1,2)],
                            "4_1k_sites" = lstAlphasBaseNew[rep(1,4)],
                            "8_1k_sites" = lstAlphasBaseNew[rep(1,8)],
                            "16_1k_sites" = lstAlphasBaseNew[rep(1,16)],
                            "2_500_sites" = lstAlphasBaseNew[rep(1,2)],
                            "4_500_sites" = lstAlphasBaseNew[rep(1,4)],
                            "8_500_sites" = lstAlphasBaseNew[rep(1,8)],
                            "16_500_sites" = lstAlphasBaseNew[rep(1,16)],
                            "covariates_small_hr0.8" = lstAlphasBaseNew,
                            "covariates_small_hr1.0" = lstAlphasBaseNew,
                            "covariates_small_hr1.2" = lstAlphasBaseNew,
                            "covariates_moderate_hr0.8" = lstAlphasBaseNew,
                            "covariates_moderate_hr1.0" = lstAlphasBaseNew,
                            "covariates_moderate_hr1.2" = lstAlphasBaseNew,
                            "covariates_large_hr0.8" = lstAlphasBaseNew,
                            "covariates_large_hr1.0" = lstAlphasBaseNew,
                            "covariates_large_hr1.2" = lstAlphasBaseNew),
    ##
    lstLstBetas      = list("base" = lstBetasBaseNew,
                            "tx25" = lstBetasBaseNew,
                            "tx10" = lstBetasBaseNew,
                            "tx05" = lstBetasBaseNew,
                            "dis1perc" = lstBetasRareDis1,
                            "dis0.1perc" = lstBetasRareDis.1,
                            "dis0.01perc" = lstBetasRareDis.01,
                            "dis_var_inc" = lstBetasVaryDisInc,
                            "harm_hr1.2" = lstBetasHarmHR1.2,
                            "harm_hr2.0" = lstBetasHarmHR2.0,
                            "2_5k_sites" = lstBetasBaseNew[rep(1,2)],
                            "4_5k_sites" = lstBetasBaseNew[rep(1,4)],
                            "8_5k_sites" = lstBetasBaseNew[rep(1,8)],
                            "vary_conf_count" = lstBetasVaryCon,
                            "noise_cov_hr2.0" = lstBetasNullHR2.0,
                            "instruments_hr2.0" = lstBetasNullHR2.0,
                            "rct_hr2.0" = lstBetasHarmHR2.0,
                            "new_tx25" = lstBetasTx25,
                            "new_tx10" = lstBetasTx10,
                            "new_tx05" = lstBetasTx05,
                            "16_5k_sites" = lstBetasBaseNew[rep(1,16)],
                            "new_harm_hr1.2_tx25" = lstBetasHarmHR1.2Tx25,
                            "new_harm_hr1.2_tx10" = lstBetasHarmHR1.2Tx10,
                            "new_harm_hr1.2_tx05" = lstBetasHarmHR1.2Tx05,
                            "new_prot_hr0.8" = lstBetasProtHR0.8,
                            "new_prot_hr0.8_tx25" = lstBetasProtHR0.8Tx25,
                            "new_prot_hr0.8_tx10" = lstBetasProtHR0.8Tx10,
                            "new_prot_hr0.8_tx05" = lstBetasProtHR0.8Tx05,
                            "new_dis1perc_hr1.2" = lstBetasRareDis1HR1.2,
                            "new_dis0.1perc_hr1.2" = lstBetasRareDis.1HR1.2,
                            "new_dis0.01perc_hr1.2" = lstBetasRareDis.01HR1.2,
                            "new_dis1perc_hr0.8" = lstBetasRareDis1HR0.8,
                            "new_dis0.1perc_hr0.8" = lstBetasRareDis.1HR0.8,
                            "new_dis0.01perc_hr0.8" = lstBetasRareDis.01HR0.8,
                            "2_10k_sites" = lstBetasBaseNew[rep(1,2)],
                            "4_10k_sites" = lstBetasBaseNew[rep(1,4)],
                            "8_10k_sites" = lstBetasBaseNew[rep(1,8)],
                            "16_10k_sites" = lstBetasBaseNew[rep(1,16)],
                            "2_20k_sites" = lstBetasBaseNew[rep(1,2)],
                            "4_20k_sites" = lstBetasBaseNew[rep(1,4)],
                            "8_20k_sites" = lstBetasBaseNew[rep(1,8)],
                            "16_20k_sites" = lstBetasBaseNew[rep(1,16)],
                            "2_1k_sites" = lstBetasBaseNew[rep(1,2)],
                            "4_1k_sites" = lstBetasBaseNew[rep(1,4)],
                            "8_1k_sites" = lstBetasBaseNew[rep(1,8)],
                            "16_1k_sites" = lstBetasBaseNew[rep(1,16)],
                            "2_500_sites" = lstBetasBaseNew[rep(1,2)],
                            "4_500_sites" = lstBetasBaseNew[rep(1,4)],
                            "8_500_sites" = lstBetasBaseNew[rep(1,8)],
                            "16_500_sites" = lstBetasBaseNew[rep(1,16)],
                            "covariates_small_hr0.8" = lstBetasProtHR0.8,
                            "covariates_small_hr1.0" = lstBetasBaseNew,
                            "covariates_small_hr1.2" = lstBetasHarmHR1.2,
                            "covariates_moderate_hr0.8" = lstBetasProtHR0.8,
                            "covariates_moderate_hr1.0" = lstBetasBaseNew,
                            "covariates_moderate_hr1.2" = lstBetasHarmHR1.2,
                            "covariates_large_hr0.8" = lstBetasProtHR0.8,
                            "covariates_large_hr1.0" = lstBetasBaseNew,
                            "covariates_large_hr1.2" = lstBetasHarmHR1.2),
    ##
    lstLstSurvParams = list("base" = lstSurvParamsBaseNew,
                            "tx25" = lstSurvParamsBaseNew,
                            "tx10" = lstSurvParamsBaseNew,
                            "tx05" = lstSurvParamsBaseNew,
                            "dis1perc" = lstSurvParamsRareDis1,
                            "dis0.1perc" = lstSurvParamsRareDis.1,
                            "dis0.01perc" = lstSurvParamsRareDis.01,
                            "dis_var_inc" = lstSurvParamsVaryDisInc,
                            "harm_hr1.2" = lstSurvParamsBaseNew,
                            "harm_hr2.0" = lstSurvParamsBaseNew,
                            "2_5k_sites" = lstSurvParamsBaseNew[rep(1,2)],
                            "4_5k_sites" = lstSurvParamsBaseNew[rep(1,4)],
                            "8_5k_sites" = lstSurvParamsBaseNew[rep(1,8)],
                            "vary_conf_count" = lstSurvParamsBaseNew,
                            "noise_cov_hr2.0" = lstSurvParamsBaseNew,
                            "instruments_hr2.0" = lstSurvParamsBaseNew,
                            "rct_hr2.0" = lstSurvParamsBaseNew,
                            "new_tx25" = lstSurvParamsTx25,
                            "new_tx10" = lstSurvParamsTx10,
                            "new_tx05" = lstSurvParamsTx05,
                            "16_5k_sites" = lstSurvParamsBaseNew[rep(1,16)],
                            "new_harm_hr1.2_tx25" = lstSurvParamsTx25,
                            "new_harm_hr1.2_tx10" = lstSurvParamsTx10,
                            "new_harm_hr1.2_tx05" = lstSurvParamsTx05,
                            "new_prot_hr0.8" = lstSurvParamsBaseNew,
                            "new_prot_hr0.8_tx25" = lstSurvParamsTx25,
                            "new_prot_hr0.8_tx10" = lstSurvParamsTx10,
                            "new_prot_hr0.8_tx05" = lstSurvParamsTx05,
                            "new_dis1perc_hr1.2" = lstSurvParamsRareDis1,
                            "new_dis0.1perc_hr1.2" = lstSurvParamsRareDis.1,
                            "new_dis0.01perc_hr1.2" = lstSurvParamsRareDis.01,
                            "new_dis1perc_hr0.8" = lstSurvParamsRareDis1,
                            "new_dis0.1perc_hr0.8" = lstSurvParamsRareDis.1,
                            "new_dis0.01perc_hr0.8" = lstSurvParamsRareDis.01,
                            "2_10k_sites" = lstSurvParamsBaseNew[rep(1,2)],
                            "4_10k_sites" = lstSurvParamsBaseNew[rep(1,4)],
                            "8_10k_sites" = lstSurvParamsBaseNew[rep(1,8)],
                            "16_10k_sites" = lstSurvParamsBaseNew[rep(1,16)],
                            "2_20k_sites" = lstSurvParamsBaseNew[rep(1,2)],
                            "4_20k_sites" = lstSurvParamsBaseNew[rep(1,4)],
                            "8_20k_sites" = lstSurvParamsBaseNew[rep(1,8)],
                            "16_20k_sites" = lstSurvParamsBaseNew[rep(1,16)],
                            "2_1k_sites" = lstSurvParamsBaseNew[rep(1,2)],
                            "4_1k_sites" = lstSurvParamsBaseNew[rep(1,4)],
                            "8_1k_sites" = lstSurvParamsBaseNew[rep(1,8)],
                            "16_1k_sites" = lstSurvParamsBaseNew[rep(1,16)],
                            "2_500_sites" = lstSurvParamsBaseNew[rep(1,2)],
                            "4_500_sites" = lstSurvParamsBaseNew[rep(1,4)],
                            "8_500_sites" = lstSurvParamsBaseNew[rep(1,8)],
                            "16_500_sites" = lstSurvParamsBaseNew[rep(1,16)],
                            "covariates_small_hr0.8" = lstSurvParamsBaseNew,
                            "covariates_small_hr1.0" = lstSurvParamsBaseNew,
                            "covariates_small_hr1.2" = lstSurvParamsBaseNew,
                            "covariates_moderate_hr0.8" = lstSurvParamsBaseNew,
                            "covariates_moderate_hr1.0" = lstSurvParamsBaseNew,
                            "covariates_moderate_hr1.2" = lstSurvParamsBaseNew,
                            "covariates_large_hr0.8" = lstSurvParamsBaseNew,
                            "covariates_large_hr1.0" = lstSurvParamsBaseNew,
                            "covariates_large_hr1.2" = lstSurvParamsBaseNew),
    mix = FALSE)
names(Scenarios) <- c(
    "base 50% treatment 5% incidence 4 sites",
    "treatment 25%",
    "treatment 10%",
    "treatment  5%",
    "disease 1%",
    "disease 0.1%",
    "disease 0.01%",
    "varying disease prevalence (5, 1, 0.1, 0.01)",
    "harmful treatment HR 1.2",
    "harmful treatment HR 2.0",
    "2 5k sites, 5% incidence",
    "4 5k sites, 5% incidence",
    "8 5k sites, 5% incidence",
    "varying confounder counts (5, 10, 20, 40)",
    "Noise covariates, harmful treatment log(2.0)",
    "Instrument only, harmful treatment log(2.0)",
    "RCT, harmful treatment log(2.0)",
    ## These maintain the same level of confounding in infrequent treatment.
    "new treatment 25%",
    "new treatment 10%",
    "new treatment  5%",
    ##
    "16 5k sites, 5% incidence",
    "25% harmful treatment HR 1.2",
    "10% harmful treatment HR 1.2",
    "5% harmful treatment HR 1.2",
    "protective treatment HR 0.8",
    "25% protective treatment HR 0.8",
    "10% protective treatment HR 0.8",
    "5% protective treatment HR 0.8",
    "harmful treatment HR 1.2, disease 1%",
    "harmful treatment HR 1.2, disease 0.1%",
    "harmful treatment HR 1.2, disease 0.01%",
    "protective treatment HR 0.8, disease 1%",
    "protective treatment HR 0.8, disease 0.1%",
    "protective treatment HR 0.8, disease 0.01%",
    ##
    "2 10k sites, 5% incidence",
    "4 10k sites, 5% incidence",
    "8 10k sites, 5% incidence",
    "16 10k sites, 5% incidence",
    "2 20k sites, 5% incidence",
    "4 20k sites, 5% incidence",
    "8 20k sites, 5% incidence",
    "16 20k sites, 5% incidence",
    ##
    "2 1k sites, 5% incidence",
    "4 1k sites, 5% incidence",
    "8 1k sites, 5% incidence",
    "16 1k sites, 5% incidence",
    "2 500 sites, 5% incidence",
    "4 500 sites, 5% incidence",
    "8 500 sites, 5% incidence",
    "16 500 sites, 5% incidence",
    ##
    "Covariates small variation, HR 0.8",
    "Covariates small variation, HR 1.0",
    "Covariates small variation, HR 1.2",
    "Covariates moderate variation, HR 0.8",
    "Covariates moderate variation, HR 1.0",
    "Covariates moderate variation, HR 1.2",
    "Covariates large variation, HR 0.8",
    "Covariates large variation, HR 1.0",
    "Covariates large variation, HR 1.2"
)

## Use this to invalidate some scenarios to make them placeholders
## GenerateDataForAllScenarios() skips invalid ones.
if (TRUE) {
    for (i in which(!grepl("Covariates", names(Scenarios)))) {
        Scenarios[i] <- list(NULL)
        cat("Dropped", names(Scenarios[i]), "\n")
    }

    for (i in seq_along(Scenarios)) {
        if (!is.null(Scenarios[[i]])) {
            cat("Kept", names(Scenarios[i]), "\n")
        }
    }
}

print(Scenarios)

cat("###  Use for script running\n")
paste(paste0("./data/ScenarioRaw0", which(grepl("Covariates", names(Scenarios))), "*"), collapse = " ")
paste(paste0("./data/ScenarioPrepared0", which(grepl("Covariates", names(Scenarios))), "*"), collapse = " ")


cat("
###
### Generate datasets
################################################################################\n")

set.seed(201605130)
parts <- 10
R <- 50

## Check for the data subfolder if not available
if (!("data" %in% dir())) {
    stop("data subfolder is not available. Create one.")
}
## Move to the data subfolder
setwd("./data")

## Parallelization at the scenario x parts level.
## If there are fewer scenario x parts combos than cores requested, CPU time is wasted.
GenerateDataForAllScenarios(Scenarios = Scenarios,
                            parts = parts,
                            R = R)


################################################################################
cat("\n### Record package versions\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("### Started ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
