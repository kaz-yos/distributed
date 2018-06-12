#!/usr/local/bin/Rscript

################################################################################
### Prepare datasets for further analyses
##
## Created on: 2016-05-27
## Author: Kazuki Yoshida
################################################################################


###
### Capture data filename argument
################################################################################

## Specify data file as the first argument
dataFileName <- commandArgs(trailingOnly = TRUE)[1]
## Specify the core count as the second argument
nCores <- as.numeric(commandArgs(trailingOnly = TRUE)[2])

## Execution not allowed without data file
stopifnot(!is.na(dataFileName))
## Execution not allowed without nCores
stopifnot(!is.na(nCores))

## Check it is a scenario file
if (!grepl("Scenario", dataFileName)) {
    stop("Not a scenario file")
}


###
### Prepare environment
################################################################################

## sink() if being run non-interactively
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    ## Include data file name
    sink(file = paste0("./log/",
                       basename(..scriptFileName..),
                       "_",
                       ## Data file name without folder or extension
                       tools::file_path_sans_ext(basename(dataFileName)),
                       ".txt"), split = TRUE)
    options(width = 100)
}

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
library(tidyverse)
library(distributed)

set.seed(201607140)


cat("
###
### Load data file and prepare
################################################################################\n")

## Load
load(dataFileName)

## Prepare data for analysis readiness
lstIterReady <- foreach::foreach(i = seq_along(lstIter)) %dorng% {
    cat("### Iteration", i, "\n")
    out <- try(RequestSiteDataPreparation(lstIter[[i]]))
    ## Report failure
    if (is.error(out)) {
        cat("### Iteration", i, "failure \n")
        print(out)
    }
    out
}

## New file name
newDataFileName <- gsub("ScenarioRaw", "ScenarioPrepared", dataFileName)

## Save as a new file
save(ScenarioDistResNet, scenarioCount, partCount, R, lstIterReady,
     file = newDataFileName)


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("### Started  ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
