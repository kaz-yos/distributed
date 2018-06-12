#!/usr/local/bin/Rscript

################################################################################
### Process analysis file and aggregate into one file
##
## Created on: 2018-03-19
## Author: Kazuki Yoshida
################################################################################


###
### Prepare environment
################################################################################

## sink() if being run non-interactively
if (sink.number() != 0) {sink()}
.script_name. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(.script_name.) == 1) {
    sink(file = paste0(.script_name., ".txt"), split = TRUE)
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
## Detect core count
n_cores <- min(parallel::detectCores(), 8)
## Used by parallel::mclapply() as default
options(mc.cores = n_cores)
## Used by doParallel as default
options(cores = n_cores)
## Register doParallel as the parallel backend for foreach
## http://stackoverflow.com/questions/28989855/the-difference-between-domc-and-doparallel-in-r
doParallel::registerDoParallel(cores = n_cores)
## Report multicore use
cat("### Using", foreach::getDoParWorkers(), "cores\n")
cat("### Using", foreach::getDoParName(), "as backend\n")

## Load packages
library(tidyverse)
library(distributed)
library(openxlsx)


cat("
###
### Load data
################################################################################\n")

dirName <- "./data/"
analysisFileNames <- paste0(dirName,
                            Filter(f = function(elt) {
                                grepl(".*ScenarioAnalyzed.*", elt)
                            }, x = dir(dirName)))


### Loop over all analysis result files and load them into a list
lstAnalysis <- lapply(analysisFileNames, function(file) {

    ## Load file corresponding to a part within a scenario
    load(file)

    ## Binary indicator for valid iterations
    logi_valid_iter <- !sapply(lstIterAnalyzed, (distributed:::is.error))
    ## Valid iteration numbers
    valid_iter <- seq_along(lstIterAnalyzed)[logi_valid_iter]

    ## report results
    cat("## Working on", file, "\n")
    cat("##  Valid iterations:", sum(logi_valid_iter), "out of", length(logi_valid_iter), "\n")
    if (any(!logi_valid_iter)) {
        cat("##  Failed iterations are\n")
        print(which(!logi_valid_iter))
        cat("##  Failure reasons are\n")
        print(unique(lstIterAnalyzed[sapply(lstIterAnalyzed, (distributed:::is.error))]))
    }

    ## Loop over iterations within a part
    lst <- lapply(valid_iter, function(i) {

        cbind(scenario = scenarioCount,
              part = partCount,
              iter = i,
              ## Iteration-level df
              lstIterAnalyzed[[i]],
              stringsAsFactors = FALSE)
    })

    ## Combine as df at part level
    df <- dplyr::bind_rows(lst)

    ## Add scenario as an attribute
    attr(df, "ScenarioDistResNet") <- ScenarioDistResNet

    df
})

cat("
###  Show first file, first part data\n")
head(head(lstAnalysis, 1)[[1]])
cat("
###  Show last file, last part data\n")
tail(tail(lstAnalysis, 1)[[1]])


### Create a data frame of all analysis results
dfAnalysis <- dplyr::bind_rows(lstAnalysis) %>%
    as_data_frame
dfAnalysis

###  Split "method" variable into several metrics
dfAnalysis <- dfAnalysis %>%
    mutate(method = as.character(dfAnalysis$method),
           ## Summary score
           score = SeqGsub(method,
                           c(".*ePs.*", "PS",
                             ".*eDrs.*", "DRS",
                             "^unadj.*", "None")),
           score = factor(score,
                          levels = c("PS","DRS","None")),
           ## Adjustment method
           adjust = SeqGsub(method,
                            c(".*Match.*", "Matching",
                              ".*Strata.*", "Stratification",
                              ".*SIptw.*", "IPTW",
                              ".*Mw.*", "MW",
                              "^unadj.*", "Unadjusted")),
           adjust = factor(adjust,
                           levels = c("Matching","Stratification","IPTW","MW","Unadjusted")),
           ## Approximation method
           approx = SeqGsub(method,
                            c(".*E|.*E_true", "Efron",
                              ".*B|.*B_true", "Breslow")),
           approx = if_else(approx %in% c("Efron","Breslow"),
                            approx,
                            ""),
           ## When data sharing method is "truth", these are results from AnalyzeSiteTruth()
           data = factor(data,
                         levels = c("meta", "summary", "risksets", "dataset", "truth"),
                         labels = c("meta", "summary", "risk set", "individual", "truth")),
           true = (data %in% "truth")
           ) %>%
    ## Remove Efron up front
    filter(approx != "Efron")


cat("
###  Show first part of combined data frame\n")
head(dfAnalysis)
cat("
###  Show last part of combined data frame\n")
tail(dfAnalysis)


cat("
###
### Scenario name assignment
################################################################################\n")

scenario_names <- c("base",
                    "tx25",
                    "tx10",
                    "tx05",
                    "dis1perc",
                    "dis0.1perc",
                    "dis0.01perc",
                    "dis_var_inc",
                    "harm_hr1.2",
                    "harm_hr2.0",
                    "2_5k_sites",
                    "4_5k_sites",
                    "8_5k_sites",
                    "vary_conf_count",
                    "noise_cov_hr2.0",
                    "instruments_hr2.0",
                    "rct_hr2.0",
                    "new_tx25",
                    "new_tx10",
                    "new_tx05",
                    "16_5k_sites",
                    "new_harm_hr1.2_tx25",
                    "new_harm_hr1.2_tx10",
                    "new_harm_hr1.2_tx05",
                    "new_prot_hr0.8",
                    "new_prot_hr0.8_tx25",
                    "new_prot_hr0.8_tx10",
                    "new_prot_hr0.8_tx05",
                    "new_dis1perc_hr1.2",
                    "new_dis0.1perc_hr1.2",
                    "new_dis0.01perc_hr1.2",
                    "new_dis1perc_hr0.8",
                    "new_dis0.1perc_hr0.8",
                    "new_dis0.01perc_hr0.8",
                    "2_10k_sites",
                    "4_10k_sites",
                    "8_10k_sites",
                    "16_10k_sites",
                    "2_20k_sites",
                    "4_20k_sites",
                    "8_20k_sites",
                    "16_20k_sites",
                    "2_1k_sites",
                    "4_1k_sites",
                    "8_1k_sites",
                    "16_1k_sites",
                    "2_500_sites",
                    "4_500_sites",
                    "8_500_sites",
                    "16_500_sites",
                    "covariates_small_hr0.8",
                    "covariates_small_hr1.0",
                    "covariates_small_hr1.2",
                    "covariates_moderate_hr0.8",
                    "covariates_moderate_hr1.0",
                    "covariates_moderate_hr1.2",
                    "covariates_large_hr0.8",
                    "covariates_large_hr1.0",
                    "covariates_large_hr1.2")

scenario_descriptions <- c("base 50% treatment 5% incidence 4 sites",
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

df_scenario_names_desc <- data_frame(scenario = seq_along(scenario_names),
                                     scenario_name = scenario_names,
                                     scenario_desc = scenario_descriptions)

dfAnalysis <- left_join(dfAnalysis,
                        df_scenario_names_desc) %>%
    as_data_frame
dfAnalysis %>%
    select(scenario, outcome, score, adjust, approx, data, method, true)


cat("
###
### Extreme estimate handling
################################################################################\n")

cat("
###  Show extreme coefficients (coef > 5; will be dropped)\n")
dfAnalysis %>%
    filter(!is.na(coef), abs(coef) > 5) %>%
    select(scenario, part, iter, outcome, method, data, coef, var, scenario_name) %>%
    print(n = Inf)

cat("
###  Show extreme variance (var > 5; will be dropped)\n")
dfAnalysis %>%
    filter(!is.na(var), abs(var) > 5) %>%
    select(scenario, part, iter, outcome, method, data, coef, var, scenario_name) %>%
    print(n = Inf)

cat("
###  Replace invalid estimates with NA\n")
dfAnalysis <- dfAnalysis %>%
    mutate(coef = if_else(!is.na(coef) & abs(coef) > 5, as.numeric(NA), coef),
           var = if_else(!is.na(var) & abs(var) > 5, as.numeric(NA), var),
           ## Only allow observations that have valid coef AND var estimates
           coef = if_else(is.na(var), as.numeric(NA), coef),
           var = if_else(is.na(coef), as.numeric(NA), var))


cat("
###
### Assess data
################################################################################\n")

cat("###  True effect calculation\n")
dfTrue <- dfAnalysis %>%
    ## Filter true effects only
    filter(true) %>%
    group_by(scenario, scenario_name, scenario_desc,
             outcome, score, adjust, approx, data, method) %>%
    summarize(true_coef = mean(coef, na.rm = TRUE))
## Only Breslow remains here because truth was calculated for IPD only.
stopifnot(all(dfTrue$approx == "Breslow"))
dfTrue


cat("###  Construct operational characteristics\n")
dfSummary <- dfAnalysis %>%
    ## Drop true effects
    filter(!true) %>%
    ## Collapsing over part & iter (all interations)
    group_by(scenario, scenario_name, scenario_desc,
             outcome, score, adjust, approx, data, method) %>%
    summarize(mean_coef = mean(coef, na.rm = TRUE),
              mean_var  = mean(var, na.rm = TRUE),
              true_var  = var(coef, na.rm = TRUE),
              ## ratio of means
              ratio_mean_true_var = mean_var / true_var,
              ratio_mean_true_se  = sqrt(mean_var) / sqrt(true_var),
              ## mean of ratios
              ## This is mean(vector / scalar) structure
              mean_ratio_est_true_se = mean((sqrt(var) / sqrt(true_var)), na.rm = TRUE),
              ## Examine iteration failure
              na_coef   = mean(is.na(coef)),
              na_var    = mean(is.na(var))) %>%
    ## Left join true coefficient values
    left_join(.,
              dfTrue %>% ungroup() %>% select(-scenario_name, -scenario_desc, -approx, -data, -method),
              by = c("scenario", "outcome", "score", "adjust")) %>%
    ## Calculate additional measures that require truth
    mutate(bias_coef = mean_coef - true_coef,
           mse_coef = true_var + bias_coef^2)
## Peek results.
dfSummary %>%
    select(scenario, outcome, score, adjust, approx, data, method, mean_coef, true_coef)

cat("###  Add true_var to iteration level dataset\n")
dfAnalysis <- left_join(x = dfAnalysis,
                        y = select(dfSummary,
                                   scenario, data, outcome, score, adjust, approx, method,
                                   true_var))


## SE_hat / SE ratio for each iteration
dfAnalysis <- dfAnalysis %>%
    group_by(scenario, data, outcome, score, adjust, approx, method) %>%
    mutate(ratio_var = var / true_var,
           ratio_se = sqrt(var) / sqrt(true_var))


cat("
###  Examine SE ratios for table and graph\n")

dfSummary %>% filter(outcome == "survival") %>%
    select(mean_ratio_est_true_se) %>%
    mutate(mean_ratio_est_true_se = round(mean_ratio_est_true_se, 5))

dfAnalysis %>% filter(outcome == "survival") %>%
    group_by(scenario, data, outcome, score, adjust, approx, method) %>%
    summarize(mean_se_ratio = round(mean(ratio_se), 5))


cat("
###  Assess weighting methods\n")
## Check if meta-analysis IPTW are bad
filter(dfSummary, adjust %in% c("IPTW","MW") & data == "meta")[,c("scenario","outcome","method","mean_coef","mean_var","true_var")] %>% as.data.frame

filter(dfAnalysis, adjust %in% c("IPTW","MW") & iter %in% c(1,5) & outcome == "survival" & data == "meta" & scenario == 5)[,c("scenario","iter","outcome","method","coef","var","coef_site1", "coef_site2", "coef_site3", "coef_site4",
"var_site1", "var_site2", "var_site3", "var_site4")] %>% as.data.frame


cat("
###  Assess scenarios with coef very close to zero in weighted meta-analysis\n")
dfAnalysis5 <- filter(dfAnalysis, scenario == 5 &
                                  data == "meta" &
                                  adjust %in% c("IPTW","MW") &
                                  approx != "Efron" &
                                  outcome == "binary" &
                                  !grepl("truth", as.character(data)))
dfAnalysis5[abs(dfAnalysis5$coef) < 0.0001,]


cat("
###  Try no event in either group data in logistic regression\n")
data1 <- data.frame(outcome = rep(0,10), exposure = rep(c(0,1),5))
glm1 <- glm(formula = outcome ~ exposure,
        family  = binomial(link = "logit"),
        data    = data1)
summary(glm1)


cat("
###  Write out summary\n")

## Write out full summary
openxlsx::write.xlsx(as.data.frame(dfSummary),
                     file = "./summary/summary.xlsx")


cat("
###  Assess NA fractions\n")
filter(dfSummary, na_coef > 0) %>%
    select(scenario, data, score, adjust, na_coef)

table(filter(dfSummary, na_coef > 0)$scenario)


cat("
###  Least biased methods\n")
filter(dfSummary, outcome == "binary") %>%
    arrange(scenario, abs(bias_coef))
filter(dfSummary, outcome == "survival") %>%
    arrange(scenario, abs(bias_coef))


cat("
###  Smallest variance methods\n")
filter(dfSummary, outcome == "binary") %>%
    arrange(scenario, abs(true_var))
filter(dfSummary, outcome == "survival") %>%
    arrange(scenario, abs(true_var))


cat("
###  Smallest MSE methods\n")
filter(dfSummary, outcome == "binary") %>%
    arrange(scenario, abs(mse_coef))
filter(dfSummary, outcome == "survival") %>%
    arrange(scenario, abs(mse_coef))

cat("
###  Equivalent methods (up to 7 decimals)\n")

## Sort
dfSummary2 <- dfSummary[,c("mean_coef","mean_var","data","outcome","score","adjust","approx")] %>%
    arrange(outcome, mean_coef, mean_var)

## Only duplicated rows are wanted
dupFromAbove <- duplicated(round(dfSummary2[c("mean_coef","mean_var")], 7))
dupFromBelow <- duplicated(round(dfSummary2[c("mean_coef","mean_var")], 7), fromLast = TRUE)
dfSummary2[(dupFromAbove + dupFromBelow) > 0,] %>%
    as.data.frame


cat("
###  Average coefficient and average SE ratio (Survival)\n")

dfSummary %>%
    filter(approx != "Efron", outcome == "survival") %>%
    select(mean_coef, mean_ratio_est_true_se) %>%
    as.data.frame


cat("
###
### Save as data file
################################################################################\n")

save(dfAnalysis, dfSummary,
     scenario_names, scenario_descriptions,
     file = paste0(dirName,"/analysis_summary_data.RData"))


################################################################################
cat("
###
### Record package versions etc
################################################################################\n")
print(sessionInfo())
## Record execution time
end_time <- Sys.time()
cat("\n### Started  ", as.character(start_time), "\n")
cat("### Finished ", as.character(end_time), "\n")
print(end_time - start_time)
## Stop sinking to a file if active
if (sink.number() != 0) {sink()}
