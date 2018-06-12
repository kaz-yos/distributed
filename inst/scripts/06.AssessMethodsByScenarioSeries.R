#!/usr/local/bin/Rscript

################################################################################
### Assess result files by scenario series
##
## Created on: 2017-10-18
## Author: Kazuki Yoshida
################################################################################


###
### Prepare environment
################################################################################

## sink() if being run non-interactively
if (sink.number() != 0) {sink()}
..scriptFileName.. <- gsub("^--file=", "", Filter(function(x) {grepl("^--file=", x)}, commandArgs()))
if (length(..scriptFileName..) == 1) {
    sink(file = paste0("./log/", basename(..scriptFileName..), ".txt"), split = TRUE)
    options(width = 100)
    options(max.print = 3000)
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
nCores <- min(parallel::detectCores(), 8)
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
library(openxlsx)

dirName <- "./data/"
load(paste0(dirName,"/analysis_summary_data.RData"))


cat("
###
### Graphing by scenario series
################################################################################\n")

## Figure dimensions
width <- 8
height <- 5.6

## Color configuration
## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999")

theme_clean <- theme_bw() + theme(legend.key = element_blank(),
                                  axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5),
                                  strip.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5))

my_collapse_labels_lines <- function (labels)  {
    out <- do.call("Map", c(list(paste, sep = " "), labels))
    list(unname(unlist(out)))
}

df_summary_scenario_series <- dfSummary %>%
    filter(
        ## Drop unadjusted results
        !(adjust %in% "Unadjusted"),
        ## Restrict to survival outcome (do not mix with binary ones)
        outcome == "survival",
        ## Do not use Efron where applicable
        approx != "Efron")

## HR data_frame for horizontal lines
df_hr <- data_frame(hr = c("HR 0.8","HR 1.0","HR 1.2"),
                    hr_value = c(0.8, 1.0, 1.2)) %>%
    mutate(hr = factor(hr, levels = c("HR 1.2", "HR 1.0", "HR 0.8")))

cat("
###  Construct datasets\n")
###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence <-
    df_summary_scenario_series %>%
    filter(grepl("^base|^new_tx|harm_hr1.2|^new_prot_hr0.8", scenario_name)) %>%
    mutate(prevalence = case_when(grepl("tx25", scenario_name) ~ 25,
                                  grepl("tx10", scenario_name) ~ 10,
                                  grepl("tx05", scenario_name) ~ 5,
                                  TRUE ~ 50),
           hr = case_when(grepl("hr1.2", scenario_name) ~ "HR 1.2",
                          grepl("hr0.8", scenario_name) ~ "HR 0.8",
                          TRUE ~ "HR 1.0"),
           hr = factor(hr, levels = c("HR 1.2", "HR 1.0", "HR 0.8")))
df_summary_scenario_series_vary_exposure_prevalence %>%
    group_by(scenario, scenario_name, scenario_desc) %>%
    summarize(n = n())

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence <-
    df_summary_scenario_series %>%
    filter(grepl("base|dis.*perc|^harm_hr1.2$|^new_prot_hr0.8$", scenario_name)) %>%
    mutate(incidence = case_when(grepl("dis1perc", scenario_name) ~ 1,
                                 grepl("dis0.1perc", scenario_name) ~ 0.1,
                                 grepl("dis0.01perc", scenario_name) ~ 0.01,
                                 TRUE ~ 5),
           hr = case_when(grepl("hr1.2", scenario_name) ~ "HR 1.2",
                          grepl("hr0.8", scenario_name) ~ "HR 0.8",
                          TRUE ~ "HR 1.0"),
           hr = factor(hr, levels = c("HR 1.2", "HR 1.0", "HR 0.8")))
df_summary_scenario_series_vary_disease_incidence %>%
    group_by(scenario, scenario_name, scenario_desc) %>%
    summarize(n = n())

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count <-
    df_summary_scenario_series %>%
    filter(grepl("k_sites", scenario_name) | grepl("500_sites", scenario_name)) %>%
    mutate(site_count = case_when(grepl("^2", scenario_name) ~ "2 sites",
                                  grepl("^4", scenario_name) ~ "4 sites",
                                  grepl("^8", scenario_name) ~ "8 sites",
                                  grepl("^16", scenario_name) ~ "16 sites"),
           site_count = factor(site_count,
                               levels = c("2 sites", "4 sites", "8 sites", "16 sites")),
           site_size = case_when(grepl("500", scenario_name) ~ 500,
                                 grepl("1k", scenario_name) ~ 1000,
                                 grepl("5k", scenario_name) ~ 5000,
                                 grepl("10k", scenario_name) ~ 10000,
                                 grepl("20k", scenario_name) ~ 20000))
df_summary_scenario_series_vary_site_count %>%
    group_by(scenario, scenario_name, scenario_desc) %>%
    summarize(n = n())

###   Vary covariate distribution All HR
df_summary_scenario_series_vary_covariate_distribution <-
    df_summary_scenario_series %>%
    filter(grepl("^base$|^covariates_|^harm_hr1.2$|^new_prot_hr0.8$", scenario_name)) %>%
    mutate(covariates = case_when(grepl("covariates_small", scenario_name) ~ "Small",
                                  grepl("covariates_moderate", scenario_name) ~ "Moderate",
                                  grepl("covariates_large", scenario_name) ~ "Large",
                                  TRUE ~ "None"),
           covariates = factor(covariates,
                               ## Reverse to be consistent with other series
                               levels = rev(c("None","Small","Moderate","Large"))),
           hr = case_when(grepl("hr1.2", scenario_name) ~ "HR 1.2",
                          grepl("hr0.8", scenario_name) ~ "HR 0.8",
                          TRUE ~ "HR 1.0"),
           hr = factor(hr, levels = c("HR 1.2", "HR 1.0", "HR 0.8")))
df_summary_scenario_series_vary_covariate_distribution %>%
    group_by(scenario, scenario_name, scenario_desc) %>%
    summarize(n = n())


cat("
###  Average point estimates\n")

pdf(file = "./summary/series_summary_mean_coef.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(exp_mean_coef = exp(mean_coef),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = exp_mean_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_point(mapping = aes(y = exp(true_coef), shape = NULL), color = "red") +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 0.20 - 1/70 * as.numeric(interaction(data, hr)))) +
    geom_hline(data = df_hr, mapping = aes(yintercept = hr_value), size = 0.2) +
    ## scale_y_continuous(breaks = seq(from = 0.7, to = 1.3, by = 0.05)) +
    scale_y_continuous(limit = c(0.7, 1.4)) +
    facet_grid(hr ~ score + adjust, scale = "free_y") +
    labs(title = NULL,
         x = "Exposure Prevalence (%)",
         y = "Mean HR",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(exp_mean_coef = exp(mean_coef),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(incidence), y = exp_mean_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_point(mapping = aes(y = exp(true_coef), shape = NULL), color = "red") +
    geom_label(mapping = aes(label = fail_prop, y = 3.0 - 0.2 * as.numeric(data)), size = 2) +
    geom_hline(data = df_hr, mapping = aes(yintercept = hr_value), size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Outcome Incidence (%)",
         y = "Mean HR",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(exp_mean_coef = exp(mean_coef),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(site_size), y = exp_mean_coef, group = interaction(data, site_count),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_point(mapping = aes(y = exp(true_coef), shape = NULL), color = "red") +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 1.5 - 1/30 * as.numeric(data)), size = 2) +
    geom_hline(yintercept = 1.0, size = 0.2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         x = "Site size",
         y = "Mean Hazard Ratio",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distribution All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(exp_mean_coef = exp(mean_coef),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(covariates), y = exp_mean_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_point(mapping = aes(y = exp(true_coef), shape = NULL), color = "red") +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 0.20 - 1/70 * as.numeric(interaction(data, hr)))) +
    geom_hline(data = df_hr, mapping = aes(yintercept = hr_value), size = 0.2) +
    ## scale_y_continuous(breaks = seq(from = 0.7, to = 1.3, by = 0.05)) +
    scale_y_continuous(limit = c(0.7, 1.4)) +
    facet_grid(hr ~ score + adjust, scale = "free_y") +
    labs(title = NULL,
         x = "Covariate variability across sites",
         y = "Mean HR",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()


cat("
###  Bias\n")

pdf(file = "./summary/series_summary_bias.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = bias_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 0.20 - 1/70 * as.numeric(interaction(data, hr)))) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust, scale = "free_y") +
    labs(title = NULL,
         x = "Exposure Prevalence (%)",
         y = "Bias on log HR scale",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(incidence), y = bias_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_label(mapping = aes(label = fail_prop, y = 1.5 - 0.2 * as.numeric(data)), size = 2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Outcome Incidence (%)",
         y = "Bias on log HR scale",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(site_size), y = bias_coef, group = interaction(data, site_count),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 1.5 - 1/7 * as.numeric(data)), size = 2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         x = "Site size",
         y = "Bias on log HR scale",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distributions All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(covariates), y = bias_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    ## Commented out after confirming no issues.
    ## geom_label(mapping = aes(label = fail_prop, y = 0.20 - 1/70 * as.numeric(interaction(data, hr)))) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust, scale = "free_y") +
    labs(title = NULL,
         x = "Covariate variability across sites",
         y = "Bias on log HR scale",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()


cat("
###  Simulation-based SD of point estimates\n")

pdf(file = "./summary/series_summary_se_coef.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(se_coef = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Simulation SD of log HR estimates",
         x = "Exposure Prevalence (%)",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(se_coef = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(incidence), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Simulation SD of log HR estimates",
         x = "Outcome Incidence (%)",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(se_coef = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(site_size), y = se_coef, group = interaction(data, site_count),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         y = "Simulation SD of log HR estimates",
         x = "Site size",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distributions All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(se_coef = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(covariates), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Simulation SD of log HR estimates",
         x = "Covariate variability across sites",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()


cat("
###  Mean estimated SE\n")

pdf(file = "./summary/series_summary_se_est.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(se_coef = sqrt(mean_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Mean SE of log HR estimates",
         x = "Exposure Prevalence (%)",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(se_coef = sqrt(mean_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(incidence), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Mean SE of log HR estimates",
         x = "Outcome Incidence (%)",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(se_coef = sqrt(mean_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(site_size), y = se_coef, group = interaction(data, site_count),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         y = "Mean SE of log HR estimates",
         x = "Site size",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distributions All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(se_coef = sqrt(mean_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(covariates), y = se_coef, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 0, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         y = "Mean SE of log HR estimates",
         x = "Covariate variability across sites",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()


cat("
###  sqrt(mean(est_var)) / Simulation-based SD of point estimates\n")

pdf(file = "./summary/series_summary_ratio_se.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(se_ratio = sqrt(mean_var) / sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = se_ratio, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 1, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Exposure Prevalence (%)",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(se_ratio = sqrt(mean_var) / sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(incidence), y = se_ratio, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 1, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Outcome Incidence (%)",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(se_ratio = sqrt(mean_var) / sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(site_size), y = se_ratio, group = interaction(data, site_count),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 1, size = 0.2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         x = "Site size",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distributions All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(se_ratio = sqrt(mean_var) / sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    ggplot(mapping = aes(x = factor(covariates), y = se_ratio, group = interaction(data, hr),
                         linetype = data, shape = data)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    geom_hline(yintercept = 1, size = 0.2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Covariate variability across sites",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()


cat("
###  sqrt(mean(est_var)) vs Simulation-based SD of point estimates (overplot)\n")

pdf(file = "./summary/series_summary_overplot_se.pdf", width = width, height = height, family = "sans")

###   Vary exposure prevalence (better confounding setting) All HR
df_summary_scenario_series_vary_exposure_prevalence %>%
    mutate(mean_se = sqrt(mean_var),
           true_se = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    gather(key = key, value = value, mean_se, true_se) %>%
    ggplot(mapping = aes(x = factor(prevalence), y = value, group = interaction(data, hr, key),
                         linetype = data, shape = data, color = key)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Exposure Prevalence (%)",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary disease incidence All HR
df_summary_scenario_series_vary_disease_incidence %>%
    mutate(mean_se = sqrt(mean_var),
           true_se = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    gather(key = key, value = value, mean_se, true_se) %>%
    ggplot(mapping = aes(x = factor(incidence), y = value, group = interaction(data, hr, key),
                         linetype = data, shape = data, color = key)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Outcome Incidence (%)",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary site count All site sizes
df_summary_scenario_series_vary_site_count %>%
    mutate(mean_se = sqrt(mean_var),
           true_se = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    gather(key = key, value = value, mean_se, true_se) %>%
    ggplot(mapping = aes(x = factor(site_size), y = value, group = interaction(data, key),
                         linetype = data, shape = data, color = key)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    facet_grid(site_count ~ score + adjust) +
    labs(title = NULL,
         x = "Site size",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

###   Vary covariate distributions All HR
df_summary_scenario_series_vary_covariate_distribution %>%
    mutate(mean_se = sqrt(mean_var),
           true_se = sqrt(true_var),
           fail_prop = round(na_coef * 100, 1)) %>%
    gather(key = key, value = value, mean_se, true_se) %>%
    ggplot(mapping = aes(x = factor(covariates), y = value, group = interaction(data, hr, key),
                         linetype = data, shape = data, color = key)) +
    geom_line() +
    geom_point(size = 3, alpha = 1/2) +
    facet_grid(hr ~ score + adjust) +
    labs(title = NULL,
         x = "Covariate variability across sites",
         y = "Average SE / Simulation SD",
         linetype = NULL, shape = NULL) +
    theme_clean

dev.off()

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
