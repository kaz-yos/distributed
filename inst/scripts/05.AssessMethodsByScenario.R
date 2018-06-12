#!/usr/local/bin/Rscript

################################################################################
### Assess result files (simplified for December 2016 reporting)
##
## Created on: 2016-12-07
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
### Create a table summarizing the base case
################################################################################\n")

base_table <- dfSummary %>%
    filter(approx != "Efron",
           !grepl("truth", as.character(data)),
           outcome == "survival",
           scenario == 1) %>%
    group_by(score, adjust, data) %>%
    rename(Score = score,
           Adjustment = adjust,
           Data = data) %>%
    transmute(`Mean log HR`      = mean_coef,
              `Simulation SE`    = sqrt(true_var),
              `Mean SE estimate` = sqrt(mean_var),
              `SE ratio`         = `Mean SE estimate` / `Simulation SE`) %>%
    mutate(`Mean log HR`      = sprintf("%.3f", `Mean log HR`),
           `Simulation SE`    = sprintf("%.3f", `Simulation SE`),
           `Mean SE estimate` = sprintf("%.3f", `Mean SE estimate`),
           `SE ratio`         = sprintf("%.3f", `SE ratio`)) %>%
    print(n = Inf)


openxlsx::write.xlsx(base_table,
                     file = "./summary/base_scenario.xlsx")


cat("
###
### Graphing
################################################################################\n")

## Color configuration
## http://www.cookbook-r.com/Graphs/Colors_(ggplot2)/#a-colorblind-friendly-palette
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#CC79A7", "#F0E442", "#0072B2", "#D55E00", "#999999")

cat("
###  Drop unadjusted and truth\n")
dfAnalysis <- dfAnalysis %>%
    filter(!(adjust %in% c("Truth","Unadjusted")))

theme_clean <- theme_bw() + theme(legend.key = element_blank(),
                                  axis.text.x = element_text(angle = 90, hjust = 0.5, vjust = 0.5),
                                  strip.background = element_blank(),
                                  plot.title = element_text(hjust = 0.5))

my_collapse_labels_lines <- function (labels)  {
    out <- do.call("Map", c(list(paste, sep = " "), labels))
    list(unname(unlist(out)))
}


cat("
###  All methods by non-NA estimates (survival)\n")
pdf(file = "./summary/summary_nonNA_survival.pdf", width = 7.5, height = 6, family = "sans")

dfAnalysis %>%
    filter(approx != "Efron",
           outcome == "survival") %>%
    group_by(score, adjust, data, scenario) %>%
    summarize(perc_success = sum(!is.na(coef)) / length(coef) * 100) %>%
    group_by(scenario) %>%
    nest() %>%
    left_join(data_frame(scenario = seq_along(scenario_names),
                         scenario_name = scenario_names,
                         scenario_description = scenario_descriptions)) %>%
    mutate(gg = map(scenario, function(i) {
        ggplot(data = data[[i]],
               mapping = aes(x    = data,
                             y    = perc_success,
                             fill = data)) +
            geom_bar(stat = "identity") +
            facet_wrap(~ score + adjust, scale = "fixed", nrow = 2,
                       labeller = my_collapse_labels_lines) +
            scale_y_continuous(limit = c(0,100)) +
            scale_fill_manual(values = cbPalette) +
            labs(title = paste0("Survival analysis successful iterations (%).\n Scenario ",
                                i,
                                ": ",
                                .$scenario_description[i]),
                 y = "% successful iterations") +
            guides(fill = FALSE) +
            theme_clean

    })) %>%
    magrittr::extract2("gg")

dev.off()


cat("
###  All methods by coefficient points (survival)\n")
pdf(file = "./summary/summary_coefs_survival.pdf", width = 7.5, height = 4, family = "sans")

dfAnalysis %>%
    filter(approx != "Efron",
           outcome == "survival") %>%
    group_by(scenario) %>%
    nest() %>%
    left_join(data_frame(scenario = seq_along(scenario_names),
                         scenario_name = scenario_names,
                         scenario_description = scenario_descriptions)) %>%
    mutate(gg = map(scenario, function(i) {
        max_exp_coef <- max(exp(data[[i]]$coef), na.rm = TRUE)
        if (max_exp_coef < 10) {
            breaks <- c(1/20, 1/10, 1/5, 1/2, 0.8, 0.9, 1, 1.1, 1.2, 2, 5, 10, 20)
        } else {
            breaks <- c(1/20, 1/10, 1/5, 1/2, 0.8, 1, 1.2, 2, 5, 10, 20)
        }
        truth_data <- data[[i]] %>%
            filter(data == "truth") %>%
            group_by(score, adjust) %>%
            summarize(coef_true = mean(coef, na.rm = TRUE))
        ggplot(data = data[[i]] %>% filter(data != "truth"),
               mapping = aes(x = data,
                             y = exp(coef))) +
            geom_hline(data = truth_data,
                       mapping = aes(yintercept = exp(coef_true))) +
            geom_boxplot() +
            ## stat_summary(fun.y = mean, colour = "darkred", geom = "point",
            ##              shape = 4, size = 1, show.legend = FALSE) +
            scale_y_log10(breaks = breaks) +
            facet_grid(. ~ score + adjust, scale = "fixed") +
            labs(title = paste0("Survival analysis log HR.\n Scenario ",
                                i,
                                ": ",
                                .$scenario_description[i]),
                 x = NULL, y = "HR") +
            guides(color = FALSE) +
            theme_clean
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./summary/summary_coefs_survival_base.pdf", width = 6, height = 3, family = "sans")

dfAnalysis %>%
    filter(approx != "Efron",
           outcome == "survival",
           scenario == 1,
           adjust != "Unadjusted") %>%
    group_by(scenario) %>%
    nest() %>%
    left_join(data_frame(scenario = seq_along(scenario_names),
                         scenario_name = scenario_names,
                         scenario_description = scenario_descriptions)) %>%
    mutate(gg = map(scenario, function(i) {
        ## breaks <- seq(from = 0.5, to = 1.5, by = 0.05)
        truth_data <- data[[i]] %>%
            filter(data == "truth") %>%
            group_by(score, adjust) %>%
            summarize(coef_true = mean(coef, na.rm = TRUE))
        ggplot(data = data[[i]] %>% filter(data != "truth"),
               mapping = aes(x = data,
                             y = coef)) +
            geom_hline(data = truth_data,
                       mapping = aes(yintercept = coef_true)) +
            geom_boxplot() +
            ## stat_summary(fun.y = mean, colour = "darkred", geom = "point",
            ##              shape = 4, size = 1, show.legend = FALSE) +
            facet_grid(. ~ score + adjust, scale = "fixed") +
            labs(title = NULL, x = NULL, y = "log HR") +
            guides(color = FALSE) +
            theme_clean
    })) %>%
    magrittr::extract2("gg")

dev.off()



cat("
###  All methods by SE/true SE ratio points (survival)\n")
pdf(file = "./summary/summary_se_survival.pdf", width = 7.5, height = 4, family = "sans")

dfAnalysis %>%
    filter(approx != "Efron",
           !grepl("truth", as.character(data)),
           outcome == "survival") %>%
    summarize(min_ratio_se = min(ratio_se, na.rm = TRUE),
              max_ratio_se = max(ratio_se, na.rm = TRUE),
              mean_ratio_se = mean(ratio_se, na.rm = TRUE))

## Use Breslow or non-approximation only
dfAnalysis %>%
    filter(approx != "Efron",
           !grepl("truth", as.character(data)),
           outcome == "survival") %>%
    group_by(scenario) %>%
    nest() %>%
    left_join(data_frame(scenario = seq_along(scenario_names),
                         scenario_name = scenario_names,
                         scenario_description = scenario_descriptions)) %>%
    mutate(gg = map(scenario, function(i) {
        ggplot(data = data[[i]],
               mapping = aes(x = data,
                             y = ratio_se)) +
            geom_hline(yintercept = 1, size = 0.2) +
            geom_boxplot() +
            ## stat_summary(fun.y = mean, colour = "darkred", geom = "point",
            ##              shape = 4, size = 1, show.legend = FALSE) +
            facet_grid(. ~ score + adjust, scale = "fixed") +
            scale_y_log10(limit = c(NA, NA), breaks = c(1/5,1/2,0.75,0.9,1,1.1,1.5,2,5)) +
            labs(title = paste0("Survival analysis SE estimate accuracy.\n Scenario ",
                                i,
                                ": ",
                                .$scenario_description[i]),
                 x = NULL, y = "Estimated SE(log HR) / simulation SE(log HR)") +
            guides(color = FALSE) +
            theme_clean
    })) %>%
    magrittr::extract2("gg")

dev.off()


pdf(file = "./summary/summary_se_survival_base.pdf", width = 6, height = 3, family = "sans")

## Use Breslow or non-approximation only
dfAnalysis %>%
    filter(approx != "Efron",
           !grepl("truth", as.character(data)),
           outcome == "survival",
           scenario == 1,
           adjust != "Unadjusted") %>%
    group_by(scenario) %>%
    nest() %>%
    left_join(data_frame(scenario = seq_along(scenario_names),
                         scenario_name = scenario_names,
                         scenario_description = scenario_descriptions)) %>%
    mutate(gg = map(scenario, function(i) {
        ggplot(data = data[[i]],
               mapping = aes(x = data,
                             y = ratio_se)) +
            geom_hline(yintercept = 1, size = 0.2) +
            geom_boxplot() +
            ## stat_summary(fun.y = mean, colour = "darkred", geom = "point",
            ##              shape = 4, size = 1, show.legend = FALSE) +
            facet_grid(. ~ score + adjust, scale = "fixed") +
            scale_y_log10(limit = c(NA, NA), breaks = seq(from = 0.9, to = 1.3, by = 0.05)) +
            labs(title = NULL, x = NULL, y = "Ratio of Estimated to Simulation SD") +
            guides(color = FALSE) +
            theme_clean
    })) %>%
    magrittr::extract2("gg")

dev.off()


cat("
###  Individual site analysis (survival)\n")

###   Success rate
pdf(file = "./summary/summary_individual_sites_success_survival.pdf", width = 10, height = 7, family = "sans")

dfAnalysis %>%
    filter(outcome == "survival",
           score %in% c("PS","DRS"),
           approx %in% c("","Breslow"),
           data == "meta") %>%
    select(starts_with("coef_site"),
           scenario_name,
           scenario_desc) %>%
    gather(key = site, value = coef_site, starts_with("coef_site")) %>%
    mutate(site = as.factor(as.numeric(gsub("coef_site","", site)))) %>%
    group_by(scenario) %>%
    nest() %>%
    mutate(gg = map(seq_along(data), function(i) {
        ggplot(data = data[[i]], mapping = aes(x = site,
                                               y = mean(!is.na(coef_site)),
                                               color = site)) +
            geom_bar(stat = "identity") +
            facet_grid(. ~ score + adjust) +
            labs(title = sprintf("Site-Specific Success Rate.\nScenario %d: %s",
                                 i,
                                 data[[i]]$scenario_desc[1])) +
            theme_bw() + theme(legend.key = element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1),
                               plot.title = element_text(hjust = 0.5))
    })) %>% magrittr::extract2("gg")

dev.off()


###   Coefficients
pdf(file = "./summary/summary_individual_sites_coef_survival.pdf", width = 10, height = 7, family = "sans")

dfAnalysis %>%
    filter(outcome == "survival",
           score %in% c("PS","DRS"),
           approx %in% c("","Breslow"),
           data == "meta") %>%
    select(starts_with("coef_site"),
           scenario_name,
           scenario_desc) %>%
    gather(key = site, value = coef_site, starts_with("coef_site")) %>%
    mutate(site = as.factor(as.numeric(gsub("coef_site","", site)))) %>%
    group_by(scenario) %>%
    nest() %>%
    mutate(gg = map(seq_along(data), function(i) {
        ggplot(data = data[[i]], mapping = aes(x = site,
                                               y = coef_site,
                                               color = site)) +
            geom_boxplot() +
            facet_grid(. ~ score + adjust) +
            labs(title = sprintf("Site-Specific Coefficients.\nScenario %d: %s",
                                 i,
                                 data[[i]]$scenario_desc[1])) +
            theme_bw() + theme(legend.key = element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1),
                               plot.title = element_text(hjust = 0.5))
    })) %>% magrittr::extract2("gg")

dev.off()


###   Simulation SE
pdf(file = "./summary/summary_individual_sites_se_survival.pdf", width = 10, height = 7, family = "sans")

dfAnalysis %>%
    filter(outcome == "survival",
           score %in% c("PS","DRS"),
           approx %in% c("","Breslow"),
           data == "meta") %>%
    select(starts_with("coef_site"),
           scenario_name,
           scenario_desc) %>%
    gather(key = site, value = coef_site, starts_with("coef_site")) %>%
    mutate(site = as.factor(as.numeric(gsub("coef_site","", site)))) %>%
    group_by(site, score, adjust, scenario) %>%
    mutate(var_site = var(coef_site, na.rm = TRUE)) %>%
    group_by(scenario) %>%
    nest() %>%
    mutate(gg = map(seq_along(data), function(i) {
        ggplot(data = data[[i]], mapping = aes(x = site,
                                               y = var_site,
                                               color = site)) +
            geom_point() +
            facet_grid(. ~ score + adjust) +
            labs(title = sprintf("Site-Specific Sim Variance.\nScenario %d: %s",
                                 i,
                                 data[[i]]$scenario_desc[1])) +
            theme_bw() + theme(legend.key = element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1),
                               plot.title = element_text(hjust = 0.5))
    })) %>% magrittr::extract2("gg")

dev.off()


###   Estimated SE

pdf(file = "./summary/summary_individual_sites_est_se_survival.pdf", width = 10, height = 7, family = "sans")

dfAnalysis %>%
    filter(outcome == "survival",
           score %in% c("PS","DRS"),
           approx %in% c("","Breslow"),
           data == "meta") %>%
    select(starts_with("var_site"),
           scenario_name,
           scenario_desc) %>%
    gather(key = site, value = var_site, starts_with("var_site")) %>%
    mutate(site = as.factor(as.numeric(gsub("var_site","", site)))) %>%
    group_by(site, score, adjust, scenario) %>%
    group_by(scenario) %>%
    nest() %>%
    mutate(gg = map(seq_along(data), function(i) {
        ggplot(data = data[[i]], mapping = aes(x = site,
                                               y = var_site,
                                               color = site)) +
            geom_boxplot() +
            facet_grid(. ~ score + adjust) +
            labs(title = sprintf("Site-Specific Sim Variance.\nScenario %d: %s",
                                 i,
                                 data[[i]]$scenario_desc[1])) +
            theme_bw() + theme(legend.key = element_blank(),
                               axis.text.x = element_text(angle = 90, hjust = 1),
                               plot.title = element_text(hjust = 0.5))
    })) %>% magrittr::extract2("gg")

dev.off()


cat("
###  Weighting and meta-analysis assessment\n")

pdf(file = "./summary/summary_weighted_metaanalysis_assessment.pdf", width = 10, height = 7, family = "sans")

dfAnalysis$site_NA <- factor(rowSums(is.na(dfAnalysis[paste0("coef_site", 1:4)])))

split(dfAnalysis, dfAnalysis$scenario) %>%
    ## Loop over scenarios
    sapply(., function(df) {

        scenarioCount <- unique(df$scenario)[1]

        df <- filter(df,
                     data == "meta" &
                     adjust %in% c("IPTW","MW") &
                     approx != "Efron" &
                     !grepl("truth", as.character(data)))

        gg1 <- ggplot(data = df,
                      mapping = aes(x = coef)) +
            geom_density() +
            scale_x_continuous(limits = c(-1,+1)) +
            facet_grid(outcome ~ adjust) +
            labs(title = sprintf("Scenario %s weighted meta-analysis", scenarioCount)) +
            guides(color = FALSE) +
            theme_clean

        gg2 <- ggplot(data = df,
                      mapping = aes(x = site_NA, y = coef, color = site_NA)) +
            geom_boxplot() +
            scale_x_discrete(drop = FALSE) +
            facet_grid(outcome ~ adjust, scales = "fixed") +
            labs(title = sprintf("Scenario %s weighted meta-analysis estimates by # of NA sites", scenarioCount)) +
        guides(color = FALSE) +
            theme_clean

        gg3 <- ggplot(data = df,
                      mapping = aes(x = site_NA, fill = site_NA)) +
            geom_bar() +
            scale_x_discrete(drop = FALSE) +
            scale_y_continuous(limits = c(0,500)) +
            facet_grid(outcome ~ adjust, scales = "fixed") +
            labs(title = sprintf("Scenario %s weighted meta-analysis # of NA sites", scenarioCount)) +
        guides(color = FALSE) +
            theme_clean

        ## Print
        out <- lapply(list(gg1, gg2, gg3), function(gg) {try(print(gg))})
        ## Check status
        out <- Filter(f = function(elt) {is.error(elt)}, x = out)
        ifelse(length(out) > 0, out, "success")
    })

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
