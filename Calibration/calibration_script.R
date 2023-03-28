## ---Packages------------------------------------------------------------------------ ##
library(tidyverse)
library(rbi)
library(rbi.helpers)
library(patchwork)
library(png)
library(binom)
library(bayesplot)
library(ggpubr)


## ---Load Files/Data----------------------------------------------------------------- ##
data_file = "included_data.csv"
prep_file = "results_data_prep.R"
libbi_model = "results_model.bi"  

source(paste("Calibration/R", prep_file, sep = "/"))
libbi_file <- paste("Calibration/LibBi", libbi_model, sep = "/")
TB_model <- rbi::bi_model(file = libbi_file)
inits <- list(min_sub = 0,
              sub_min = 0,        
              sub_clin = 0, 
              clin_sub = 0) 


## ---Run LibBi Fit------------------------------------------------------------------- ##
initial_fit <- rbi::sample(TB_model, target = "posterior", nsamples = 1000, nparticles = 1, 
                           init = inits, obs = obs, verbose = FALSE, input = input,
                           start_time = 0, end_time = 20, noutputs = 240)

adapted <- initial_fit %>%
    adapt_proposal(min = 0.25, max = 0.35, adapt = "both",
                   max_iter = 20, truncate = TRUE, verbose = FALSE)

posterior <- adapted %>%
    rbi::sample(nsamples = 10000, verbose = FALSE)

## ---Output Data--------------------------------------------------------------------- ##
param_table <- summary(posterior, type = "param", quantiles = c(0.025, 0.975))

traj_sum <- summary(posterior, type = "state", quantiles = c(0.025, 0.975))

fitting_traj <- traj_sum[which(traj_sum$var %in% c("sub_over_clin_p",
                                                   "min_over_infect_p",
                                                   "duration") &
                                   traj_sum$time == 2),]

param_table <- rbind(param_table, fitting_traj[,-2])

traces = get_traces(posterior)

write.csv(traces, file = "Results/Calibration/latest_traces.csv", row.names = FALSE)
write.csv(traj_sum, file = "Results/Calibration/latest_trajectories.csv", row.names = FALSE)
write.csv(param_table, file = "Results/Calibration/latest_param_table.csv", row.names = FALSE)

## ---Plots--------------------------------------------------------------------------- ##
source("Calibration/R/calibration_plots.R")


trace_plot = mcmc_trace(traces)
corr_plot = mcmc_pairs(traces, diag_fun = "dens")
hist_plot = mcmc_hist(traces)

traj_sum[,2:8] = lapply(traj_sum[,2:8],as.numeric)

fit_plots = total_fit_plots(data, traj_sum)

png(filename = "Figures/Calibration-plots.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
fit_plots
dev.off()

