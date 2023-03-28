library(reshape2)
library(ggforce)

source("Analysis/cohort_model_functions.R")
source("Analysis/analysis_functions.R")


traces <- read.csv("Results/Calibration/latest_traces.csv")

## ---Cohort Setup---------------------------------------------------------------- ##
nreps <- 1000
nsteps <- 121
npeople <- 10000
treatment <- 0
threshold_months <- 9
threshold_changes <- 2

## ---No treatment (main analysis)------------------------------------------------ ##

IBM_sim_sub <- full_simulation(nsteps,npeople,1,traces[,1:6],"one",treatment)
IBM_sim_clin <- full_simulation(nsteps,npeople,0,traces[,1:6],"one", treatment)
IBM_sim <- full_simulation(nsteps,npeople,0.5,traces[,1:6],"one",treatment)
traj_sub <- report_trajectories(IBM_sim_sub, npeople, 5, threshold_months, threshold_changes)
traj_clin <- report_trajectories(IBM_sim_clin, npeople, 5, threshold_months, threshold_changes)

# one model run trajectories
traj_sum_sub <- summarise_trajectories(traj_sub,5)
traj_sum_clin <- summarise_trajectories(traj_clin,5)

# 10,000 model run trajectories (slower)
traj_sum_sub <- median_trajectories(nreps, nsteps, npeople, 1, traces[,1:6], "one", treatment, 5)
traj_sum_clin <- median_trajectories(nreps, nsteps, npeople, 0, traces[,1:6], "one", treatment, 5)

write.csv(traj_sum_sub,"Results/Analysis/latest_trajectories_subclinical.csv")
write.csv(traj_sum_clin,"Results/Analysis/latest_trajectories_clinical.csv")

traj_plot_sub <- plot_trajectories(traj_sum_sub, 5, "") + ggtitle("Subclinical cohort")
traj_plot_clin <- plot_trajectories(traj_sum_clin, 5, "") + ggtitle("Clinical cohort")
traj_plot <- traj_plot_sub + traj_plot_clin + plot_layout(guides = "collect") + plot_annotation(tag_levels = 'A')
sankey_plot <- sankey_nat(IBM_sim)

png(filename = "Figures/subclinical-trajectory-plot.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
traj_plot_sub
dev.off()

png(filename = "Figures/clinical-trajectory-plot.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
traj_plot_clin
dev.off()

png(filename = "Figures/trajectory-plots.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
traj_plot
dev.off()

png(filename = "Figures/sankey-plot.png",
    width = 290,
    height = 180,
    units = "mm",
    res = 600)
sankey_plot
dev.off()


summary_stats_IBM_sub <- rep_IBM_results(nreps, nsteps,npeople,1,traces[,1:6],"one",treatment)
summary_stats_IBM_clin <- rep_IBM_results(nreps, nsteps,npeople,0,traces[,1:6],"one",treatment)
summary_stats_IBM <- rep_IBM_results(nreps, nsteps,npeople,0.5,traces[,1:6],"one",treatment)
summary_stats_IBM_mix <- rep_IBM_results(nreps, nsteps,npeople,1/7,traces[,1:6],"one",treatment,5/7)

dur_inf_dis_med <- round(summary_stats_IBM[2,"duration_infectious"], 2)
dur_inf_dis_lower <- round(summary_stats_IBM[3,"duration_infectious"], 2)
dur_inf_dis_upper <- round(summary_stats_IBM[1,"duration_infectious"], 2)

dur_symp_death_med <- round(summary_stats_IBM[2,"symp_before_death_med"], 2)
dur_symp_death_lower <- round(summary_stats_IBM[3,"symp_before_death_med"], 2)
dur_symp_death_upper <- round(summary_stats_IBM[1,"symp_before_death_med"], 2)

dur_symp_reg_med <- round(summary_stats_IBM[2,"symp_before_recover_med"], 2)
dur_symp_reg_lower <- round(summary_stats_IBM[3,"symp_before_recover_med"], 2)
dur_symp_reg_upper <- round(summary_stats_IBM[1,"symp_before_recover_med"], 2)

write.csv(summary_stats_IBM_sub,"Results/Analysis/latest_IBM_sub.csv")
write.csv(summary_stats_IBM_clin,"Results/Analysis/latest_IBM_clin.csv")
write.csv(summary_stats_IBM,"Results/Analysis/latest_IBM.csv")


## ---With treatment at 0.7 (additional analysis)--------------------------------- ##

summary_stats_IBM_t <- rep_IBM_results(nreps, nsteps,npeople,0.5,traces[,1:6],"one",0.7)
summary_stats_IBM_mix_t <- rep_IBM_results(nreps, nsteps,npeople,1/7,traces[,1:6],"one",0.7,5/7)

dur_sympT_death_med <- round(summary_stats_IBM_t[2,"symp_before_death_med"], 2)
dur_sympT_death_lower <- round(summary_stats_IBM_t[3,"symp_before_death_med"], 2)
dur_sympT_death_upper <- round(summary_stats_IBM_t[1,"symp_before_death_med"], 2)

dur_sympT_treat_med <- round(summary_stats_IBM_t[2,"symp_before_treat_med"], 2)
dur_sympT_treat_lower <- round(summary_stats_IBM_t[3,"symp_before_treat_med"], 2)
dur_sympT_treat_upper <- round(summary_stats_IBM_t[1,"symp_before_treat_med"], 2)

dur_sympT_reg_med <- round(summary_stats_IBM_t[2,"symp_before_recover_med"], 2)
dur_sympT_reg_lower <- round(summary_stats_IBM_t[3,"symp_before_recover_med"], 2)
dur_sympT_reg_upper <- round(summary_stats_IBM_t[1,"symp_before_recover_med"], 2)

write.csv(summary_stats_IBM_t,"Results/Analysis/latest_IBM_t.csv")
