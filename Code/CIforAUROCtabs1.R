################################################################################
# 
# Calculation of Forecast Intervals for the Area Under the ROC Curve
# 
# Lee Morin, Ph.D.
# Adjunct Assistant Professor
# Queen's University
# 
# January 16, 2018
# 
################################################################################
# 
# This program calculates tables summarizing simulation results of 
# forecast intervals for the area under the ROC curve, when the population is evolving.
# 
# 
# 
################################################################################


################################################################################
# Setup Workspace and Load Libraries
################################################################################

# Clear workspace.
rm(list=ls(all=TRUE))

# Set working directory.
# wd_path <- '/home/ec2-user/CIforAUC' # On AWS
# wd_path <- 'C:/Users/iky155/Documents/Fisher/CIforAUC' 
wd_path <- 'C:/Users/iky155/Documents/Repos/CIforAUROC_copy' 
setwd(wd_path)

# Set path to data folder.
doc_path <- sprintf('%s', wd_path)

# Set path to folder of simulation results.
sim_path <- sprintf('%s/data/sims', wd_path)

# Set path to figure folder.
fig_path <- sprintf('%s/Slides/Figs', wd_path)

# Set path to table folder.
tab_path <- sprintf('%s/Slides/Tabs', wd_path)

# Set path to data folder.
# fun_path <- sprintf('%s/function', wd_path)
fun_path <- sprintf('%s', wd_path)

# Load data.table package for efficient use of databases.
library(data.table)

# Load package for weighted histograms.
library(plotrix)

# Load required libraries.
bench_file_name <- 'CIforAUROCbench1.R'
bench_path_file_name <- sprintf('%s/%s', fun_path, bench_file_name)
source(bench_path_file_name)
null_file_name <- 'CIforAUROCnull2.R'
null_path_file_name <- sprintf('%s/%s', fun_path, null_file_name)
source(null_path_file_name)


################################################################################
# Set Parameters and Load Data
################################################################################



#--------------------------------------------------------------------------------
# Load table of simulated of confidence intervals
#--------------------------------------------------------------------------------


# Read the simulation results back in.
# sim_file_name <- 'sim_1.csv' # 120 reps.
sim_file_name <- 'sim_2.csv' # 1000 reps.
sim_path_file_name <- sprintf('%s/%s', sim_path, sim_file_name)
sim_results <- fread(sim_path_file_name)

summary(sim_results)
# table(sim_results[, V1])
nrow(sim_results)

# sim_results[, sim_num := as.numeric(V1)]


hist(sim_results[, dist_L_F])
hist(sim_results[, dist_H_F])


#--------------------------------------------------------------------------------
# Set parameters to match loaded data
#--------------------------------------------------------------------------------


# Number of replications.
# num_reps <- 120
num_reps <- 1000


# Set sample sizes.
n_x_list <- c(rep(1000, 4))
n_y_list <- n_x_list*0.10

# Set values for true AUROC by regime.
auroc_L_list <- c(0.68, 0.65, 0.75, 0.70)
auroc_H_list <- c(0.72, 0.75, 0.80, 0.70)





################################################################################
# Create tables of coverage rates.
################################################################################

coverage_dt <- data.table(n_x = n_x_list, 
                          n_y = n_y_list, 
                          auroc_L = auroc_L_list, 
                          auroc_H = auroc_H_list)

# model_num <- 1
for (model_num in 1:nrow(coverage_dt)) {
  
  # Extract parameters.
  model_n_x <- coverage_dt[model_num, n_x]
  model_n_y <- coverage_dt[model_num, n_y]
  model_auroc_L <- coverage_dt[model_num, auroc_L]
  model_auroc_H <- coverage_dt[model_num, auroc_H]
  
  # Select results corresponding to model parameters.
  sel_results <- which(sim_results[, n_x] == model_n_x & 
                         sim_results[, n_y] == model_n_y & 
                         sim_results[, auroc_L] == model_auroc_L & 
                         sim_results[, auroc_H] == model_auroc_H)
  
  # Average estimates.
  coverage_dt[model_num, 'mean_A_hat_L'] <- sim_results[sel_results, mean(A_hat_L)]
  coverage_dt[model_num, 'mean_A_hat_H'] <- sim_results[sel_results, mean(A_hat_H)]
  coverage_dt[model_num, 'mean_A_hat_F'] <- sim_results[sel_results, mean(A_hat_F)]
  
  
  # Coverage rates for true values (prediction intervals).
  coverage_dt[model_num, 'auroc_L_in_ci_pred'] <- sum(sim_results[sel_results, ci_l_pred] < coverage_dt[model_num, auroc_L] &
                                                   sim_results[sel_results, ci_u_pred] > coverage_dt[model_num, auroc_L]) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci_pred'] <- sum(sim_results[sel_results, ci_l_pred] < coverage_dt[model_num, auroc_H] &
                                                   sim_results[sel_results, ci_u_pred] > coverage_dt[model_num, auroc_H]) / num_reps
  
  # Coverage rates for true values (confidence intervals).
  coverage_dt[model_num, 'auroc_L_in_ci_binorm'] <- sum(sim_results[sel_results, ci_l_binorm] < coverage_dt[model_num, auroc_L] &
                                                          sim_results[sel_results, ci_u_binorm] > coverage_dt[model_num, auroc_L]) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci_binorm'] <- sum(sim_results[sel_results, ci_l_binorm] < coverage_dt[model_num, auroc_H] &
                                                          sim_results[sel_results, ci_u_binorm] > coverage_dt[model_num, auroc_H]) / num_reps
  coverage_dt[model_num, 'auroc_L_in_ci_delong'] <- sum(sim_results[sel_results, ci_l_delong] < coverage_dt[model_num, auroc_L] &
                                                          sim_results[sel_results, ci_u_delong] > coverage_dt[model_num, auroc_L]) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci_delong'] <- sum(sim_results[sel_results, ci_l_delong] < coverage_dt[model_num, auroc_H] &
                                                          sim_results[sel_results, ci_u_delong] > coverage_dt[model_num, auroc_H]) / num_reps
  coverage_dt[model_num, 'auroc_L_in_ci_boot'] <- sum(sim_results[sel_results, ci_l_boot] < coverage_dt[model_num, auroc_L] &
                                                        sim_results[sel_results, ci_u_boot] > coverage_dt[model_num, auroc_L]) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci_boot'] <- sum(sim_results[sel_results, ci_l_boot] < coverage_dt[model_num, auroc_H] &
                                                        sim_results[sel_results, ci_u_boot] > coverage_dt[model_num, auroc_H]) / num_reps
  coverage_dt[model_num, 'auroc_L_in_ci_upper'] <- sum(sim_results[sel_results, ci_l_upper] < coverage_dt[model_num, auroc_L] &
                                                        sim_results[sel_results, ci_u_upper] > coverage_dt[model_num, auroc_L]) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci_upper'] <- sum(sim_results[sel_results, ci_l_upper] < coverage_dt[model_num, auroc_H] &
                                                        sim_results[sel_results, ci_u_upper] > coverage_dt[model_num, auroc_H]) / num_reps
  
  
  
  
  # Coverage rates for estimated values (estimated with variability).
  coverage_dt[model_num, 'A_hat_L_in_ci_pred'] <- 0
  coverage_dt[model_num, 'A_hat_H_in_ci_pred'] <- 0
  coverage_dt[model_num, 'A_hat_F_in_ci_pred'] <- 0
  
  coverage_dt[model_num, 'A_hat_L_in_ci_binorm'] <- 0
  coverage_dt[model_num, 'A_hat_H_in_ci_binorm'] <- 0
  coverage_dt[model_num, 'A_hat_L_in_ci_delong'] <- 0
  coverage_dt[model_num, 'A_hat_H_in_ci_delong'] <- 0
  coverage_dt[model_num, 'A_hat_L_in_ci_boot'] <- 0
  coverage_dt[model_num, 'A_hat_H_in_ci_boot'] <- 0
  coverage_dt[model_num, 'A_hat_L_in_ci_upper'] <- 0
  coverage_dt[model_num, 'A_hat_H_in_ci_upper'] <- 0
  
  
  
  A_hat_L_list <- sim_results[sel_results, A_hat_L]
  A_hat_H_list <- sim_results[sel_results, A_hat_H]
  A_hat_F_list <- sim_results[sel_results, A_hat_F]
  
  for (rep_num in 1:num_reps) {
    
    # Correct forecast rates (prediction intervals).
    coverage_dt[model_num, 'A_hat_L_in_ci_pred'] <- coverage_dt[model_num, 'A_hat_L_in_ci_pred'] + 
      sum(sim_results[sel_results, ci_l_pred] < A_hat_L_list[rep_num] &
            sim_results[sel_results, ci_u_pred] > A_hat_L_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_H_in_ci_pred'] <- coverage_dt[model_num, 'A_hat_H_in_ci_pred'] + 
      sum(sim_results[sel_results, ci_l_pred] < A_hat_H_list[rep_num] &
            sim_results[sel_results, ci_u_pred] > A_hat_H_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_F_in_ci_pred'] <- coverage_dt[model_num, 'A_hat_F_in_ci_pred'] + 
      sum(sim_results[sel_results, ci_l_pred] < A_hat_F_list[rep_num] &
            sim_results[sel_results, ci_u_pred] > A_hat_F_list[rep_num]) / num_reps

    
    
    
    # Correct forecast rates (confidence intervals).
    coverage_dt[model_num, 'A_hat_L_in_ci_binorm'] <- coverage_dt[model_num, 'A_hat_L_in_ci_binorm'] + 
      sum(sim_results[sel_results, ci_l_binorm] < A_hat_L_list[rep_num] &
            sim_results[sel_results, ci_u_binorm] > A_hat_L_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_H_in_ci_binorm'] <- coverage_dt[model_num, 'A_hat_H_in_ci_binorm'] + 
      sum(sim_results[sel_results, ci_l_binorm] < A_hat_H_list[rep_num] &
            sim_results[sel_results, ci_u_binorm] > A_hat_H_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_L_in_ci_delong'] <- coverage_dt[model_num, 'A_hat_L_in_ci_delong'] + 
      sum(sim_results[sel_results, ci_l_delong] < A_hat_L_list[rep_num] &
            sim_results[sel_results, ci_u_delong] > A_hat_L_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_H_in_ci_delong'] <- coverage_dt[model_num, 'A_hat_H_in_ci_delong'] + 
      sum(sim_results[sel_results, ci_l_delong] < A_hat_H_list[rep_num] &
            sim_results[sel_results, ci_u_delong] > A_hat_H_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_L_in_ci_boot'] <- coverage_dt[model_num, 'A_hat_L_in_ci_boot'] + 
      sum(sim_results[sel_results, ci_l_boot] < A_hat_L_list[rep_num] &
            sim_results[sel_results, ci_u_boot] > A_hat_L_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_H_in_ci_boot'] <- coverage_dt[model_num, 'A_hat_H_in_ci_boot'] + 
      sum(sim_results[sel_results, ci_l_boot] < A_hat_H_list[rep_num] &
            sim_results[sel_results, ci_u_boot] > A_hat_H_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_L_in_ci_upper'] <- coverage_dt[model_num, 'A_hat_L_in_ci_upper'] + 
      sum(sim_results[sel_results, ci_l_upper] < A_hat_L_list[rep_num] &
            sim_results[sel_results, ci_u_upper] > A_hat_L_list[rep_num]) / num_reps
    
    coverage_dt[model_num, 'A_hat_H_in_ci_upper'] <- coverage_dt[model_num, 'A_hat_H_in_ci_upper'] + 
      sum(sim_results[sel_results, ci_l_upper] < A_hat_H_list[rep_num] &
            sim_results[sel_results, ci_u_upper] > A_hat_H_list[rep_num]) / num_reps
    
    
    
        
  }
  
  
  
  
  
}

colnames(coverage_dt)


coverage_dt[, c('A_hat_L_in_ci_pred', 'A_hat_H_in_ci_pred', 'A_hat_F_in_ci_pred')] <- 
  coverage_dt[, c('A_hat_L_in_ci_pred', 'A_hat_H_in_ci_pred', 'A_hat_F_in_ci_pred')] / num_reps
coverage_dt[, c('A_hat_L_in_ci_binorm', 'A_hat_H_in_ci_binorm', 'A_hat_F_in_ci_binorm')] <- 
  coverage_dt[, c('A_hat_L_in_ci_binorm', 'A_hat_H_in_ci_binorm')] / num_reps
coverage_dt[, c('A_hat_L_in_ci_delong', 'A_hat_H_in_ci_delong')] <- 
  coverage_dt[, c('A_hat_L_in_ci_delong', 'A_hat_H_in_ci_delong')] / num_reps
coverage_dt[, c('A_hat_L_in_ci_boot', 'A_hat_H_in_ci_boot')] <- 
  coverage_dt[, c('A_hat_L_in_ci_boot', 'A_hat_H_in_ci_boot')] / num_reps
coverage_dt[, c('A_hat_L_in_ci_upper', 'A_hat_H_in_ci_upper')] <- 
  coverage_dt[, c('A_hat_L_in_ci_upper', 'A_hat_H_in_ci_upper')] / num_reps


coverage_dt



# Save the simulation results after summary calculation.
# sim_file_name <- 'sim_1.csv' # 120 reps.
sim_out_file_name <- 'sim_2_summary.csv' # 1000 reps.
sim_out_path_file_name <- sprintf('%s/%s', sim_path, sim_out_file_name)
write.csv(coverage_dt, sim_out_path_file_name)


################################################################################
# Create tables for LaTeX documents.
################################################################################




# Calculate coverage rates for each set of confidence intervals.
coverage_dt[, 'forecast_cov'] <- (coverage_dt[, auroc_L_in_ci_pred] + coverage_dt[, auroc_H_in_ci_pred])/2
coverage_dt[, 'forecast_fore'] <- (coverage_dt[, A_hat_L_in_ci_pred] + coverage_dt[, A_hat_H_in_ci_pred])/2

coverage_dt[, 'binorm_cov'] <- (coverage_dt[, auroc_L_in_ci_binorm] + coverage_dt[, auroc_H_in_ci_binorm])/2
coverage_dt[, 'binorm_fore'] <- (coverage_dt[, A_hat_L_in_ci_binorm] + coverage_dt[, A_hat_H_in_ci_binorm])/2

coverage_dt[, 'delong_cov'] <- (coverage_dt[, auroc_L_in_ci_delong] + coverage_dt[, auroc_H_in_ci_delong])/2
coverage_dt[, 'delong_fore'] <- (coverage_dt[, A_hat_L_in_ci_delong] + coverage_dt[, A_hat_H_in_ci_delong])/2

coverage_dt[, 'boot_cov'] <- (coverage_dt[, auroc_L_in_ci_boot] + coverage_dt[, auroc_H_in_ci_boot])/2
coverage_dt[, 'boot_fore'] <- (coverage_dt[, A_hat_L_in_ci_boot] + coverage_dt[, A_hat_H_in_ci_boot])/2

coverage_dt[, 'upper_cov'] <- (coverage_dt[, auroc_L_in_ci_upper] + coverage_dt[, auroc_H_in_ci_upper])/2
coverage_dt[, 'upper_fore'] <- (coverage_dt[, A_hat_L_in_ci_upper] + coverage_dt[, A_hat_H_in_ci_upper])/2





# Assemble into coverage rates by sets of parameters.
stat_name_list <- c('Bi-normal', 'DeLong et. al.', 'Bootstrap', 'Upper Bound', 'Forecast')
stat_tag_list <- c('binorm', 'delong', 'boot', 'upper', 'forecast')
num_stats <- length(stat_tag_list)
for (model_num in 1:nrow(coverage_dt)) {
  
  # Generate table of coverage rates and correct forecast rates.
  coverage_row <- data.table(n = rep(coverage_dt[model_num, n_x + n_y], num_stats), 
                             auroc_L = rep(coverage_dt[model_num, auroc_L], num_stats), 
                             auroc_H = rep(coverage_dt[model_num, auroc_H], num_stats), 
                             method = stat_name_list, 
                             coverage = coverage_dt[model_num, sprintf('%s_cov', stat_tag_list), with = FALSE],
                             forecast = coverage_dt[model_num, sprintf('%s_fore', stat_tag_list), with = FALSE])
  
  # Print out the relevant matrix.
  print(sprintf('Coverage rates for n = %d, A_L = %f, A_H = %f', 
                coverage_dt[model_num, n_x + n_y], 
                coverage_dt[model_num, auroc_L], 
                coverage_dt[model_num, auroc_H]))
  print(coverage_row[, c('method', 'coverage', 'forecast'), with = FALSE])
  
  
}









################################################################################
# End
################################################################################

