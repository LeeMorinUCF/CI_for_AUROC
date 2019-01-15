################################################################################
# 
# Calculation of COnfidence Intervals for the Area Under the ROC Curve
# 
# Lee Morin, Ph.D.
# Adjunct Assistant Professor
# Queen's University
# 
# January 6, 2018
# 
################################################################################
# 
# This program calculates a series of confidence intervals for the area under
# the ROC curve, when the population is evolving.
# 
# This version introduces uniform notation for L,M,H AUROC.
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

# Set path to function folder.
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
# Set Parameters and Generate Data
################################################################################

#--------------------------------------------------------------------------------
# True AUROC values.
#--------------------------------------------------------------------------------

# Set tags for regimes.
reg_tag <- c('L', 'H')

# Set values for true AUROC by regime.
auroc_L <- 0.68
auroc_H <- 0.72
auroc_set <- c(auroc_L, auroc_H)

# Solve for binormal parameters at midpoint.
binorm_spec_list_L <- binorm_spec(auroc_L)
binorm_spec_list_H <- binorm_spec(auroc_H)



#--------------------------------------------------------------------------------
# Binormal model for main population.
#--------------------------------------------------------------------------------


# Typical unbalanced case (requires re-weighting).
n_x <- 1000
n_y <- 100
# n_x <- 500
# n_y <- 50


# Binormal Model for low-AUROC population.
binorm_L <- data.table(binorm_gen(n_x, n_y, 
                                  mu_x = binorm_spec_list_L$mu_x, 
                                  mu_y = binorm_spec_list_L$mu_y, 
                                  sigma_x = binorm_spec_list_L$sigma_x, 
                                  sigma_y = binorm_spec_list_L$sigma_y))
binorm_L[, segment := 'L']
binorm_L[, weights_seg := weights]

# Binormal Model for high-AUROC population.
binorm_H <- data.table(binorm_gen(n_x, n_y, 
                                  mu_x = binorm_spec_list_H$mu_x, 
                                  mu_y = binorm_spec_list_H$mu_y, 
                                  sigma_x = binorm_spec_list_H$sigma_x, 
                                  sigma_y = binorm_spec_list_H$sigma_y))
binorm_H[, segment := 'H']
binorm_H[, weights_seg := weights]

# Collect into a single dataset.
binorm_dt <- rbind(binorm_L, binorm_H)
binorm_dt[, weights := weights/2]

# Trust but verify.
# summary(binorm_dt)
# summary(binorm_dt[outcomes == TRUE, ])
# summary(binorm_dt[outcomes == FALSE, ])



# Test AUROC for each population.
A_hat_L <- auc_weighted_loops(scores = binorm_dt[segment == 'L', scores], 
                              outcomes = binorm_dt[segment == 'L', outcomes], 
                              weights = binorm_dt[segment == 'L', weights_seg])
A_hat_H <- auc_weighted_loops(scores = binorm_dt[segment == 'H', scores], 
                              outcomes = binorm_dt[segment == 'H', outcomes], 
                              weights = binorm_dt[segment == 'H', weights_seg])
A_hat_F <- auc_weighted_loops(scores = binorm_dt[, scores], 
                              outcomes = binorm_dt[, outcomes], 
                              weights = binorm_dt[, weights])

c(A_hat_L, A_hat_F, A_hat_H)


################################################################################
# Calculating Distances between Observed Distributions
################################################################################


# Quantiles for histogram for calculating distance.
num_quantiles <- 10

# Obtain the quantiles of the observed scores by outcome category.
quantile_df <- distn_quantiles(binorm_dt, num_quantiles)
# Quantiles calculated from the full distribution.

# Calculate observed distances between distributions.

# Obtain a table of relative frequencies for the observed scores.
distn_table_F <- quantile_assign(binorm_dt, quantile_cuts = quantile_df, type = 'observed')
distn_table_L <- quantile_assign(binorm_dt[segment == 'L'], 
                                     quantile_cuts = quantile_df, type = 'observed')
distn_table_H <- quantile_assign(binorm_dt[segment == 'H'], 
                                     quantile_cuts = quantile_df, type = 'observed')


# Pass the tabulated results to calculate distance from each segment to full dataset.
dist_L_F <- auroc_dist(outcomes = distn_table_F[, 'outcomes'], 
                           weights_1 = distn_table_F[, 'distn_obs'], 
                           weights_2 = distn_table_L[, 'distn_obs'], 
                           metric = 'kld0_joint')
dist_H_F <- auroc_dist(outcomes = distn_table_F[, 'outcomes'], 
                           weights_1 = distn_table_F[, 'distn_obs'], 
                           weights_2 = distn_table_H[, 'distn_obs'], 
                           metric = 'kld0_joint')

c(dist_L_F, dist_H_F)




#--------------------------------------------------------------------------------
# Finding AUROC range for given distance values.
#--------------------------------------------------------------------------------

summary(binorm_dt)

# Initialize with required parameters.
dist_avg_F <- mean(c(dist_L_F, dist_H_F))

A_lim_est_L <- A_hat_F - 6*(A_hat_F - A_hat_L)
A_lim_est_H <- A_hat_F + 4*(A_hat_H - A_hat_F)


# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
tol <- 0.00001
max_iter <- 20

# Solve for the upper root.
ci_u_solve <- uniroot(f = dist_A_lim_diff, 
                      interval = c(A_hat_F, A_lim_est_H), 
                      tol = .Machine$double.eps^0.25, maxiter = 100, 
                      D_bar = dist_avg_F, 
                      dt_in = binorm_dt, 
                      metric = 'kld0_joint', 
                      num_quantiles, min_eps = 10^(-4), 
                      eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                      display = FALSE, ylim = NA)

# Solve for the lower root.
ci_l_solve <- uniroot(f = dist_A_lim_diff, 
                      interval = c(A_lim_est_L, A_hat_F + 0.01), 
                      tol = .Machine$double.eps^0.25, maxiter = 100, 
                      D_bar = dist_avg_F, 
                      dt_in = binorm_dt, 
                      metric = 'kld0_joint', 
                      num_quantiles, min_eps = 10^(-4), 
                      eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                      display = FALSE, ylim = NA)




################################################################################
# Put it together into a simulation
################################################################################


#--------------------------------------------------------------------------------
# Set parameters
#--------------------------------------------------------------------------------

# Number of replications.
num_reps <- 10

# Quantiles for histogram for calculating distance.
num_quantiles <- 10

# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
tol <- 0.00001
max_iter <- 20

# Set tags for regimes.
reg_tag <- c('L', 'H')

# Typical unbalanced case (requires re-weighting).
n_x <- 1000
n_y <- 100

# Set values for true AUROC by regime.
auroc_L <- 0.68
auroc_H <- 0.72
auroc_set <- c(auroc_L, auroc_H)


# Solve for binormal parameters at midpoint.
binorm_spec_list_L <- binorm_spec(auroc_L)
binorm_spec_list_H <- binorm_spec(auroc_H)



#--------------------------------------------------------------------------------
# Run simulation.
#--------------------------------------------------------------------------------

# Initialize matrices for storage of confidence intervals.
conf_int_dt <- data.table(A_hat_F = numeric(num_reps),
                          A_hat_L = numeric(num_reps),
                          A_hat_H = numeric(num_reps),
                          dist_L_F = numeric(num_reps),
                          dist_H_F = numeric(num_reps),
                          cu_l = numeric(num_reps),
                          cu_u = numeric(num_reps))

# rep_num <- 1
for (rep_num in 1:num_reps) {
  
  #--------------------------------------------------------------------------------
  # Generate data.
  #--------------------------------------------------------------------------------
  
  # Binormal Model for low-AUROC population.
  binorm_L <- data.table(binorm_gen(n_x, n_y, 
                                    mu_x = binorm_spec_list_L$mu_x, 
                                    mu_y = binorm_spec_list_L$mu_y, 
                                    sigma_x = binorm_spec_list_L$sigma_x, 
                                    sigma_y = binorm_spec_list_L$sigma_y))
  binorm_L[, segment := 'L']
  binorm_L[, weights_seg := weights]
  
  # Binormal Model for high-AUROC population.
  binorm_H <- data.table(binorm_gen(n_x, n_y, 
                                    mu_x = binorm_spec_list_H$mu_x, 
                                    mu_y = binorm_spec_list_H$mu_y, 
                                    sigma_x = binorm_spec_list_H$sigma_x, 
                                    sigma_y = binorm_spec_list_H$sigma_y))
  binorm_H[, segment := 'H']
  binorm_H[, weights_seg := weights]
  
  # Collect into a single dataset.
  binorm_dt <- rbind(binorm_L, binorm_H)
  binorm_dt[, weights := weights/2]
  
  #--------------------------------------------------------------------------------
  # Set parameters for solution of confidence bounds.
  #--------------------------------------------------------------------------------
  
  
  # Test AUROC for each population.
  A_hat_L <- auc_weighted_loops(scores = binorm_dt[segment == 'L', scores], 
                                outcomes = binorm_dt[segment == 'L', outcomes], 
                                weights = binorm_dt[segment == 'L', weights_seg])
  A_hat_H <- auc_weighted_loops(scores = binorm_dt[segment == 'H', scores], 
                                outcomes = binorm_dt[segment == 'H', outcomes], 
                                weights = binorm_dt[segment == 'H', weights_seg])
  A_hat_F <- auc_weighted_loops(scores = binorm_dt[, scores], 
                                outcomes = binorm_dt[, outcomes], 
                                weights = binorm_dt[, weights])
  
  #--------------------------------------------------------------------------------
  # Calculate Distances between Observed Distributions
  #--------------------------------------------------------------------------------
  
  # Obtain the quantiles of the observed scores by outcome category.
  quantile_df <- distn_quantiles(binorm_dt, num_quantiles)
  # Quantiles calculated from the full distribution.
  
  # Calculate observed distances between distributions.
  
  # Obtain a table of relative frequencies for the observed scores.
  distn_table_F <- quantile_assign(binorm_dt, quantile_cuts = quantile_df, type = 'observed')
  distn_table_L <- quantile_assign(binorm_dt[segment == 'L'], 
                                   quantile_cuts = quantile_df, type = 'observed')
  distn_table_H <- quantile_assign(binorm_dt[segment == 'H'], 
                                   quantile_cuts = quantile_df, type = 'observed')
  
  
  # Pass the tabulated results to calculate distance from each segment to full dataset.
  dist_L_F <- auroc_dist(outcomes = distn_table_F[, 'outcomes'], 
                         weights_1 = distn_table_F[, 'distn_obs'], 
                         weights_2 = distn_table_L[, 'distn_obs'], 
                         metric = 'kld0_joint')
  dist_H_F <- auroc_dist(outcomes = distn_table_F[, 'outcomes'], 
                         weights_1 = distn_table_F[, 'distn_obs'], 
                         weights_2 = distn_table_H[, 'distn_obs'], 
                         metric = 'kld0_joint')
  
  
  # Initialize with required parameters.
  dist_avg_F <- mean(c(dist_L_F, dist_H_F))
  
  
  # Solve for the upper root.
  ci_u_solve <- uniroot(f = dist_A_lim_diff, 
                        interval = c(A_hat_F, 1.0), 
                        tol = .Machine$double.eps^0.25, maxiter = 100, 
                        D_bar = dist_avg_F, 
                        dt_in = binorm_dt, 
                        metric = 'kld0_joint', 
                        num_quantiles, min_eps = 10^(-4), 
                        eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                        display = FALSE, ylim = NA)
  
  # Solve for the upper root.
  ci_l_solve <- uniroot(f = dist_A_lim_diff, 
                        interval = c(0.5, A_hat_F), 
                        tol = .Machine$double.eps^0.25, maxiter = 100, 
                        D_bar = dist_avg_F, 
                        dt_in = binorm_dt, 
                        metric = 'kld0_joint', 
                        num_quantiles, min_eps = 10^(-4), 
                        eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                        display = FALSE, ylim = NA)
  
  
  
  #--------------------------------------------------------------------------------
  # Store values of confidence intervals
  #--------------------------------------------------------------------------------
  
  conf_int_dt[rep_num, 'A_hat_F'] <- A_hat_F
  conf_int_dt[rep_num, 'A_hat_L'] <- A_hat_L
  conf_int_dt[rep_num, 'A_hat_H'] <- A_hat_H
  conf_int_dt[rep_num, 'dist_L_F'] <- dist_L_F
  conf_int_dt[rep_num, 'dist_H_F'] <- dist_H_F
  conf_int_dt[rep_num, 'cu_l'] <- ci_l_solve$root
  conf_int_dt[rep_num, 'cu_u'] <- ci_u_solve$root
  
  print(conf_int_dt[rep_num, ])
  
}

print(conf_int_dt)


################################################################################
# Create tables of coverage rates.
################################################################################


# Read the simulation results back in.
sim_file_name <- 'sim_1.csv'
sim_path_file_name <- sprintf('%s/%s', sim_path, sim_file_name)
sim_results <- fread(sim_path_file_name)

summary(sim_results)
# table(sim_results[, V1])
nrow(sim_results)

sim_results[, sim_num := as.numeric(V1)]


hist(sim_results[, dist_L_F])
hist(sim_results[, dist_H_F])


# Compared to true values:



# Direct approach: Loop/lapply over Intervals.

# Initialize a matrix.




# Number of replications.
num_reps <- 120


# Set sample sizes.
n_x_list <- c(rep(1000, 4))
n_y_list <- n_x_list*0.10

# Set values for true AUROC by regime.
auroc_L_list <- c(0.68, 0.65, 0.75, 0.70)
auroc_H_list <- c(0.72, 0.75, 0.80, 0.70)


coverage_dt <- data.table(n_x = n_x_list, 
                          n_y = n_y_list, 
                          auroc_L = auroc_L_list, 
                          auroc_H = auroc_H_list)

model_num <- 1
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
  
  
  # Coverage rates.
  coverage_dt[model_num, 'auroc_L_in_ci'] <- sum(sim_results[sel_results, ci_l] < coverage_dt[model_num, 'auroc_L'] &
                                                   sim_results[sel_results, ci_u] > coverage_dt[model_num, 'auroc_L']) / num_reps
  coverage_dt[model_num, 'auroc_H_in_ci'] <- sum(sim_results[sel_results, ci_l] < coverage_dt[model_num, 'auroc_H'] &
                                                   sim_results[sel_results, ci_u] > coverage_dt[model_num, 'auroc_H']) / num_reps
  
  # coverage_dt[model_num, 'A_hat_F_in_ci'] <- sum(sim_results[sel_results, ci_l] < coverage_dt[model_num, 'A_hat_F'] &
  #                                                  sim_results[sel_results, ci_u] > coverage_dt[model_num, 'A_hat_F']) / num_reps
  
  
  
  
  
}


coverage_dt


# Compute (one-sided) coverage rates using AUROC.

# Create a temporary data table.


# Calculate 





################################################################################
# Previous Version with Reweighting to Determine DGP
################################################################################


# #--------------------------------------------------------------------------------
# # Finding weights for low-AUROC population.
# #--------------------------------------------------------------------------------
# 
# # Set parameters for re-weighting toward the null hypothesis.
# eta <- 4
# tol <- 0.00000001
# max_iter <- 20
# ylim <- c(0, 0.4)
# 
# 
# 
# A_0 <- auroc_set[1]
# 
# adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
#                display = TRUE, ylim = ylim)
# 
# 
# 
# # Record new weights for low-AUROC distribution.
# binorm_dt[, weights_L := opt_weights]
# 
# # Test AUROC for main population.
# A_hat_L <- auc_weighted_loops(scores = binorm_dt[, scores], 
#                                 outcomes = binorm_dt[, outcomes], 
#                                 weights = binorm_dt[, weights_L])
# A_hat_L
# 
# 
# #--------------------------------------------------------------------------------
# # Finding weights for high-AUROC population.
# #--------------------------------------------------------------------------------
# 
# A_0 <- auroc_set[3]
# 
# adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
#                display = TRUE, ylim = ylim)
# 
# 
# 
# # Record new weights for low-AUROC distribution.
# binorm_dt[, weights_H := opt_weights]
# 
# # Test AUROC for main population.
# A_hat_H <- auc_weighted_loops(scores = binorm_dt[, scores], 
#                                  outcomes = binorm_dt[, outcomes], 
#                                  weights = binorm_dt[, weights_H])
# A_hat_H


#--------------------------------------------------------------------------------
# Generate full dataset by stacking segment datasets
#--------------------------------------------------------------------------------

# # Add a segment variable and initialize the full distribution.
# binorm_dt[, segment := 2]
# binorm_dt[, weights_seg := weights_M]
# # summary(binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE])
# binorm_full <- binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE]
# 
# # Append low-AUROC dataset.
# binorm_dt[, segment := 1]
# binorm_dt[, weights_seg := weights_L]
# # summary(binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE])
# binorm_full <- rbind(binorm_full, 
#                      binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE])
# 
# # Append high-AUROC dataset.
# binorm_dt[, segment := 3]
# binorm_dt[, weights_seg := weights_H]
# # summary(binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE])
# binorm_full <- rbind(binorm_full, 
#                      binorm_dt[, c('segment', 'outcomes', 'scores', 'weights_seg'), with = FALSE])
# 
# summary(binorm_full)
# 
# 
# A_hat_full <- auc_weighted_loops(scores = binorm_full[, scores], 
#                                  outcomes = binorm_full[, outcomes], 
#                                  weights = binorm_full[, weights_seg/3])
# A_hat_full
# 



################################################################################
# Estimate Distance between Distributions
################################################################################


#--------------------------------------------------------------------------------
# Theoretical quantities using true distributions
#--------------------------------------------------------------------------------


# dist_L_M <- auroc_dist(outcomes = binorm_dt[, outcomes], 
#                        weights_1 = binorm_dt[, weights_L], 
#                        weights_2 = binorm_dt[, weights_M], 
#                        metric = 'kld0_joint', min_eps = 10^(-4))
# 
# dist_M_H <- auroc_dist(outcomes = binorm_dt[, outcomes], 
#                        weights_1 = binorm_dt[, weights_M], 
#                        weights_2 = binorm_dt[, weights_H], 
#                        metric = 'kld0_joint', min_eps = 10^(-4))
# 
# dist_L_H <- auroc_dist(outcomes = binorm_dt[, outcomes], 
#                        weights_1 = binorm_dt[, weights_L], 
#                        weights_2 = binorm_dt[, weights_H], 
#                        metric = 'kld0_joint', min_eps = 10^(-4))
# 
# c(dist_L_M, dist_M_H, dist_L_H)
# 
# 
# 
# # Calculate true average distance under known breaks between segments.
# 
# # Symmetric Intervals (unfounded)
# dist_avg_est <- 2/9*dist_L_M + 2/9*dist_M_H + 2/9*dist_L_H
# dist_avg_est
# 
# # Asymmetric intervals.
# dist_L_M_est <- 2/3*dist_L_M



#--------------------------------------------------------------------------------
# Estimation from realizations
#--------------------------------------------------------------------------------

# A_hat_full is the benchmark

# Reps for selecting samples for calculating distance.
num_dist_reps <- 10
# Quantiles for histogram for calculating distance.
num_quantiles <- 10


#--------------------------------------------------------------------------------
# Draw a realization from the DGP
#--------------------------------------------------------------------------------


# Set the number of realized datasets drawn from the DGP.
sim_dt <- NULL
for (regime_num in 1:length(auroc_set)) {
  
  # Select a sample of the observations within the selected regime.
  regime_sel_wts <- sprintf('weights_%s', reg_tag[regime_num])
  pos_sel_ind <- sample(x = 1:nrow(binorm_dt), size = n_y, replace = TRUE, 
                        prob = binorm_dt[, get(regime_sel_wts)*outcomes])
  neg_sel_ind <- sample(x = 1:nrow(binorm_dt), size = n_x, replace = TRUE, 
                        prob = binorm_dt[, get(regime_sel_wts)*(1 - outcomes)])
  pos_sel_scores <- binorm_dt[pos_sel_ind, scores]
  neg_sel_scores <- binorm_dt[neg_sel_ind, scores]
  sample_dt <- data.table(scores = c(pos_sel_scores, neg_sel_scores),
                          outcomes = c(rep(TRUE, n_y), rep(FALSE, n_x)),
                          weights = c(rep(1/n_y, n_y), rep(1/n_x, n_x)),
                          segment = rep(regime_num, n_x + n_y))
  # Note equal weighting since the observations are observed [sic].
  
  # Append this draw to the simulated dataset.
  sim_dt <- rbind(sim_dt, sample_dt)
  
}
summary(sim_dt)

# Obtain the quantiles of the observed scores by outcome category.
quantile_df <- distn_quantiles(sim_dt, num_quantiles)
# Quantiles calculated from the full distribution.

# Calculate observed distances between distributions.

# Obtain a table of relative frequencies for the observed scores.
distn_table_F_sim <- quantile_assign(sim_dt, quantile_cuts = quantile_df, type = 'observed')
distn_table_L_sim <- quantile_assign(sim_dt[segment == 1], 
                                     quantile_cuts = quantile_df, type = 'observed')
distn_table_M_sim <- quantile_assign(sim_dt[segment == 2], 
                                     quantile_cuts = quantile_df, type = 'observed')
distn_table_H_sim <- quantile_assign(sim_dt[segment == 3], 
                                     quantile_cuts = quantile_df, type = 'observed')


# Pass the tabulated results to calculate distance from each segment to full dataset.
dist_L_F_sim <- auroc_dist(outcomes = distn_table_F_sim[, 'outcomes'], 
                           weights_1 = distn_table_F_sim[, 'distn_obs'], 
                           weights_2 = distn_table_L_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')
dist_M_F_sim <- auroc_dist(outcomes = distn_table_F_sim[, 'outcomes'], 
                           weights_1 = distn_table_F_sim[, 'distn_obs'], 
                           weights_2 = distn_table_M_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')
dist_H_F_sim <- auroc_dist(outcomes = distn_table_F_sim[, 'outcomes'], 
                           weights_1 = distn_table_F_sim[, 'distn_obs'], 
                           weights_2 = distn_table_H_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')



# Pass the tabulated results to calculate distance between pairs of segments.
dist_L_M_sim <- auroc_dist(outcomes = distn_table_L_sim[, 'outcomes'], 
                           weights_1 = distn_table_L_sim[, 'distn_obs'], 
                           weights_2 = distn_table_M_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')
dist_M_H_sim <- auroc_dist(outcomes = distn_table_M_sim[, 'outcomes'], 
                           weights_1 = distn_table_M_sim[, 'distn_obs'], 
                           weights_2 = distn_table_H_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')
dist_L_H_sim <- auroc_dist(outcomes = distn_table_L_sim[, 'outcomes'], 
                           weights_1 = distn_table_L_sim[, 'distn_obs'], 
                           weights_2 = distn_table_H_sim[, 'distn_obs'], 
                           metric = 'kld0_joint')


# Distances from full dataset.
c(dist_L_F_sim, dist_M_F_sim, dist_H_F_sim)

# Distances between pairs of segments in samples.
c(dist_L_M_sim, dist_M_H_sim, dist_L_H_sim)

# Distances between pairs of segments in DGP.
c(dist_L_M, dist_M_H, dist_L_H)
# Note different units.

# GOAL: To match units of distance 
# from 
c(dist_L_F_sim, dist_H_F_sim)
# to 
c(dist_L_M, dist_M_H)

#--------------------------------------------------------------------------------
# Draw a realization from the DGP
#--------------------------------------------------------------------------------


# Initialize totals.
dist_L_M_hat <- 0
dist_M_H_hat <- 0



# Run simulation to estimate distance in each direction.
for (rep_num in 1:num_dist_reps) {
  
  # Select a sample from one of the regimes.
  regime_num <- sample(x = 1:length(auroc_set), size = 1, replace = TRUE)
  regime_sel_wts <- sprintf('weights_%s', reg_tag[regime_num])
  
  # Select a stratified sample of the observations within the selected regime
  
  # Select a sample of the observations within the selected regime.
  pos_sel_ind <- sample(x = 1:nrow(binorm_dt), size = n_y, replace = TRUE, 
                        prob = binorm_dt[, get(regime_sel_wts)*outcomes])
  neg_sel_ind <- sample(x = 1:nrow(binorm_dt), size = n_x, replace = TRUE, 
                        prob = binorm_dt[, get(regime_sel_wts)*(1 - outcomes)])
  pos_sel_scores <- binorm_dt[pos_sel_ind, scores]
  neg_sel_scores <- binorm_dt[neg_sel_ind, scores]
  sample_dt <- data.table(scores = c(pos_sel_scores, neg_sel_scores),
                          outcomes = c(rep(TRUE, n_y), rep(FALSE, n_x)),
                          weights = c(rep(1/n_y, n_y), rep(1/n_x, n_x)))
  
  # sample_dt <- binorm_dt[pos_sel_ind, ':='(scores = scores, outcomes = TRUE, weights = 1/n_y)]
  
  
  # sample_dt <- binorm_dt[sample(x = 1:.N, size = .N, replace = TRUE, prob = get(regime_sel_wts)), 
  #                        c('outcomes', 'scores', regime_sel_wts), with = FALSE]
  # # Append equal weights, as observed in sample.
  # n_x_rep <- sample_dt[outcomes == FALSE, .N]
  # n_x_rep <- sample_dt[outcomes == TRUE, .N]
  
  
  #--------------------------------------------------------------------------------
  # Note that this sample is not stratified.
  #--------------------------------------------------------------------------------
  
  
  # Calculate the AUROC.
  A_hat_rep <- auc_weighted_loops(scores = sample_dt[, scores], 
                                   outcomes = sample_dt[, outcomes], 
                                   weights = sample_dt[, weights])
  
  # Calculate the distance from the full sample.
  
  
  
  # Obtain a table of relative frequencies for the observed scores.
  distn_table <- quantile_assign(sample_dt, quantile_cuts = quantile_df, type = 'observed')
  # Pass the tabulated results to calculate distance.
  auroc_dist_sim <- auroc_dist(outcomes = distn_table[, 'outcomes'], 
                               weights_1 = distn_table[, 'distn_obs'], 
                               weights_2 = distn_table[, 'distn_sim'], 
                               metric = metric)
  
  # Check direction and record appropriately.
  if (A_hat_rep > A_hat_full) {
    
  }
  
  # Record the results for coverage and power analysis.
  
}



#--------------------------------------------------------------------------------
# 
#--------------------------------------------------------------------------------





################################################################################
# End
################################################################################


