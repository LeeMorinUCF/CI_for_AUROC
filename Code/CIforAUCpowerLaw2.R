################################################################################
# 
# Calculation of COnfidence Intervals for the Area Under the ROC Curve
# 
# Lee Morin, Ph.D.
# Assistant Professor
# Department of Economics
# College of Business Administration
# University of Central Florida
# 
# January 14, 2019
# 
################################################################################
# 
# This program calculates the area under the ROC curve, 
# under a power law distribution of scores.
# 
################################################################################



################################################################################
# Setup Workspace and Load Libraries
################################################################################

# Clear workspace.
rm(list=ls(all=TRUE))

# Set working directory.
# wd_path <- '/home/ec2-user/CIforAUC' # On AWS
wd_path <- 'C:/Users/le279259/Documents/Research/CIforAUROC/CI_for_AUROC/Code' # On Windows

setwd(wd_path)


# Set path to function folder.
# fun_path <- sprintf('%s/function', wd_path)
fun_path <- sprintf('%s', wd_path)


# Load required libraries.
bench_file_name <- 'CIforAUROCbench1.R'
bench_path_file_name <- sprintf('%s/%s', fun_path, bench_file_name)
source(bench_path_file_name)
null_file_name <- 'CIforAUROCnull2.R'
null_path_file_name <- sprintf('%s/%s', fun_path, null_file_name)
source(null_path_file_name)

# Load data.table package for efficient use of databases.
library(data.table)


################################################################################
# Set Parameters and Generate Data
################################################################################

# Number of replications for simulation of distribution.
num_sims <- 500

#--------------------------------------------------------------------------------
# Power Law as true model.
#--------------------------------------------------------------------------------


# Typical unbalanced case (requires re-weighting).
# n_x <- 1000
# n_y <- 100
# n_x <- 500
# n_y <- 500
n_x <- 1000
n_y <- 1000

# Set parameters for positive score distribution. 
# y_min <- 2.5
y_min <- 2.0
gamma_y <- 1.2

# Set parameters for negative score distribution. 
x_min <- 2.0
alpha_x <- 2.5

# Select continuous or discrete power law.
# disc_or_cts <- 'disc' # Default
disc_or_cts <- 'cts'

# Specify method of handling ties in ranking for AUROC.
ties_method_list <- c("average", "first", "last", 
                      "random", "max", "min")



################################################################################
# Calculate summary statistics for AUROC Statistics
################################################################################


# Calculate theoretical AUROC.
A_hat_bi_power <- auc_bi_power(x_min, alpha_x, 
                               y_min, gamma_y)
# A_hat_bi_power

# Calculate empirical AUROC and variance from a realization. 
bipower_df <- bipower_gen(n_x, n_y, 
                          x_min, alpha_x, 
                          y_min, gamma_y)

A_hat_rank <- auc_rank(scores = bipower_df[, 'scores'], 
                       outcomes = bipower_df[, 'outcomes'], 
                       ties_method = 'average')
# A_hat_rank
A_hat_var <- var_auroc_fast_calc(dt_in = as.data.table(bipower_df))
# A_hat_var

# Calculate confidence interval by inverting z-statistic. 
A_hat_CI <- confidence_interval(mean = A_hat_bi_power, 
                                var = A_hat_var, 
                                alpha = 0.05)
# A_hat_CI

# # Plot the pair of distributions.
# hist(bipower_df[!bipower_df[, 'outcomes'], 'scores'],
#      breaks = 1000000,
#      main = 'Histogram of Scores',
#      xlab = 'Scores',
#      col='blue', 
#      xlim = c(0, quantile(bipower_df[, 'scores'], 0.25)))
# 
# hist(bipower_df[bipower_df[, 'outcomes'], 'scores'], 
#      main = 'Histogram of Scores', 
#      xlab = 'Scores', 
#      col='blue', xlim = c(0, quantile(bipower_df[, 'scores'], 0.99)))
# hist(bipower_df[!bipower_df[, 'outcomes'], 'scores'], 
#      col='red', add = T)
# Completely uninformative. 

################################################################################
# Simulate distribution of AUROC Statistics
################################################################################


# Simulate distribution of AUROC Statistics
auroc_table <- data.frame(A_hat = rep(NA, num_sims), 
                          A_hat_plug = rep(NA, num_sims), 
                          x_min_hat = rep(NA, num_sims),
                          y_min_hat = rep(NA, num_sims),
                          alpha_x_hat = rep(NA, num_sims), 
                          gamma_y_hat = rep(NA, num_sims))

for (sim_num in 1:num_sims) {
  
  if (round(sim_num/num_sims*10) == sim_num/num_sims*10) {
    print(sprintf('Now calculating AUROC number %d of %d.', 
                  sim_num, num_sims))
  }
  
  # Generate the sample.
  bipower_df <- bipower_gen(n_x, n_y, 
                            x_min, alpha_x, 
                            y_min, gamma_y,
                            disc_or_cts)
  
  
  # # Power law estimation: Unknown xmin estimated
  # 
  # # Estimate the power law exponents (positive distribution).
  # m_pl <- displ$new(bipower_df[bipower_df[, 'outcomes'], 'scores'])
  # est_pl <- estimate_xmin(m_pl)
  # # est_pl$setXmin(est_pl)
  # gamma_y_hat <- est_pl$pars
  # y_min_hat <- est_pl$xmin
  # # Estimate the power law exponents (negative distribution).
  # m_pl <- displ$new(bipower_df[!bipower_df[, 'outcomes'], 'scores'])
  # est_pl <- estimate_xmin(m_pl)
  # # est_pl$setXmin(est_pl)
  # alpha_x_hat <- est_pl$pars
  # x_min_hat <- est_pl$xmin
  
  # Power law estimation: known xmin imposed on MLE of alpha_x
  
  # Estimate the power law exponents (positive distribution).
  m_pl <- displ$new(bipower_df[bipower_df[, 'outcomes'], 'scores'])
  m_pl$setXmin(y_min)
  est_pl <- estimate_pars(m_pl)
  gamma_y_hat <- est_pl$pars
  y_min_hat <- y_min # Plug in true value
  
  # Estimate the power law exponents (negative distribution).
  m_pl <- displ$new(bipower_df[!bipower_df[, 'outcomes'], 'scores'])
  m_pl$setXmin(x_min)
  est_pl <- estimate_pars(m_pl)
  alpha_x_hat <- est_pl$pars
  x_min_hat <- x_min # Plug in true value
  
  
  # Calculate the theoretical AUROC with a plug-in estimator. 
  # Use known true values of location parameter. 
  A_hat_plug <- auc_bi_power(x_min_hat, alpha_x_hat, 
                             y_min_hat, gamma_y_hat)
  
  
  # Calculate and store empirical AUROC statistic.
  for (tie_method_num in 1:length(ties_method_list)) {
    
    
    # Calculate the empirical AUROC statistic.
    A_hat_rank <- auc_rank(scores = bipower_df[, 'scores'], 
                           outcomes = bipower_df[, 'outcomes'], 
                           ties_method = ties_method_list[tie_method_num])
    
    # Store results from AUROC with ranking options.
    # auroc_table[sim_num, 'A_hat'] <- A_hat_rank
    A_hat_col_name <- sprintf('A_hat_%s', ties_method_list[tie_method_num])
    auroc_table[sim_num, A_hat_col_name] <- A_hat_rank
  
  }
  
  # Store default A_hat calculation using average rank for ties.
  auroc_table[sim_num, 'A_hat'] <- auroc_table[sim_num, 'A_hat_average']
  
  
  
  # Store results from known power law estimates. 
  auroc_table[sim_num, 'A_hat_plug'] <- A_hat_plug
  auroc_table[sim_num, 'x_min_hat'] <- x_min_hat
  auroc_table[sim_num, 'y_min_hat'] <- y_min_hat
  auroc_table[sim_num, 'alpha_x_hat'] <- alpha_x_hat
  auroc_table[sim_num, 'gamma_y_hat'] <- gamma_y_hat
  
}

# Select a particular A_hat for plotting.
tie_method_num <- 1
A_hat_col_name <- sprintf('A_hat_%s', ties_method_list[tie_method_num])
auroc_table[, 'A_hat'] <- auroc_table[, A_hat_col_name]


# Plot comparison of Empirical and plug-in estimators.
min_A_hat <- min(auroc_table[, c('A_hat', 'A_hat_plug')])
max_A_hat <- max(auroc_table[, c('A_hat', 'A_hat_plug')])

row_sel <- auroc_table[, 'x_min_hat'] == x_min & 
  auroc_table[, 'y_min_hat'] == y_min



plot(auroc_table[row_sel, 'A_hat'], 
     auroc_table[row_sel, 'A_hat_plug'], 
     col = 'blue', 
     main = 'AUROC Estimators', 
     xlab = 'Empirical AUROC',
     ylab = 'Plug-in Estimator',
     xlim = c(min_A_hat, max_A_hat), 
     ylim = c(min_A_hat, max_A_hat))
lines(c(min_A_hat, max_A_hat), c(min_A_hat, max_A_hat), 
      col = 'black', lwd = 3)
lines(c(A_hat_bi_power, A_hat_bi_power), c(min_A_hat, max_A_hat), 
      col = 'black', lwd = 2)
lines(c(min_A_hat, max_A_hat), c(A_hat_bi_power, A_hat_bi_power), 
      col = 'black', lwd = 2)


hist(auroc_table[, 'A_hat_plug'], 
     breaks = 100,
     main = 'Histogram of AUROC',
     xlab = 'AUROC',
     col='red',
     xlim = c(min_A_hat, max_A_hat))
hist(auroc_table[, 'A_hat'], 
     breaks = 100,
     col='blue',
     add = TRUE)
lines(c(A_hat_bi_power, A_hat_bi_power), c(0, num_sims), 
      col = 'black', lwd = 3)


summary(auroc_table)
A_hat_bi_power


################################################################################
# End
################################################################################
