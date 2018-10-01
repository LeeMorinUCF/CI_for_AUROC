################################################################################
# 
# Calculation of COnfidence Intervals for the Area Under the ROC Curve
# 
# Lee Morin, Ph.D.
# Adjunct Assistant Professor
# Queen's University
# 
# December 2, 2017
# 
################################################################################
# 
# This program calculates a series of confidence intervals for the area under
# the ROC curve, comparing a null version with a number of existing methods.
# It serves as a benchmark aginst which to compar new methods.
# 
# This version is modified to include data.table versions of the functions.
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
wd_path <- 'C:/Users/iky155/Documents/Fisher/CIforAUC' 
setwd(wd_path)

# Set path to data folder.
# doc_path <- sprintf('%s', wd_path)

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
null_file_name <- 'CIforAUROCnull1.R'
null_path_file_name <- sprintf('%s/%s', fun_path, null_file_name)
source(null_path_file_name)


################################################################################
# Set Parameters and Generate Data
################################################################################

#--------------------------------------------------------------------------------
# Binormal model.
#--------------------------------------------------------------------------------

# Set parameters.
mu_x <- 1
sigma_x <- 1
mu_y <- 2
sigma_y <- 1

# Typical unbalanced case (requires re-weighting).
# n_x <- 1000
# n_y <- 100
# Balanced case, to start.
n_x <- 500
n_y <- 50


# Binormal Model.
binorm_dt <- data.table(binorm_gen(n_x, n_y, mu_x, mu_y, sigma_x, sigma_y))

# Trust but verify.
summary(binorm_dt)
summary(binorm_dt[outcomes == TRUE, ])
summary(binorm_dt[outcomes == FALSE, ])



################################################################################
# Optimize Distribution
################################################################################


#--------------------------------------------------------------------------------
# Measure AUROC
#--------------------------------------------------------------------------------


# True value.
auc_binorm_true <- auc_binorm(mu_x, sigma_x, mu_y, sigma_y)

# Known true model.
auc_binorm_binorm <- auc_binorm(mu_0 = mean(binorm_dt[outcomes == FALSE, scores]), 
                                sigma2_0 = sd(binorm_dt[outcomes == FALSE, scores]), 
                                mu_1 = mean(binorm_dt[outcomes == TRUE, scores]), 
                                sigma2_1 = sd(binorm_dt[outcomes == TRUE, scores]))

# Distribution free estimate.
binorm_dt[, rank_obs := rank(scores)]
A_hat_curr_binorm <- binorm_dt[outcomes == TRUE, ( sum(rank_obs) - 0.5*n_y*(n_y + 1) )/(n_y*n_x)]

# Weighted version (with initially equal weights).
# A_hat_current_dt <- auc_weighted_loops_dt(binorm_dt)
A_hat_current_dt <- auc_weighted_loops(scores = binorm_dt[, scores], 
                                       outcomes = binorm_dt[, outcomes], 
                                       weights = binorm_dt[, weights])

# Compare results for consistency.
auc_binorm_true
auc_binorm_binorm
A_hat_curr_binorm
A_hat_current_dt


#--------------------------------------------------------------------------------
# Optimize distribution to match AUROC
#--------------------------------------------------------------------------------

# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
A_0 <- 0.85
tol <- 0.00000001
max_iter <- 20
ylim <- c(0, 0.4)

# Data frame version.
new_weights <- adj_weights(scores = binorm_dt[, scores], 
                           outcomes = binorm_dt[, outcomes], 
                           weights = binorm_dt[, weights], 
                           eta, A_0, A_hat = A_hat_curr_binom, tol, max_iter, 
                           display = TRUE, ylim = ylim)

# Data table version.
adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
               display = TRUE, ylim = ylim)

summary(binorm_dt)


#--------------------------------------------------------------------------------
# Test summing weights by quantile.
# For distance calculations afte optimization.
#--------------------------------------------------------------------------------

# Calculate distribution for actual observations (uniform).
num_quantiles <- 10
quantile_df <- distn_quantiles(binorm_dt, num_quantiles)
binorm_dt[outcomes == TRUE , 
          quantile := cut(scores,
                          breaks = quantile_df[, 'pos_quantile'],
                          include.lowest = TRUE, labels = 1:num_quantiles) ]
table(binorm_dt[, quantile], useNA = 'ifany')
binorm_dt[outcomes == FALSE , 
          quantile := cut(scores,
                          breaks = quantile_df[, 'neg_quantile'],
                          include.lowest = TRUE, labels = 1:num_quantiles) ]
table(binorm_dt[, quantile], useNA = 'ifany')

# Replicate the table of distributions.
distn_table <- data.frame(quantiles = c(sprintf('Q_1_%d', 1:num_quantiles), 
                                        sprintf('Q_0_%d', 1:num_quantiles)), 
                          outcomes = c(rep(TRUE, num_quantiles), 
                                       rep(FALSE, num_quantiles)), 
                          distn_obs = c(as.numeric(table(binorm_dt[outcomes == TRUE, 
                                                                   quantile]) / n_y)[1:num_quantiles],
                                        as.numeric(table(binorm_dt[outcomes == FALSE, 
                                                                   quantile]) / n_x)[1:num_quantiles]))

# Calculate for optimized weights.
distn_table[distn_table[, 'outcomes'] == TRUE, 'distn_sim'] <- 
  binorm_dt[outcomes == TRUE, sum(opt_weights), by = quantile][, V1]
distn_table[distn_table[, 'outcomes'] == FALSE, 'distn_sim'] <- 
  binorm_dt[outcomes == FALSE, sum(opt_weights), by = quantile][, V1]



test <- binorm_dt[outcomes == TRUE, sum(opt_weights), by = quantile]
colnames(test)
test[, V1]

sum(test[, V1])

################################################################################
# Run Simulations
################################################################################

num_boots <- 100
num_quantiles <- 10

auroc_sim_df <- auroc_sim(binorm_dt, metric = 'kld0_joint', num_quantiles, num_boots)

# auroc_sim_df*100
head(auroc_sim_df)
summary(auroc_sim_df)



#--------------------------------------------------------------------------------
# Analyze Simulation Results
#--------------------------------------------------------------------------------

plot(auroc_sim_df[, 'auroc'], auroc_sim_df[, 'distance'])
hist(auroc_sim_df[, 'distance'])
hist(auroc_sim_df[, 'auroc'])

# lm_fit <- lm(formula = distance ~ auroc, 
#              data = auroc_sim_df[!is.infinite(auroc_sim_df[, 'distance']), ])
lm_fit <- lm(formula = distance ~ auroc, data = auroc_sim_df)
summary(lm_fit)


#--------------------------------------------------------------------------------
# Calculate Confidence Intervals
#--------------------------------------------------------------------------------



################################################################################
# End
################################################################################

