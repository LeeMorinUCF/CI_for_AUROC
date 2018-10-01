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

# Set path to folder for figures.
fig_path <- sprintf('%s/Figures', wd_path)

# Load data.table package for efficient use of databases.
library(data.table)

# Load competing package for ROC analysis.
library(pROC)

# Load package for weighted histograms.
library(plotrix)

# Load library for series of colors.
library(colorspace)


# Load required libraries.
bench_file_name <- 'CIforAUROCbench1.R'
bench_path_file_name <- sprintf('%s/%s', fun_path, bench_file_name)
source(bench_path_file_name)
null_file_name <- 'CIforAUROCnull1.R'
null_path_file_name <- sprintf('%s/%s', fun_path, null_file_name)
source(null_path_file_name)


################################################################################
# Set Parameters
################################################################################




################################################################################
# Run Simulations
################################################################################


auroc_vec <- c(0.75)
n_x <- 1000
n_y <- 100

num_sims <- 100

num_quantiles <- 10
num_boots <- 200
alpha <- 0.05

metric <- 'kld0_joint'
eta <- 4
tol <- 0.00000001
max_iter <- 20


binorm_coverage <- binorm_coverage_sim(n_x, n_y, auroc_vec, num_sims,
                            metric, eta, tol, max_iter, 
                            num_quantiles, num_boots, alpha)
binorm_coverage


################################################################################
# Testing functions
################################################################################


# Create a sample data table to test.

# Binormal Model.
binorm_spec_list <- binorm_spec(auroc_vec)
auroc_num <- 1
A_0 <- auroc_vec[auroc_num]
mu_x = binorm_spec_list[auroc_num, 'mu_x']
sigma_x = binorm_spec_list[auroc_num, 'sigma_x']
mu_y = binorm_spec_list[auroc_num, 'mu_y']
sigma_y = binorm_spec_list[auroc_num, 'sigma_y']
binorm_dt <- data.table(binorm_gen(n_x, n_y, mu_x, mu_y, sigma_x, sigma_y))


# Test the confidence intervals inside.
binorm_ci_df <- binorm_ci_calc(binorm_dt, 
                               mu_x, sigma_x, mu_y, sigma_y,
                               metric, eta, A_0, tol, max_iter, 
                               num_quantiles, num_boots, alpha)
binorm_ci_df


#--------------------------------------------------------------------------------
# Testing alternative AUROC calculations.
#--------------------------------------------------------------------------------


# Test fast AUROC variance.
var_auroc_fast <- var_auroc_fast_calc(binorm_dt)
nrow(binorm_dt)

# Compare with pROC package.
auc(binorm_dt[, outcomes], binorm_dt[, scores])


# Compare with rank calculation.
n_y <- binorm_dt[outcomes == TRUE, .N]
n_x <- binorm_dt[outcomes == FALSE, .N]
binorm_dt[, rank_z := rank(scores)]
auroc_hat_est <- binorm_dt[outcomes == TRUE, ( sum(rank_z) - 0.5*n_y*(n_y + 1) )/(n_y*n_x)]
# Checks out.

#--------------------------------------------------------------------------------
# Testing pROC package.
#--------------------------------------------------------------------------------


# Consider standard errors.

rocobj <- roc(binorm_dt[, outcomes], binorm_dt[, scores])

var_auroc_delong <- var(roc = rocobj, method = "delong",
                        boot.n = num_boots, boot.stratified = TRUE)
var_auroc_boot <- var(roc = rocobj, method = "bootstrap",
                      boot.n = num_boots, boot.stratified = TRUE)
var_auroc_obuch <- var(roc = rocobj, method = "obuchowski",
                       boot.n = num_boots, boot.stratified = TRUE)

auroc_proc <- auc(binorm_dt[, outcomes], binorm_dt[, scores])

auc_ci_delong <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "delong", 
                        boot.n = num_boots, boot.stratified = TRUE)
auc_ci_boot <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "bootstrap", 
                      boot.n = num_boots, boot.stratified = TRUE)

c(auroc_proc - qnorm(1 - alpha/2)*sqrt(var_auroc_delong), 
  auroc_proc + qnorm(1 - alpha/2)*sqrt(var_auroc_delong))
c(auroc_proc - qnorm(1 - alpha/2)*sqrt(var_auroc_boot), 
  auroc_proc + qnorm(1 - alpha/2)*sqrt(var_auroc_boot))

# cu_l and cu_u are in elements 1 and 3.
auc_ci_delong[1]
auc_ci_delong[2]
auc_ci_delong[3]

auc_ci_boot[1]
auc_ci_boot[2]
auc_ci_boot[3]


#--------------------------------------------------------------------------------
# Checking simulated auroc and distance.
#--------------------------------------------------------------------------------


colnames(binorm_dt)
summary(binorm_dt)

A_hat_new_binom_dt <- auc_weighted_loops(scores = binorm_dt[, scores], 
                                         outcomes = binorm_dt[, outcomes], 
                                         weights = binorm_dt[, opt_weights])


auroc_sim_df <- auroc_sim(binorm_dt, metric, num_quantiles, num_boots)
summary(auroc_sim_df)
nrow(auroc_sim_df)
head(auroc_sim_df)


#--------------------------------------------------------------------------------
# Solving for the limiting A_0 a set distance away.
#--------------------------------------------------------------------------------


# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
# A_0 <- 0.85
tol <- 0.00000001
max_iter <- 20
ylim <- c(0, 0.4)

adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
               display = FALSE, ylim = NA)

A_hat_new_binom_dt <- auc_weighted_loops(scores = binorm_dt[, scores], 
                                         outcomes = binorm_dt[, outcomes], 
                                         weights = binorm_dt[, opt_weights])


summary(auroc_sim_df)
max(auroc_sim_df[, 'auroc'])
quantile(auroc_sim_df[, 'distance'], 1)
quantile(auroc_sim_df[, 'distance'], 1 - alpha)

# distance_A_lim <- dist_A_lim_calc(dt_in, metric, min_eps = 10^(-4), 
#                                   eta, A_lim, tol_A_0, max_iter_A_0, 
#                                   display = FALSE, ylim = NA)


# Solve for upper root.
# uniroot(f, interval, lower = min(interval), upper = max(interval),
#         tol = .Machine$double.eps^0.25, maxiter = 1000, ...)



# uniroot(f = dist_A_lim_calc, 
#         interval = c(0, max(auroc_sim_df[, 'auroc'])), 
#         tol = .Machine$double.eps^0.25, maxiter = 1000, 
#         dt_in = binorm_dt, metric, min_eps = 10^(-4), 
#         eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
#         display = FALSE, ylim = NA)
# Fail: Need parameter of interest first.

# Solve for root of distance difference as function of AUROC limit. 
D_bar <- quantile(auroc_sim_df[, 'distance'], 1 - alpha)

# Test the function first.
A_lim_test <- 0.76
dist_diff <- dist_A_lim_diff(A_lim = A_lim_test, 
                             D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                             eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                             display = FALSE, ylim = NA)


# Test the end points.
dist_diff <- dist_A_lim_diff(A_lim = auroc_hat_est, 
                             D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                             eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                             display = FALSE, ylim = NA)
dist_diff
dist_diff <- dist_A_lim_diff(A_lim = 1.25*max(auroc_sim_df[, 'auroc']), 
                             D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                             eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                             display = FALSE, ylim = NA)
dist_diff
dist_diff <- dist_A_lim_diff(A_lim = 0.0*min(auroc_sim_df[, 'auroc']), 
                             D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                             eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                             display = FALSE, ylim = NA)
dist_diff



# Solve for the upper root.
ci_u_solve <- uniroot(f = dist_A_lim_diff, 
                      # interval = c(auroc_hat_est, max(auroc_sim_df[, 'auroc'])), 
                      interval = c(auroc_hat_est, 1), 
                      tol = .Machine$double.eps^0.25, maxiter = 100, 
                      D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                      eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                      display = FALSE, ylim = NA)

# Solve for the lower root.
ci_l_solve <- uniroot(f = dist_A_lim_diff, 
                      # interval = c(min(auroc_sim_df[, 'auroc']), auroc_hat_est), 
                      interval = c(0, auroc_hat_est), 
                      tol = .Machine$double.eps^0.25, maxiter = 100, 
                      D_bar, dt_in = binorm_dt, metric, num_quantiles, min_eps = 10^(-4), 
                      eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                      display = FALSE, ylim = NA)
auc_ci_root <- c(ci_l_solve$root, ci_u_solve$root)
auc_ci_root



# Finding an initial interval.
summary(auroc_sim_df)
quantile(auroc_sim_df[, 'distance'], 1 - alpha)


max(auroc_sim_df[, 'auroc'])

D_bar <- quantile(auroc_sim_df[, 'distance'], 1 - alpha)

summary(auroc_sim_df)
summary(auroc_sim_df[auroc_sim_df[, 'distance'] < D_bar, ])


max(auroc_sim_df[auroc_sim_df[, 'distance'] < D_bar, 'auroc'])


plot(auroc_sim_df[, 'auroc'], auroc_sim_df[, 'distance'])



################################################################################
# Power Analysis
################################################################################

auroc_vec <- c(0.65, 0.75, 0.85)
auroc_grid <- seq(-0.10, 0.10, by = 0.02)
stats_list <- c('binorm_true', 'binorm_est', 'upper', # 'distn_free', 
                'auroc_delong', 'auroc_boot', 'auroc_root',
                'dist_sim', 'auroc_sim', 'dist_sim_adj', 'auroc_sim_adj')
binorm_power <- expand.grid(auroc_grid, auroc_vec, stats_list)

binorm_power <- data.frame(A_0 = binorm_power[, 'Var2'],
                           A_true = binorm_power[, 'Var2'] + binorm_power[, 'Var1'],
                           stat = binorm_power[, 'Var3'],
                           power = rep(NA,nrow(binorm_power)))

head(binorm_power, 100)


num_sims <- 2
binorm_power <- binorm_power_sim(n_x, n_y, auroc_vec, auroc_grid, num_sims,
                                 metric, eta, tol, max_iter, 
                                 num_quantiles, num_boots, alpha)


# Save the result.

# Fill it with a set of dummy power curves.
binorm_power[, 'power'] <- (binorm_power[, 'A_true'] - binorm_power[, 'A_0'])^2
binorm_power[, 'power'] <- binorm_power[, 'power'] / (0.1 + binorm_power[, 'power']) + 
  0.05*binorm_power[, 'A_true'] + 
  rep(seq(1, length(stats_list))/25, rep(length(auroc_vec)*length(auroc_grid), length(stats_list)))
summary(binorm_power)

# Create separate figures for each A_0 value.
# For those, put power curves for each statistic.
# Set color palette for each statistic.
color_list <- rainbow_hcl(length(stats_list))
line_wd <- 2


# Plot and save the graphs.
num_null_auroc <- length(auroc_vec)
# null_auroc_num <- 1
for (null_auroc_num in 1:num_null_auroc) {
  
  A_0 <- auroc_vec[null_auroc_num] 
  A_true_grid <- A_0 + auroc_grid
  
  
  x_min <- min(A_true_grid)
  x_max <- max(A_true_grid)
  
  # Set the path and file name for each figure.
  fig_file_name <- sprintf('/Power/fig_power_A_0_%d.png', round(A_0*100))
  fig_path_file_name <- sprintf('%s/%s', fig_path, fig_file_name)
  
  png(fig_path_file_name)
  plot(NA,
       main = sprintf('Power Curves for A_0 = %3.2f', A_0),
       xlab = 'True Value of the AUROC',
       ylab = 'Power',
       ylim = c(0, 1),
       xlim = c(x_min, x_max))
  
  # Add power curves one at a time.
  # stat_num <- 1
  for (stat_num in 1:length(stats_list)) {
    
    
    
    power_curve <- binorm_power[binorm_power[, 'A_0'] == A_0 & 
                                  binorm_power[, 'stat'] == stats_list[stat_num], 'power']
    lines(A_true_grid, power_curve, col = color_list[stat_num], lwd = line_wd)
    abline(h = 0.0, col = 'black', lwd = 1)
    abline(h = 0.05, col = 'black', lwd = 3)
    
  }
  # Close the figure for this A_0.
  dev.off()
  
}


rainbow_hcl(3)
library("colorspace")
pal <- choose_palette()


################################################################################
# Variability of AUROC
################################################################################

# Plot range of AUROC as a function of distance from true value.
# Plot for a particular dataset.
# See that AUROC values a particular distance from the true value 
# have more variation than the distance itself.
# Should show a sideways-hourglass shape that maximizes range at zero and one.
# Plot for multiple sample sizes.



################################################################################
# End
################################################################################
