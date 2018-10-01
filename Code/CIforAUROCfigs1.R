################################################################################
# 
# Calculation of Confidence Intervals for the Area Under the ROC Curve
# 
# Lee Morin, Ph.D.
# Adjunct Assistant Professor
# Queen's University
# 
# December 2, 2017
# 
################################################################################
# 
# This program creates figures for documents involving confidence intervals 
# for the area under the ROC curve.
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
# wd_path <- 'C:/Users/iky155/Documents/GitRepositories/CIforAUROC_copy' 
wd_path <- 'C:/Users/iky155/Documents/Repos/CIforAUROC_copy' 
setwd(wd_path)

# Set path to figure folder.
fig_path <- sprintf('%s/Slides/Figs', wd_path)

# Set path to data folder.
# fun_path <- sprintf('%s/function', wd_path)
fun_path <- sprintf('%s', wd_path)

# Load data.table package for efficient use of databases.
library(data.table)

# Load package for weighted histograms.
library(plotrix)

# Load package for calculations related to the ROC curve.
library(pROC)

# Load required libraries.
bench_file_name <- 'CIforAUROCbench1.R'
bench_path_file_name <- sprintf('%s/%s', fun_path, bench_file_name)
source(bench_path_file_name)
null_file_name <- 'CIforAUROCnull2.R'
null_path_file_name <- sprintf('%s/%s', fun_path, null_file_name)
source(null_path_file_name)




################################################################################
# Binormal KLD Calculation
################################################################################

# Set parameters.
mu_x <- 1
sigma_x <- 1
mu_y <- 2
sigma_y <- 1

# Set grid for calculation of density.
x_grid <- seq(-2, 5, by = 0.1)

f_x <- dnorm(x = x_grid, mean = mu_x, sd = sigma_x, log = FALSE)
f_y <- dnorm(x = x_grid, mean = mu_y, sd = sigma_y, log = FALSE)
g_x <- dnorm(x = x_grid, mean = mu_x, sd = sigma_x, log = TRUE)
g_y <- dnorm(x = x_grid, mean = mu_y, sd = sigma_y, log = TRUE)

# Set file name for plot.
kld_file_name <- 'KLD_calc_1.png'
fig_path_file_name <- sprintf('%s/%s', fig_path, kld_file_name)


# Plot grid with only first term.
png(fig_path_file_name)
plot(x_grid, 
     f_y, 
     type = 'n',
     # main = 'AR(1) Process', 
     xlab = 'Classification variable', 
     ylab = 'Density, Distance',
     xlim = c(-0.5, 3.5),
     ylim = c(-0.3, 0.45),
     cex.lab = 1.5
)
lines(x_grid, 
      0*f_y,
      type = 'l', 
      lwd = 2,
      col = 'black')
lines(x_grid, 
      f_y,
      type = 'l', 
      lwd = 5,
      col = 'red')
lines(x_grid, 
      f_x,
      type = 'l', 
      lwd = 5,
      col = 'blue')
lines(x_grid, 
      f_x - f_y,
      type = 'l', 
      lwd = 5,
      col = 'magenta')
lines(x_grid, 
      g_x - g_y,
      type = 'l', 
      lwd = 5,
      col = 'orange')
dev.off()



# Plot grid with scale large enough for log term.
# Not informative

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
# Figure for Optimizing Distributions
################################################################################


#--------------------------------------------------------------------------------
# Find nearest distribution.
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Plot figure.
#--------------------------------------------------------------------------------



################################################################################
# Figures for Power Curves
################################################################################


#--------------------------------------------------------------------------------
# Obtain power statistics.
#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Plot figure.
#--------------------------------------------------------------------------------



################################################################################
# Figures for Confidence Interval Widths
################################################################################

# For a set AUROC, display confidence interval as a function of tuning parameter.

#--------------------------------------------------------------------------------
# Calculate matrices of confidence intervals.
#--------------------------------------------------------------------------------

# Set parameters and build dataset.
alpha <- 0.05
n_x <- 500
n_y <- 50
auroc <- 0.7
binorm_spec_list <- binorm_spec(auroc)

mu_x <- binorm_spec_list$mu_x
mu_y <- binorm_spec_list$mu_y
sigma_x <- binorm_spec_list$sigma_x
sigma_y <- binorm_spec_list$sigma_y


binorm_dt <- data.table(binorm_gen(n_x, n_y, 
                                   mu_x = binorm_spec_list$mu_x, 
                                   mu_y = binorm_spec_list$mu_y, 
                                   sigma_x = binorm_spec_list$sigma_x, 
                                   sigma_y = binorm_spec_list$sigma_y))





# Calculate AUROC
A_hat <- auc_weighted_loops(scores = binorm_dt[, scores], 
                            outcomes = binorm_dt[, outcomes], 
                            weights = binorm_dt[, weights])


# Calculate a variety of confidence intervals.
auc_ci_binorm_true <- pnorm(confidence_interval(mean = (mu_y - mu_x) / 
                                                  sqrt(sigma_x^2 + sigma_y^2), 
                                                var = (sigma_x^2/n_x + sigma_y^2/n_y) / 
                                                  (sigma_x^2 + sigma_y^2), 
                                                alpha))

auc_ci_binorm <- auc_ci_binorm(scores = binorm_dt[, scores], 
                               outcomes = binorm_dt[, outcomes], 
                               alpha = alpha)

auc_ci_upper_bound <- auc_ci_upper_bound(auc = A_hat, n_0 = n_x, n_1 = n_y, alpha)

# From package pROC.
num_boots <- 1000
rocobj <- roc(binorm_dt[, outcomes], binorm_dt[, scores])
auc_ci_delong <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "delong", 
                        boot.n = num_boots, boot.stratified = TRUE, progress = 'none')
auc_ci_boot <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "bootstrap", 
                      boot.n = num_boots, boot.stratified = TRUE, progress = 'none')

auc_ci_binorm_true
auc_ci_binorm
auc_ci_delong
auc_ci_boot
auc_ci_upper_bound

#--------------------------------------------------------------------------------
# Calculate confidence intervals with fixed error rate.
#--------------------------------------------------------------------------------


# Fixed Error Rate.
error_rate_list <- seq(0, 0.1, by = 0.05)
ci_err_list <- data.frame(error_rate = numeric(length(error_rate_list)), 
                          ci_u = numeric(length(error_rate_list)), 
                          ci_l = numeric(length(error_rate_list)))
for (error_rate in error_rate_list) {
  
  error_count <- error_rate*(n_x +n_y)
  
  # Calculate CI
  mean_auc_fix_err <- mean_auc_fixed_error(n_0 = n_x, n_1 = n_y, k = error_count, max_n = 1000)
  auc_ci_fix_err <- auc_ci_fixed_error(n_0 = n_x, n_1 = n_y, k = error_count, max_n = 1000, 
                                      mean_auc = mean_auc_fix_err, alpha)
  
  # Store values.
  ci_err_list[ci_err_list[, 'error_rate'] == error_rate, 
              c('ci_l', 'ci_u')] <- auc_ci_fix_err
  
}



# Fixed Distance metric.


# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
tol <- 0.00001
max_iter <- 20
num_quantiles <- 10

# Initialize.
dist_list <- seq(0, 4.0, by = 0.05)
ci_dist_list <- data.frame(dist = dist_list, 
                          ci_u = numeric(length(dist_list)), 
                          ci_l = numeric(length(dist_list)))

for (dist in dist_list) {
  
  # Progress report.
  print(sprintf('Calculating CI for distance %f.', dist))
  
  
  # Calculate CI
  # Solve for the upper root.
  ci_u_solve <- uniroot(f = dist_A_lim_diff, 
                        interval = c(A_hat, 1), 
                        tol = .Machine$double.eps^0.25, maxiter = 100, 
                        D_bar = dist, 
                        dt_in = binorm_dt, 
                        metric = 'kld0_joint', 
                        num_quantiles, min_eps = 10^(-4), 
                        eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                        display = FALSE, ylim = NA)
  
  # Solve for the lower root.
  ci_l_solve <- uniroot(f = dist_A_lim_diff, 
                        interval = c(0, A_hat), 
                        tol = .Machine$double.eps^0.25, maxiter = 100, 
                        D_bar = dist, 
                        dt_in = binorm_dt, 
                        metric = 'kld0_joint', 
                        num_quantiles, min_eps = 10^(-4), 
                        eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                        display = FALSE, ylim = NA)
  
  
  # Store values.
  ci_dist_list[ci_dist_list[, 'dist'] == dist, 
              c('ci_l', 'ci_u')] <- c(ci_l_solve$root, ci_u_solve$root)
  
  print(sprintf('CI for distance %f:', dist))
  print(c(ci_l_solve$root, ci_u_solve$root))
  
}


# Set file name for plot.
fcst_file_name <- 'Forecast_int_1.png'
fig_path_file_name <- sprintf('%s/%s', fig_path, fcst_file_name)



# Plot the confidence bounds.
png(fig_path_file_name)
plot(NA,
     type = 'n',
     main = 'Forecasting and Confidence Intervals', 
     xlab = 'Distance', 
     ylab = 'AUROC',
     xlim = c(min(dist_list), max(dist_list))*(num_quantiles - 1),
     ylim = c(0.2, 1.0),
     cex.lab = 1.5)
# Fill in with a line plot with settings for line style.
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      ci_dist_list[, 'ci_l'],
      type = 'l', 
      lwd = 3,
      col = 'black')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      ci_dist_list[, 'ci_u'],
      type = 'l', 
      lwd = 3,
      col = 'black')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(A_hat, length(dist_list)),
      type = 'l', 
      lwd = 2,
      col = 'black')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(1.0, length(dist_list)),
      type = 'l', 
      lwd = 2,
      col = 'black')
# lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
#       rep(0.5, length(dist_list)),
#       type = 'l', 
#       lwd = 2,
#       col = 'black')

# Add confidence bounds.
auc_ci_binorm_true
auc_ci_binorm
auc_ci_delong
auc_ci_boot
auc_ci_upper_bound

lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_binorm_true[1], length(dist_list)), type = 'l', lwd = 2, col = 'blue')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_binorm_true[2], length(dist_list)), type = 'l', lwd = 2, col = 'blue')

lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_binorm[1], length(dist_list)), type = 'l', lwd = 2, col = 'green')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_binorm[2], length(dist_list)), type = 'l', lwd = 2, col = 'green')

lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_delong[1], length(dist_list)), type = 'l', lwd = 2, col = 'red')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_delong[3], length(dist_list)), type = 'l', lwd = 2, col = 'red')

lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_boot[1], length(dist_list)), type = 'l', lwd = 2, col = 'orange')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_boot[3], length(dist_list)), type = 'l', lwd = 2, col = 'orange')

lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_upper_bound[1], length(dist_list)), type = 'l', lwd = 2, col = 'yellow')
lines(ci_dist_list[, 'dist']*(num_quantiles - 1), 
      rep(auc_ci_upper_bound[2], length(dist_list)), type = 'l', lwd = 2, col = 'yellow')


legend(x = 'bottomleft', y = NULL, 
       legend = c('True model', 'Binormal', 'DeLong, et. al.', 'Bootstrap', 'Upper bound', 'Forecast'), 
       fill = c('blue', 'green', 'red', 'orange', 'yellow', 'black')
       # col = c('blue', 'green', 'red', 'orange', 'yellow')
       )
dev.off()



#--------------------------------------------------------------------------------
# Plot figure.
#--------------------------------------------------------------------------------


################################################################################
# Figures for Confidence Interval Widths
################################################################################

# For a grid of AUROC values, display a variety of confidence intervals.

#--------------------------------------------------------------------------------
# Calculate matrices of confidence intervals.
#--------------------------------------------------------------------------------

auroc_grid <- seq(0.5, 1.0, by = 0.005)


#--------------------------------------------------------------------------------
# Plot figure.
#--------------------------------------------------------------------------------





################################################################################
# End
################################################################################

