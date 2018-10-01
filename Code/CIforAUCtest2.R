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
binorm_df <- binorm_gen(n_x, n_y, mu_x, mu_y, sigma_x, sigma_y)

# Trust but verify.
summary(binorm_df)
summary(binorm_df[binorm_df[, 'outcomes'], ])
summary(binorm_df[!binorm_df[, 'outcomes'], ])


# Set arguments for troubleshooting functions.
scores = binorm_df[, 'scores']
outcomes = binorm_df[, 'outcomes']
weights = binorm_df[, 'weights']


#--------------------------------------------------------------------------------
# Binormal model - Data Table
#--------------------------------------------------------------------------------

# Create a matching data table.
binorm_dt <- data.table(binorm_df)
# Prime optimal weights with observed weights.
binorm_dt[, opt_weights := weights]

colnames(binorm_dt)
summary(binorm_dt)

binorm_dt[, sum(outcomes)]
binorm_dt[, sum(!outcomes)]

head(binorm_dt)
tail(binorm_dt)
binorm_dt[549]
binorm_dt[546]
binorm_dt[549, scores]
which(binorm_dt[, outcomes])
which(binorm_dt[, !outcomes])

# Check quantiles.
quantile_df <- distn_quantiles(binorm_dt, num_quantiles = 10)

summary(binorm_dt[outcomes, scores])
summary(binorm_dt[!outcomes, scores])
quantile_df

# Map scores to quantiled bins.
binorm_dt[outcomes , 
          quantile := cut(scores,
                          breaks = quantile(scores,
                                            probs = seq(0, 1, by = 0.1), na.rm = T),
                          include.lowest = TRUE, labels = 1:10) ]
binorm_dt[!outcomes , 
          quantile := cut(scores,
                          breaks = quantile(scores,
                                            probs = seq(0, 1, by = 0.1), na.rm = T),
                          include.lowest = TRUE, labels = 1:10) ]
# Uniform by definition.
table(binorm_dt[outcomes, quantile])
table(binorm_dt[!outcomes, quantile])
# Percentages.
table(binorm_dt[outcomes, quantile])/binorm_dt[outcomes, .N]
table(binorm_dt[!outcomes, quantile])/binorm_dt[!outcomes, .N]




#--------------------------------------------------------------------------------
# Biexponential model.
#--------------------------------------------------------------------------------

# Set parameters.
lambda_x <- 3
lambda_y <- 1

# Biexponential model.
biexp_df <- biexp_gen(n_x, n_y, lambda_x, lambda_y)


# Trust but verify.
summary(biexp_df[, 'scores'])
summary(biexp_df[biexp_df[, 'outcomes'], 'scores'])
summary(biexp_df[!biexp_df[, 'outcomes'], 'scores'])

ylim = c(0, 0.4)
plot_distn_shift(scores = binorm_df[, 'scores'], 
                 outcomes = binorm_df[, 'outcomes'], 
                 weights = binorm_df[, 'weights'], 
                 ylim = ylim)


################################################################################
# Test Calculation of AUROC
################################################################################


auc_loops_binorm <- auc_loops(scores = binorm_df[, 'scores'], outcomes = binorm_df[, 'outcomes'])

auc_rank_binorm <- auc_rank(scores = binorm_df[, 'scores'], outcomes = binorm_df[, 'outcomes'])

auc_binorm_binorm <- auc_binorm(mu_0 = mu_x, sigma2_0 = sigma_x, mu_1 = mu_y, sigma2_1 = sigma_y)


A_hat_curr_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                       outcomes = binorm_df[, 'outcomes'], 
                                       weights = binorm_df[, 'weights'])


A_hat_current_dt <- auc_weighted_loops_dt(binorm_dt)

# Compare results for consistency.
auc_loops_binorm
auc_rank_binorm
auc_binorm_binorm
A_hat_curr_binom
A_hat_current_dt




################################################################################
# Test Calculation of Null Hypothesis
################################################################################


# Set parameters for re-weighting toward the null hypothesis.
eta <- 0.5
A_0 <- 0.75
tol <- 0.01
max_iter <- 10

neg_gradient <- pos_cond_shift(scores = binorm_df[, 'scores'], 
                               outcomes = binorm_df[, 'outcomes'], 
                               weights = binorm_df[, 'weights'])
pos_gradient <- neg_cond_shift(scores = binorm_df[, 'scores'], 
                               outcomes = binorm_df[, 'outcomes'], 
                               weights = binorm_df[, 'weights'])
next_weights <- iter_weights(scores = binorm_df[, 'scores'], 
                             outcomes = binorm_df[, 'outcomes'], 
                             weights = binorm_df[, 'weights'], eta, A_0, A_hat = A_hat_curr_binom)


#--------------------------------------------------------------------------------
# Check data table versions.
#--------------------------------------------------------------------------------

# Reset data table. 
binorm_dt <- data.table(binorm_df)
# Prime optimal weights with observed weights.
binorm_dt[, opt_weights := weights]

iter_weights_dt(binorm_dt, eta, A_0, A_hat, 
                display = FALSE, ylim = NA)

# colnames(binorm_dt)
# summary(binorm_dt)

summary(next_weights)
summary(binorm_dt[, get('opt_weights')])

summary(next_weights[outcomes])
summary(binorm_dt[outcomes, get('opt_weights')])

summary(next_weights[!outcomes])
summary(binorm_dt[!outcomes, get('opt_weights')])



#--------------------------------------------------------------------------------


# Recalculate AUROC with next weights.
A_hat_next_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                       outcomes = binorm_df[, 'outcomes'], 
                                       weights = next_weights)
A_hat_curr_binom
A_hat_next_binom


# Compare with gradient approach.
sum(pos_gradient)
sum(neg_gradient)

sum(pos_gradient * weights[outcomes])
sum(neg_gradient * weights[!outcomes])
# Not sure what this equality means.



# Calculate distance.
chi2_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = next_weights, metric = 'chi2')
kld0_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = next_weights, metric = 'kld0')
kld1_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = next_weights, metric = 'kld1')

chi2_dist
kld0_dist
kld1_dist



# Set parameters for re-weighting toward the null hypothesis.
eta <- 4
A_0 <- 0.85
tol <- 0.00000001
max_iter <- 20
ylim <- c(0, 0.4)

new_weights <- adj_weights(scores = binorm_df[, 'scores'], 
                           outcomes = binorm_df[, 'outcomes'], 
                           weights = binorm_df[, 'weights'], 
                           eta, A_0, A_hat = A_hat_curr_binom, tol, max_iter, 
                           display = TRUE, ylim = ylim)

#--------------------------------------------------------------------------------
# Check data table versions.
#--------------------------------------------------------------------------------

# Reset inputs.
rm(dt_in)
summary(dt_in)

# Reset the data table.
binorm_dt <- data.table(binorm_df)
# Prime optimal weights with observed weights.
# binorm_dt[, opt_weights := weights]
summary(binorm_dt)

# new_weights_dt <- adj_weights_dt(dt_in, eta, A_0, tol, max_iter, 
#                                  display = TRUE, ylim = ylim)
adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
               display = TRUE, ylim = ylim)


summary(new_weights[outcomes])
summary(new_weights[!outcomes])

summary(binorm_dt[outcomes, opt_weights])
summary(binorm_dt[!outcomes, opt_weights])

A_hat_new_binom_dt <- auc_weighted_loops(scores = binorm_dt[, scores], 
                                         outcomes = binorm_dt[, outcomes], 
                                         weights = binorm_dt[, opt_weights])


# Trust but verify.
A_hat_new_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                      outcomes = binorm_df[, 'outcomes'], 
                                      weights = new_weights)
A_hat_curr_binom
A_hat_new_binom
A_hat_new_binom_dt


#--------------------------------------------------------------------------------






kld0_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = new_weights, metric = 'kld0_joint')
kld0_dist
# weights <- new_weights




# Calculate distance.
chi2_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = new_weights, metric = 'chi2')
kld0_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = new_weights, metric = 'kld0')
kld1_dist <- auroc_dist(outcomes = binorm_df[, 'outcomes'], 
                        weights_1 = binorm_df[, 'weights'], 
                        weights_2 = new_weights, metric = 'kld1')

chi2_dist
kld0_dist
kld1_dist


# Compare them statistically.
summary(binorm_df[binorm_df[, 'outcomes'], 'weights'] - new_weights[binorm_df[, 'outcomes']])
summary(binorm_df[!binorm_df[, 'outcomes'], 'weights'] - new_weights[!binorm_df[, 'outcomes']])
# Note the mean difference is zero, by definition, since probabilities sum to one.
# Bigger changes to the positive outcomes, which have more weight (fewer outcomes).
summary((binorm_df[binorm_df[, 'outcomes'], 'weights'] - new_weights[binorm_df[, 'outcomes']]) / 
          binorm_df[binorm_df[, 'outcomes'], 'weights'])
summary((binorm_df[!binorm_df[, 'outcomes'], 'weights'] - new_weights[!binorm_df[, 'outcomes']]) / 
          binorm_df[!binorm_df[, 'outcomes'], 'weights'])
# Still bigger proportional changes to the positive outcomes.
# Is there a need for optimal weighting?
# Or modify for distance between joint distributions.


source(bench_path_file_name)
source(null_path_file_name)


################################################################################
# End
################################################################################
