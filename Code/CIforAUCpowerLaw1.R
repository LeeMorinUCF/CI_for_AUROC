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


#--------------------------------------------------------------------------------
# Power Law as true model.
#--------------------------------------------------------------------------------


# Typical unbalanced case (requires re-weighting).
n_x <- 1000
n_y <- 100
# n_x <- 500
# n_y <- 50

# Set parameters for positive score distribution. 
y_min <- 3
gamma_y <- 2

# Set parameters for negative score distribution. 
x_min <- 2
alpha_x <- 3


bipower_df <- bipower_gen(n_x, n_y, 
                          x_min, alpha_x, 
                          y_min, gamma_y)


summary(bipower_df)
summary(bipower_df[bipower_df[, 'outcomes'] == 0, ])
summary(bipower_df[bipower_df[, 'outcomes'] == 1, ])

hist(bipower_df[bipower_df[, 'outcomes'] == 0, 'scores'])
hist(bipower_df[bipower_df[, 'outcomes'] == 1, 'scores'])


head(bipower_df, 10)
tail(bipower_df, 10)



################################################################################
# Calculate AUROC Statistics
################################################################################


# Calculate AUROC by loops.
A_hat_loops <- auc_weighted_loops(scores = bipower_df[, 'scores'], 
                              outcomes = bipower_df[, 'outcomes'], 
                              weights = bipower_df[, 'weights'])
A_hat_loops


A_hat_rank <- auc_rank(scores = bipower_df[, 'scores'], 
                       outcomes = bipower_df[, 'outcomes'])
A_hat_rank

# Compare with theoretical calculation.
A_hat_bi_power <- auc_bi_power(x_min, alpha_x, 
                               y_min, gamma_y)
A_hat_bi_power



#--------------------------------------------------------------------------------
# Testing calculations of AUROC statistics
#--------------------------------------------------------------------------------



# Calculate densities on a grid and integrate. 

x_step_1 <- 0.1
y_step_1 <- 0.1
x_grid_1 <- seq(x_min, 10, by = x_step_1)
x_grid_2 <- seq(11, 100, by = 1)
y_grid_1 <- seq(y_min, 10, by = y_step_1)
y_grid_2 <- seq(11, 100, by = 1)
x_grid <- c(x_grid_1, x_grid_2)
y_grid <- c(y_grid_1, y_grid_2)


# Create a dataset with observations weighted by densities. 

bipower_test <- data.frame(outcomes = c(rep(FALSE, length(x_grid)),
                                        rep(TRUE, length(y_grid))),
                         scores = c(x_grid, y_grid))


# Add weights by outcome using midpoint rule.
bipower_test[bipower_test[, 'outcomes'], 'weights'] <- 
  (gamma_y - 1) * y_min^(gamma_y - 1) * 
  c( (y_grid_1 + 0.5*y_step_1)^(-gamma_y)*y_step_1, 
     (y_grid_2 + 0.5)^(-gamma_y) )


bipower_test[!bipower_test[, 'outcomes'], 'weights'] <- 
  (alpha_x - 1) * x_min^(alpha_x - 1) * 
  c( (x_grid_1 + 0.5*x_step_1)^(-alpha_x)*x_step_1, 
     (x_grid_2 + 0.5)^(-alpha_x) )


# Check for accuracy.
head(bipower_test, 10)
tail(bipower_test, 10)


summary(bipower_test[bipower_test[, 'outcomes'], ])
summary(bipower_test[!bipower_test[, 'outcomes'], ])

plot(bipower_test[bipower_test[, 'outcomes'], 'weights'])
plot(bipower_test[!bipower_test[, 'outcomes'], 'weights'])



sum(bipower_test[bipower_test[, 'outcomes'], 'weights'])
sum(bipower_test[!bipower_test[, 'outcomes'], 'weights'])


# Adjust weights on tails of distributions to normalize. 
# pdf_adj <- 1 - sum(bipower_test[bipower_test[, 'outcomes'], 'weights'])
# bipower_test[bipower_test[, 'outcomes'], 'weights'] <- 
#   bipower_test[bipower_test[, 'outcomes'], 'weights'] * 
#   c(rep(1, length(x_grid_1)), rep(pdf_adj, length(x_grid_1)))


# Calculate AUROC by loops.
A_hat_loops_test <- auc_weighted_loops(scores = bipower_test[, 'scores'], 
                                       outcomes = bipower_test[, 'outcomes'], 
                                       weights = bipower_test[, 'weights'])


# Compare with other calculations.
A_hat_loops
A_hat_rank
A_hat_bi_power
A_hat_loops_test


################################################################################
# Simulate distribution of AUROC Statistics
################################################################################


# Typical unbalanced case (requires re-weighting).
# n_x <- 1000
# n_y <- 100
n_x <- 500
n_y <- 500

# Set parameters for positive score distribution. 
y_min <- 2.5
gamma_y <- 1.2

# Set parameters for negative score distribution. 
x_min <- 2
alpha_x <- 1.5


# Calculate theoretical confidence interval.
A_hat_bi_power <- auc_bi_power(x_min, alpha_x, 
                               y_min, gamma_y)
A_hat_bi_power

# Calculate variance. 
bipower_df <- bipower_gen(n_x, n_y, 
                          x_min, alpha_x, 
                          y_min, gamma_y)
A_hat_var <- var_auroc_fast_calc(dt_in = as.data.table(bipower_df))
A_hat_var

# Calculate confidence interval by inverting z-statistic. 
A_hat_CI <- confidence_interval(mean = A_hat_bi_power, 
                                var = A_hat_var, 
                                alpha = 0.05)
A_hat_CI

# Simulate distribution of AUROC Statistics
num_sims <- 500

auroc_table <- rep(NA, num_sims)

for (sim_num in 1:num_sims) {
  
  if (round(sim_num/num_sims*10) == sim_num/num_sims*10) {
    print(sprintf('Now calculating AUROC number %d of %d.', 
                  sim_num, num_sims))
  }
  
  bipower_df <- bipower_gen(n_x, n_y, 
                            x_min, alpha_x, 
                            y_min, gamma_y)
  
  
  A_hat_rank <- auc_rank(scores = bipower_df[, 'scores'], 
                         outcomes = bipower_df[, 'outcomes'])
  
  auroc_table[sim_num] <- A_hat_rank
  
}

hist(auroc_table)
summary(auroc_table)



################################################################################
# Simulation Analysis of Coverage Rates
################################################################################

# Construct a table with the parameters required. 

# Equal sample sizes is most efficient.
n_y_list <- c(100, 500, 1000)
n_x_list <- c(100, 500, 1000)
# Minimum values are relative to each other.
y_min_list <- c(1.0, 1.5, 2.0)
x_min_list <- c(1.0, 1.5, 2.0)
# Positive observations have longer right tails.
gamma_y_list <- c(1.25, 1.5, 1.75, 2.0)
alpha_x_list <- c(1.5, 2.0, 2.5)
# Still want the tails thick to exaggerate asymmetry. 


# initialize storage of of AUROC statistics
num_sims <- 10

auroc_table <- expand.grid(sim_num = 1:num_sims,
                           n_y = n_y_list, 
                           # n_x = n_x_list, # Impose symmetric sample
                           y_min = y_min_list, 
                           x_min = x_min_list,
                           gamma_y = gamma_y_list, 
                           alpha_x = alpha_x_list)
# Impose symmetric sample sizes.
auroc_table[, 'n_x'] <- auroc_table[, 'n_y']
# Rearrange columns. 
auroc_table <- auroc_table[, c('n_y', 'n_x', 
                               'y_min', 'x_min',
                               'gamma_y', 'alpha_x')]

head(auroc_table, 20)

# Replicate table for required number of replications. 


for (sim_num in 1:num_sims) {
  
  # Print progress report. 
  if (round(sim_num/num_sims*10) == sim_num/num_sims*10) {
    print(sprintf('Now performing replication number %d of %d.', 
                  sim_num, num_sims))
  }
  
  # Extract parameters for this replication. 
  n_x <- auroc_table
  
  bipower_df <- bipower_gen(n_x, n_y, 
                            x_min, alpha_x, 
                            y_min, gamma_y)
  
  
  A_hat_rank <- auc_rank(scores = bipower_df[, 'scores'], 
                         outcomes = bipower_df[, 'outcomes'])
  
  auroc_table[sim_num] <- A_hat_rank
  
  
  
  # Save the table 
  
}



################################################################################
# End
################################################################################
