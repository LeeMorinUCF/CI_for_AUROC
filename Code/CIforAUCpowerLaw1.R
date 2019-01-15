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

n_x_list <- c(100, 500, 1000)
n_y_list <- c(100, 500, 1000)
y_min_list <- c(1.0, 1.25, 1.5, 1.75, 2.0)
x_min_list <- c(1.0, 1.25, 1.5, 1.75, 2.0)
gamma_y_list <- c(1.0, 1.25, 1.5, 1.75, 2.0)


# Set parameters for positive score distribution. 
y_min <- 2.5
gamma_y <- 1.2

# Set parameters for negative score distribution. 
x_min <- 2
alpha_x <- 1.5




################################################################################
# End
################################################################################
