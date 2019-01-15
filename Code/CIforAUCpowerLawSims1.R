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




################################################################################
# Simulation Analysis 
# First: Distributions of AUROC Statistics
# Later: Coverage Rates
################################################################################


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
