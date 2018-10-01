################################################################################
# 
# Calculation of Critical Values for Kullback-Leibler Divergence
# 
# Lee Morin, Ph.D.
# Senior Data Scientist
# Capital One Bank Canada
# 
# October 12, 2017
# 
################################################################################
# 
# This program calculates critical values for the Kullback-Leibler divergence
# criterion, for a few examples of benchmark distributions.
# 
# It shows that the critical values coincide with those of the chi-squared
# distribution, with degrees of freedom as the number of probability bins, 
# less one, for the restriction that the probabilities sum to zero.
# 
# Note:
# This result can also be derived theoretically, by using the central limit theorem,
# the continuous mapping theorem and by taking a second order Taylor expansion of 
# the log function, in which the higher order terms vanish as the sample size increases.
# 
################################################################################



################################################################################
# Clearing Workspace and Declaring Packages
################################################################################

# Clear workspace.
rm(list=ls(all=TRUE))

# No packages required.


################################################################################
# Set parameters for benchmark distributions and samples.
################################################################################

# Set sample sizes.
num_bench <- 10000
num_current <- 1000

# Set distributions.
# unif_distn <- rep(0.1,10)
unif_sample <- rep(1000,10)
unif_distn <- unif_sample/sum(unif_sample)


# 'Linearly' skewed distribution.
lin_sample <- seq(7500, 12000, by = 500)
lin_distn <- lin_sample/sum(lin_sample)
lin_distn

# 'Exponentially' skewed distribution.
exp_sample <- c(rep(2000,7),3000,5000,10000)
exp_distn <- exp_sample/sum(exp_sample)
exp_distn

# Collect the examples into a list.
example_distn <- list(unif_distn = unif_distn,
                      lin_distn = lin_distn,
                      exp_distn = exp_distn)
example_distn_label <- c('uniform', 'linearly skewed', 'highly skewed')


################################################################################
# Functions.
################################################################################


# Define a function for calculating the distance between 
# the benchmark and current distributions.
# 
# Inputs: 
#   curr_distn is a numeric vector of probabilities in the curent distribution
#     that is being compared to the benchmark.
#   bench_distn is a numeric vector of probabilities in the curent distribution.
#   scale_factor is a numeric scalar to scale the KLD/PSI statistic. The default 
#     of 1 implies no normalization, as it is commonly computed at Capital One.
#     When this is set to the current sample size, the statistic takes on a 
#     chi-squared distribution, as described above. 
# 
KL_divergence <- function(curr_distn, bench_distn, scale_factor = 1) {
  
  # Calculate the statistic itself, with scaling as desired.
  sum((curr_distn - bench_distn)*log(curr_distn/bench_distn))*scale_factor

  
}



################################################################################
# Calculation of critical values.
################################################################################

#--------------------------------------------------------------------------------
# Set parameters for distribution.
#--------------------------------------------------------------------------------


# Compare to a particular chi-square distribution.
chi_2_df <- 9

# Compare on selected quantiles.
quantile_list <- c(0.50,0.75,0.90,0.95,0.975)

# Set number of realizations for the simulation exercise.
num_sims <- 1000000

# Set sample sizes.
num_bench <- 100000
num_current <- 10000
# num_current <- num_bench # Optional for symmetry (not required).

# Set the normalization factor.
scale_factor <- num_current # Imposes scale for chi-squared distribution.


#--------------------------------------------------------------------------------
# Simulation.
#--------------------------------------------------------------------------------

# Perform the simulation on the list of examples determined above.
for (example_num in 1:length(example_distn)) {
  
  # Select the benchmark distribution.
  bench_distn <- example_distn[[example_num]]
  distn_label <-example_distn_label[example_num]
  
  # Generate a sample from the multinomial distribution.
  # Generate num_sims realizations of num_current drawn from the multinomial
  # distribution with specified probabilities, normalized by the sample size.
  current_multi_distn <- t(rmultinom(n = num_sims, size = num_current, prob = bench_distn))/num_current
  
  # Repeat for the benchmark sample.
  bench_multi_distn <- t(rmultinom(n = num_sims, size = num_bench, prob = bench_distn))/num_bench
  
  #--------------------------------------------------------------------------------
  # Critical Values: Treating the benchmark as known.
  #--------------------------------------------------------------------------------
  
  # Calculate a list of KLD statistics for the current samples.
  KLD_sample_known <- apply(current_multi_distn, 1, 
                            function(curr_distn_row) {
                              KL_divergence(curr_distn_row, bench_distn, scale_factor = scale_factor)
                            })
  
  # Calculate some statistics and critical values.
  print(sprintf('Summary statistics for the %s distribution (with known benchmark distribution):', distn_label))
  print(summary(KLD_sample_known))
  print(sprintf('Quantiles of the %s distribution:', distn_label))
  print(round(quantile(KLD_sample_known, quantile_list), 5))
  
  
  #--------------------------------------------------------------------------------
  # Critical Values: Treating the benchmark as a draw from the distribution.
  #--------------------------------------------------------------------------------
  
  # Calculate a list of KLD statistics for the current samples.
  KLD_sample_draws <- numeric(num_sims)
  for (row_num in seq(num_sims)) {
    KLD_sample_draws[row_num] <- KL_divergence(current_multi_distn[row_num, ], 
                                               bench_multi_distn[row_num, ], 
                                               scale_factor = scale_factor)
  }
  
  
  # Calculate some statistics and critical values.
  print(sprintf('Summary statistics for the %s distribution (with drawn benchmark distribution):', distn_label))
  print(summary(KLD_sample_draws))
  print(sprintf('Quantiles of the %s distribution:', distn_label))
  print(round(quantile(KLD_sample_draws, quantile_list), 5))
  
}


# Compare with quantiles of the chi-squared distribution.
print(qchisq(quantile_list, df = chi_2_df))


# Notice rough agreement with the theoretical distribution.
# Also notice the slight oversize of the distributions for the skewed distributions
# These tests would reject the null of identical distributions slightly more often
# than the stated level of significance. 



################################################################################
# End
################################################################################
