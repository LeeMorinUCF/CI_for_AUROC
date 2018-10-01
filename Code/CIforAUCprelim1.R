################################################################################
# 
# Calculation of COnfidence Intervals for the Area Under the AUC
# 
# Lee Morin, Ph.D.
# Senior Data Scientist
# Capital One
# 
# October 16, 2017
# 
################################################################################
# 
# This program calculates a series of confidence intervals for the area under
# the ROC curve, for a number of proposed methods.
# 
# 
################################################################################


################################################################################
# Clearing Workspace and Declaring Packages
################################################################################

# Clear workspace.
rm(list=ls(all=TRUE))

# Set working directory.
# wd_path <- '/prod/user/sam/can/card/canada/stats/non_npi/Lee/fisher/somersd' # On unix
wd_path <- 'C:/Users/iky155/Documents/Fisher/CIforAUC' # On windows
setwd(wd_path)

# Set path to data folder.
# dataPath <- sprintf('%s/your_data_folder', wd_path)
# data_path <- sprintf('%s/data', wd_path)
not_skynet_path <- 'C:/Users/iky155/Documents/Fisher/LearningMachineLearning/NotSkyNetShowdown' # On windows
data_in_path <- sprintf('%s', not_skynet_path)

# Set path for backups (same location, for now).
# img_path <- sprintf('%s', wd_path)




################################################################################
# Loading packages
################################################################################

library(data.table)



################################################################################
# Data Loading and Verification
################################################################################


#--------------------------------------------------------------------------------
# Dataset of predictions from 'This is Not SkyNet' competition.
#--------------------------------------------------------------------------------

# Set parameters for data loading.
dataFileName <- 'SkyNetPredictions1.csv'
dataPathFileName <- sprintf('%s/%s', data_in_path, dataFileName)

# Read file.
comp_data <- read.csv(dataPathFileName)


nrow(comp_data)
colnames(comp_data)
summary(comp_data)



# Select some scores.
score_names <- colnames(comp_data)[6:7]




################################################################################
# Data Generation
################################################################################


#--------------------------------------------------------------------------------
# Binormal model.
#--------------------------------------------------------------------------------

mu_x <- 1
sigma_x <- 1
mu_y <- 2
sigma_y <- 1

n_x <- 1000
n_y <- 100

binorm_df <- data.frame(outcomes = c(rep(FALSE, n_x),rep(TRUE, n_y)),
                        scores = c(rnorm(n = n_x, mean = mu_x, sd = sigma_x), 
                                   rnorm(n = n_y, mean = mu_y, sd = sigma_y)))
# Add weights by outcome.
binorm_df[binorm_df[, 'outcomes'], 'weights'] <- 1/n_y
binorm_df[!binorm_df[, 'outcomes'], 'weights'] <- 1/n_x

# Sort by scores.
binorm_df <- binorm_df[order(binorm_df[, 'scores']), ]

# Trust but verify.
summary(binorm_df)
summary(binorm_df[binorm_df[, 'outcomes'], ])
summary(binorm_df[!binorm_df[, 'outcomes'], ])



auc_loops_binorm <- auc_loops(scores = binorm_df[, 'scores'], outcomes = binorm_df[, 'outcomes'])

auc_rank_binorm <- auc_rank(scores = binorm_df[, 'scores'], outcomes = binorm_df[, 'outcomes'])

auc_binorm_binorm <- auc_binorm(mu_0 = mu_x, sigma2_0 = sigma_x, mu_1 = mu_y, sigma2_1 = sigma_y)


A_hat_curr_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                       outcomes = binorm_df[, 'outcomes'], 
                                       weights = binorm_df[, 'weights'])

# Compare results for consistency.
auc_loops_binorm
auc_rank_binorm
auc_binorm_binorm
A_hat_curr_binom

# Set parameters for re-weighting toward the null hypothesis.
eta <- 0.1
A_0 <- 0.75
tol <- 0.01
max_iter <- 10

# df_in <- binorm_df

neg_gradient <- pos_cond_shift(scores = binorm_df[, 'scores'], 
                               outcomes = binorm_df[, 'outcomes'], 
                               weights = binorm_df[, 'weights'])
pos_gradient <- neg_cond_shift(scores = binorm_df[, 'scores'], 
                               outcomes = binorm_df[, 'outcomes'], 
                               weights = binorm_df[, 'weights'])
next_weights <- iter_weights(scores = binorm_df[, 'scores'], 
                             outcomes = binorm_df[, 'outcomes'], 
                             weights = binorm_df[, 'weights'], eta, A_0, A_hat = A_hat_curr_binom)
# Recalculate AUROC with next weights.
A_hat_next_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                       outcomes = binorm_df[, 'outcomes'], 
                                       weights = next_weights)
A_hat_curr_binom
A_hat_next_binom

# Set parameters for re-weighting toward the null hypothesis.
eta <- 0.25
A_0 <- 0.85
tol <- 0.001
max_iter <- 25


new_weights <- adj_weights(scores = binorm_df[, 'scores'], 
                           outcomes = binorm_df[, 'outcomes'], 
                           weights = binorm_df[, 'weights'], eta, A_0, A_hat = A_hat_curr_binom, tol, max_iter)

# Trust but verify.
A_hat_new_binom <- auc_weighted_loops(scores = binorm_df[, 'scores'], 
                                      outcomes = binorm_df[, 'outcomes'], 
                                      weights = new_weights)
A_hat_curr_binom
A_hat_new_binom



#--------------------------------------------------------------------------------
# Biexponential model.
#--------------------------------------------------------------------------------

lambda_x <- 3
lambda_y <- 1

biexp_df <- data.frame(outcomes = c(rep(0, n_x),rep(1, n_y)),
                       scores = c(rexp(n = n_x, rate = lambda_x), 
                                  rexp(n = n_y, rate = lambda_y)))

# Add weights by outcome.
biexp_df[outcomes, 'weights'] <- 1/n_y
biexp_df[!outcomes, 'weights'] <- 1/n_x

summary(biexp_df[, 'scores'])
summary(biexp_df[biexp_df[, 'outcomes'], 'scores'])
summary(biexp_df[!biexp_df[, 'outcomes'], 'scores'])


auc_loops_biexp <- auc_loops(scores = biexp_df[, 'scores'], outcomes = biexp_df[, 'outcomes'])

auc_rank_biexp <- auc_rank(scores = biexp_df[, 'scores'], outcomes = biexp_df[, 'outcomes'])

A_hat_curr_biexp <- auc_weighted_loops(scores = biexp_df[, 'scores'], 
                                       outcomes = biexp_df[, 'outcomes'], 
                                       weights = biexp_df[, 'weights'])

auc_loops_biexp
auc_rank_biexp
A_hat_curr_biexp




################################################################################
# Functions
################################################################################



#--------------------------------------------------------------------------------
# Plot of positive and negative score distributions
#--------------------------------------------------------------------------------

# Plots two histograms of scores by outcome. 

nBins <- 100

colTransparency <- 1/4
colObs <- rgb(0,0,1,colTransparency)
colSim <- rgb(1,0,0,colTransparency)

# Plot histogram excluding zero values.
hist(seriesObs, 
     freq = displayFreq,
     breaks = nBins, 
     ylim = c(0, maxY),
     xlim = c(0, maxX),
     col = colObs,
     main = title,
     xlab = xLabel)


#--------------------------------------------------------------------------------
# AUC calculated - arbitrary weights.
#--------------------------------------------------------------------------------

# Calculate AUROC estimate under new weights.
# A_hat_current <- auc_weighted_loops(scores, outcomes, weights)

auc_weighted_loops <- function(scores, outcomes, weights) {
  
  # Calculate weighted AUROC.
  
  # Split into positives and negatives.
  pos_scores <- scores[outcomes]
  neg_scores <- scores[!outcomes]
  
  # Separate positive and negative weights.
  pos_wts <- weights[outcomes]
  neg_wts <- weights[!outcomes]
  
  # Calculate constants.
  n_1 <- length(pos_scores)
  n_0 <- length(neg_scores)
  
  # Initialize.
  auc_loops <- 0
  
  # Test rank-ordering of each pair of observations.
  for (i in 1:n_0) {
    
    for (j in 1:n_1) {
      
      auc_loops <- auc_loops + pos_wts[j]*neg_wts[i]*(pos_scores[j] > neg_scores[i])
      
    }
    
  }
  
  # Calculate probability from number of permutations.
  # auc_loops <- auc_loops/n_0/n_1
  auc_loops
  
}




#--------------------------------------------------------------------------------
# AUC Shift - Move to nearest distribution consistent with the null hypothesis.
#--------------------------------------------------------------------------------

# Iteratively adjusts weights untill distribution is the nearest distribution 
# that is consistent with the null hypothesis.

# new_weights <- adj_weights(scores, outcomes, weights, eta, A_0, A_hat_current)

adj_weights <- function(scores, outcomes, weights, eta, A_0, A_hat, tol, max_iter) {
  
  # Initialize current state with estimated AUROC and associated weights.
  A_hat_current <- A_hat
  new_weights <- weights
  num_iter <- 0
  while (abs(A_hat_current - A_0) > tol & num_iter <= max_iter) {
    
    num_iter <- num_iter + 1
    
    # Adjust weights.
    new_weights <- iter_weights(scores, outcomes, new_weights, eta, A_0, A_hat_current)
    
    # Recalculate AUROC estimate under new weights.
    # A_hat_current <- auc_weighted(scores, outcomes, weights)
    A_hat_current <- auc_weighted_loops(scores, outcomes, new_weights)
    
    print(A_hat_current)
    
  }
  
  return(new_weights)
  
}


#--------------------------------------------------------------------------------
# AUC Shift - Iteration on distributions
#--------------------------------------------------------------------------------

# Computes a shift in positive and negative distributions toward desired AUROC.


#
# NOTE: Clean up business of positive and negative scores.
# 

# next_weights <- iter_weights(scores, outcomes, weights, eta, A_0, A_hat)

# scores, outcomes
iter_weights <- function(scores, outcomes, weights, eta, A_0, A_hat) {
  
  # Calculate step size.
  # lambda <- eta*(A_0 - A_hat)
  lambda <- - eta*(A_0 - A_hat)
  
  # Separate positive and negative scores.
  pos_scores <- scores[outcomes]
  neg_scores <- scores[!outcomes]
  
  # Separate positive and negative weights.
  pos_wts <- weights[outcomes]
  neg_wts <- weights[!outcomes]
  
  # Calculate gradient in weights.
  # Note that this assumes that data are sorted by scores (increasing order).
  pos_gradient <- neg_cond_shift(scores, outcomes, weights)
  neg_gradient <- pos_cond_shift(scores, outcomes, weights)
  
  # Calculate updated weights.
  pos_wts_next <- pos_wts*exp(lambda*pos_gradient)
  neg_wts_next <- neg_wts*exp(lambda*neg_gradient)
  
  # Normalize for unit probability mass.
  pos_wts_next <- pos_wts_next/sum(pos_wts_next)
  neg_wts_next <- neg_wts_next/sum(neg_wts_next)
  
  # Place these weights in the data frame.
  next_weights <- 0*weights
  next_weights[outcomes] <- pos_wts_next
  next_weights[!outcomes] <- neg_wts_next
  
  return(next_weights)
  
}



#--------------------------------------------------------------------------------
# AUC Shift - Positive conditional distribution
#--------------------------------------------------------------------------------

# Computes a reweighting gradient to adjust the distribution of the scores 
# corresponding to positive observations.

# pos_gradient <- neg_cond_shift(scores, outcomes, weights)

# Vector version first.
neg_cond_shift <- function(scores, outcomes, weights) {
  
  pos_gradient <- cumsum(!outcomes[order(-scores)] * weights[order(-scores)])[outcomes[order(-scores)]]
  
  # Need to reverse order, since sorting by score is in reverse.
  pos_gradient <- pos_gradient[order(seq(length(pos_gradient), 1, by = -1))]
  
  
  # Normalize.
  pos_gradient <- pos_gradient/sum(outcomes)
  
}


# neg_cond_shift <- function(scores, outcomes, weights) {
#   
#   # Separate negative weights, in reverse order.
#   neg_wts <- weights[!outcomes][order(-scores[!outcomes])]
#   
#   # Start with the opposite of outcomes sorted by scores, in reverse order.
#   # rank_neg_wts <- !outcomes[order(-scores)]
#   rank_neg_wts <- 0*weights
#   
#   # Multiply negative outcomes by weights.
#   rank_neg_wts[!outcomes[order(-scores)]] <- neg_wts
#   
#   # cumsum to get rank-weights, where ranks would normally appear.
#   rank_neg_wts <- cumsum(rank_neg_wts)
#   
#   # Re-reverse order of rank-weights.
#   # rank_neg_wts <- rank_neg_wts[order()] # How to sort.
#   
#   # Select the rank-weights corresponding to the negative outcomes.
#   pos_gradient <- rank_neg_wts[outcomes[order(-scores)]]
#   
#   # Re-reverse order of rank-weights.
#   pos_gradient <- pos_gradient[order()] # WTF?
#   
# }




#--------------------------------------------------------------------------------
# AUC Shift - Negative conditional distribution
#--------------------------------------------------------------------------------

# Computes a reweighting gradient to adjust the distribution of the scores 
# corresponding to negative observations.

# neg_gradient <- pos_cond_shift(scores, outcomes, weights)

# Vector version first.
pos_cond_shift <- function(scores, outcomes, weights) {
  
  neg_gradient <- cumsum(outcomes[order(scores)] * weights[order(scores)])[!outcomes[order(scores)]]
  
  # Normalize.
  neg_gradient <- neg_gradient/sum(!outcomes)
  
}



# pos_cond_shift <- function(scores, outcomes, weights) {
#   
#   # Separate positive weights.
#   pos_wts <- weights[outcomes]
#   
#   # Start with outcomes sorted by scores.
#   # rank_pos_wts <- outcomes[order(scores)]
#   rank_pos_wts <- 0*weights
#   
#   # Multiply positive outcomes by weights.
#   rank_pos_wts[outcomes[order(scores)]] <- pos_wts[order(scores[outcomes])]
#   
#   # cumsum to get rank-weights, where ranks would normally appear.
#   rank_pos_wts <- cumsum(rank_pos_wts)
#   
#   # Select the rank-weights corresponding to the negative outcomes.
#   neg_gradient <- rank_pos_wts[!outcomes[order(scores)]]
#   
# }
# 
# 
# # df_in <- binorm_df
# pos_cond_shift <- function(df_in) {
#   
#   # Expecting these columns:
#   # 'scores', 'outcomes', 'weights'
#   
#   # Reorder observations.
#   df_in <- df_in[order(df_in[, 'scores']), ]
#   
#   # Multiply positive outcomes by weights and accumulate sums.
#   # df_in[, 'rank_pos_wts'] <- cumsum(df_in[, 'outcomes']*df_in[, 'weights'])
#   
#   # Select the rank-weights corresponding to the negative outcomes.
#   # neg_gradient <- df_in[!df_in[, 'outcomes'], 'rank_pos_wts']
#   
#   # One-line nerd contest winner.
#   df_in[, 'neg_gradient'] <- cumsum(df_in[order(df_in[, 'scores']), 'outcomes'] *
#                                       df_in[order(df_in[, 'scores']), 'weights']) * 
#     (!df_in[order(df_in[, 'scores']), 'outcomes'])
#   
# }



#--------------------------------------------------------------------------------
# AUC Calculation 0 - Direct Computation through Loops (Illustration Only)
#--------------------------------------------------------------------------------

# Runs through permutations of pairs of (unequal) scores from the positive and 
# negative populations.

auc_loops <- function(scores, outcomes) {
  
  # Split into positives and negatives.
  pos_scores <- scores[outcomes]
  neg_scores <- scores[!outcomes]
  
  # Calculate constants.
  n_1 <- length(pos_scores)
  n_0 <- length(neg_scores)
  
  # Initialize.
  auc_loops <- 0
  
  # Test rank-ordering of each pair of observations.
  for (i in 1:n_0) {
    
    for (j in 1:n_1) {
      
      auc_loops <- auc_loops + (neg_scores[i] < pos_scores[j])
      
    }
    
  }
  
  # Calculate probability from number of permutations.
  auc_loops <- auc_loops/n_0/n_1
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 1 - Direct Computation
#--------------------------------------------------------------------------------

# Runs through permutations of pairs of (unequal) scores from the positive and 
# negative populations.

auc_slow <- function(scores, outcomes) {
  
  # Calculate constants.
  n_1 <- sum(outcomes, na.rm = FALSE)
  n_0 <- sum(!outcomes, na.rm = FALSE)
  
  # Test rank-ordering of each pair of observations.
  auc_slow <- sum(sapply(scores[!outcomes], 
                         function(x_i) {
                           sapply(scores[outcomes], 
                                  function(y_j) {
                                    x_i < y_j
                                  })
                         }))
  
  # Calculate probability from number of permutations.
  auc_slow <- auc_slow/n_0/n_1
  
}

#--------------------------------------------------------------------------------
# AUC Calculation 2 - Shortcut Using Ranks
#--------------------------------------------------------------------------------

# Follows the methodology for the equivalent Wilcoxon Signed-Rank Test.

auc_rank <- function(scores, outcomes) {
  
  # Calculate constants.
  n_1 <- sum(outcomes, na.rm = FALSE)
  n_0 <- sum(!outcomes, na.rm = FALSE)
  
  # Obtain ranks of scores.
  ranks <- rank(scores)
  
  # Calculate probability from number of permutations.
  auc_rank <- (sum(ranks[outcomes]) - n_1*(n_1 + 1)/2)/n_0/n_1
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 3 - Binormal distribution
#--------------------------------------------------------------------------------

# Closed-form expression for the special case of the binorml classification model.

auc_binorm <- function(mu_0, sigma2_0, mu_1, sigma2_1) {
  
  auc <- pnorm((mu_1 - mu_0)/sqrt(sigma2_1 + sigma2_0))
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 4 - Bi-exponential distribution
#--------------------------------------------------------------------------------

# Closed-form expression for the special case of the binorml classification model.

auc_bi_exp <- function(lambda_0, lambda_1) {
  
  auc <- 7
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 5 - Mean AUC from distribution-free, with fixed error rate
#--------------------------------------------------------------------------------

# Closed-form expression for the special case of the binorml classification model.

mean_auc_fixed_error <- function(n_0, n_1, k) {
  
  mean_auc <- 7
  
}


#--------------------------------------------------------------------------------
# P_yyx Calculation 0 - Direct Computation through Loops (Illustration Only)
#--------------------------------------------------------------------------------

# Runs through permutations of triplets of (unequal) scores from the positive and 
# negative populations.

P_yyx_loops <- function(scores, outcomes) {
  
  # Split into positives and negatives.
  pos_scores <- scores[outcomes]
  neg_scores <- scores[!outcomes]
  
  # Calculate constants.
  n_1 <- length(pos_scores)
  n_0 <- length(neg_scores)
  
  # Initialize.
  P_yyx_loops <- 0
  
  # Test rank-ordering of each pair of observations.
  for (i in 1:n_0) {
    
    for (j in 1:n_1) {
      
      for (l in 1:n_1) {
        
        P_yyx_loops <- P_yyx_loops + (neg_scores[i] < pos_scores[j] & neg_scores[i] < pos_scores[l])
        
      }
      
    }
    
  }
  
  # Calculate probability from number of permutations.
  P_yyx_loops <- P_yyx_loops/n_0/n_1/n_1
  
}


#--------------------------------------------------------------------------------
# P_yyx Calculation 1 - Direct Computation - Vectorized
#--------------------------------------------------------------------------------

# Runs through permutations of triplets of (unequal) scores from the positive and 
# negative populations.

P_yyx_vec <- function(scores, outcomes) {
  
  # Calculate constants.
  n_1 <- sum(outcomes, na.rm = FALSE)
  n_0 <- sum(!outcomes, na.rm = FALSE)
  
  # Test rank-ordering of each pair of observations.
  P_yyx_vec <- sum(sapply(scores[!outcomes], 
                          function(x_i) {
                            sapply(scores[outcomes], 
                                   function(y_j) {
                                     sapply(scores[outcomes],
                                            function(y_l) {
                                              x_i < y_j & x_i < y_l
                                            })
                                   })
                          }))
  
  # Calculate probability from number of permutations.
  P_yyx_vec <- P_yyx_vec/n_0/n_1/n_1
  
}

#--------------------------------------------------------------------------------
# P_yxx Calculation 0 - Direct Computation through Loops (Illustration Only)
#--------------------------------------------------------------------------------

# Runs through permutations of triplets of (unequal) scores from the positive and 
# negative populations.

P_yxx_loops <- function(scores, outcomes) {
  
  # Split into positives and negatives.
  pos_scores <- scores[outcomes]
  neg_scores <- scores[!outcomes]
  
  # Calculate constants.
  n_1 <- length(pos_scores)
  n_0 <- length(neg_scores)
  
  # Initialize.
  P_yxx_loops <- 0
  
  # Test rank-ordering of each pair of observations.
  for (i in 1:n_0) {
    
    for (j in 1:n_1) {
      
      for (l in 1:n_0) {
        
        P_yxx_loops <- P_yxx_loops + (neg_scores[i] < pos_scores[j] & neg_scores[l] < pos_scores[j])
        
      }
      
    }
    
  }
  
  # Calculate probability from number of permutations.
  P_yxx_loops <- P_yxx_loops/n_0/n_1/n_1
  
}


#--------------------------------------------------------------------------------
# P_yxx Calculation 1 - Direct Computation - Vectorized
#--------------------------------------------------------------------------------

# Runs through permutations of triplets of (unequal) scores from the positive and 
# negative populations.

P_yxx_vec <- function(scores, outcomes) {
  
  # Calculate constants.
  n_1 <- sum(outcomes, na.rm = FALSE)
  n_0 <- sum(!outcomes, na.rm = FALSE)
  
  # Test rank-ordering of each pair of observations.
  P_yxx_vec <- sum(sapply(scores[!outcomes], 
                          function(x_i) {
                            sapply(scores[outcomes], 
                                   function(y_j) {
                                     sapply(scores[!outcomes],
                                            function(x_l) {
                                              x_i < y_j & x_l < y_j
                                            })
                                   })
                          }))
  
  # Calculate probability from number of permutations.
  P_yxx_vec <- P_yxx_vec/n_0/n_1/n_1
  
}


#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 0 - From Variance
#--------------------------------------------------------------------------------

# Usual symmetric confidence interval.


# Calculate the bounds on the confidence interval.
confidence_interval <- function(mean, var, alpha) {
  
  z_critical <- qnorm(1 - alpha/2)
  auc_ci_binorm <- c(pnorm(mean - z_critical*var), 
                     pnorm(mean + z_critical*var))
  
}



#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 1 - Upper Bound
#--------------------------------------------------------------------------------

# Upper bound over all possible distributions.

auc_ci_upper_bound <- function(auc, n_0, n_1, alpha) {
  
  
  # Upper bound on the variance.
  var_upper_bound <- auc*(1-auc)/min(n_0, n_1)
  
  # Confidence interval.
  auc_ci_upper_bound <- confidence_interval(auc, var_upper_bound, alpha)
  
  
  
}


#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 2 - Binormal Model
#--------------------------------------------------------------------------------

# Assumes that scores come from two normal distributions.

auc_ci_binorm <- function(scores, outcomes, alpha = 0.05) {
  
  # Extract parameters from scores.
  n <- sum(outcomes)
  m <- sum(!outcomes)
  y_bar <- mean(scores[outcomes], na.rm = TRUE)
  x_bar <- mean(scores[!outcomes], na.rm = TRUE)
  s2_y <- var(scores[outcomes], na.rm = TRUE)
  s2_x <- var(scores[!outcomes], na.rm = TRUE)
  
  # Calculate the difference between distributions.
  delta <- (y_bar - x_bar)/sqrt(s2_y + s2_x)
  
  # Calculate the variance of the difference between distributions.
  var_delta <- (s2_x/m + s2_y/n)/(s2_x + s2_y) + 
    (y_bar - x_bar)^2/2/(s2_x + s2_y)^3 * (s2_x^2/(m - 1) + s2_y^2/(n - 1))
  
  # Calculate the bounds on the confidence interval.
  auc_ci_binorm <- pnorm(confidence_interval(delta, var_delta, alpha))
  
}

#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 3 - Bi-exponential Model
#--------------------------------------------------------------------------------

# Assumes that scores come from two exponential distributions.
# This is the classic version developed by Hanley and McNeil.

auc_ci_biexp <- function(auc, n_0, n_1, alpha) {
  
  # Calculate the inputs. 
  P_yyx <- auc/(2 - auc) 
  P_yxx <- 2*auc^2/(1 + auc)
  
  # Calculate the variance.
  sigma2_auc <- (auc*(1-auc) + (n_1 - 1)*(P_yyx - auc^2) + (n_0 - 1)*(P_yxx - auc^2))/n_0/n_1
  
  # Confidence interval.
  auc_ci_biexp <- confidence_interval(auc, sigma2_auc, alpha)
  
  
}

#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 4 - Distribution-free (unrestricted)
#--------------------------------------------------------------------------------

# Makes no distributional assumptions, in the sense that every possible outcome 
# is given equal weight.

auc_ci_distn_free <- function(auc, P_yyx, P_yxx, n_0, n_1, alpha) {
  
  # Calculate the variance.
  sigma2_auc <- (auc*(1-auc) + (n_1 - 1)*(P_yyx - auc^2) + (n_0 - 1)*(P_yxx - auc^2))
  
  # Confidence interval.
  auc_ci_distn_free <- confidence_interval(auc, sigma2_auc, alpha)
  
  
}



#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 5 - Distribution-free (fixed error rate)
#--------------------------------------------------------------------------------

# Makes no distributional assumptions, in the sense that every possible outcome 
# is given equal weight, subject to samples with the same aggregate number of errors.



################################################################################
# Calculation of standard errors.
################################################################################

summary(comp_data)

sel_obsns <- 1:100000

target <- comp_data[sel_obsns, 'Target']
# scores <- 7

z_1 <- comp_data[sel_obsns, 'Darren_Lin']
z_2 <- comp_data[sel_obsns, 'Jon_S_Lee']

# y_1 <- comp_data[comp_data[, 'Target'] == 1, 'Darren_Lin']
# x_1 <- comp_data[comp_data[, 'Target'] == 0, 'Darren_Lin']
# y_2 <- comp_data[comp_data[, 'Target'] == 1, 'Jon_S_Lee']
# x_2 <- comp_data[comp_data[, 'Target'] == 0, 'Jon_S_Lee']


y_1 <- z_1[target == 1]
x_1 <- z_1[target == 0]
y_2 <- z_2[target == 1]
x_2 <- z_2[target == 0]


auc_1 <- auc_loops(scores = z_1, outcomes = target)
auc_2 <- auc_loops(scores = z_2, outcomes = target)







################################################################################
# End
################################################################################
