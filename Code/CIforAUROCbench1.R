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
# the ROC curve, for a number of existing methods.
# It serves as a benchmark aginst which to compare new methods.
# 
# TODO: Fix biexponential example.
# TODO: Complete calculation of remaining variances and confidence intervals. Done.
# TODO: COde up method with fixed error rate.
# TODO: Update calculation of higher-order terms, such as Pyyx and Pyxx. Done.
# 
################################################################################



################################################################################
# Data Generation
################################################################################

#--------------------------------------------------------------------------------
# Binormal Model Specification
#--------------------------------------------------------------------------------

# Specifies parameters to generate a required series of true AUROC values.

# binorm_spec_list <- binorm_spec(auroc_vec)

binorm_spec <- function(auroc_vec) {
  
  num_auroc <- length(auroc_vec)
  
  # Calculate mu_y with required AUROC.
  mu_y_auroc <- qnorm(auroc_vec)
  
  binorm_spec_list <- data.frame(mu_x = rep(0, num_auroc),
                                 mu_y = mu_y_auroc,
                                 sigma_x = rep(1/sqrt(2), num_auroc),
                                 sigma_y = rep(1/sqrt(2), num_auroc))
  
  return(binorm_spec_list)
  
}

# # Test: 
# auroc_vec <- c(0.5, 0.6, 0.75, 0.8, 0.95, 0.975)
# binorm_spec_list <- binorm_spec(auroc_vec)
# auc_binorm_true <- rep(NA, length(auroc_vec))
# for (auroc_num in 1:length(auroc_vec)) {
#   auc_binorm_true[auroc_num] <- auc_binorm(mu_0 = binorm_spec_list[auroc_num, 'mu_x'], 
#                                            sigma_0 = binorm_spec_list[auroc_num, 'sigma_x'], 
#                                            mu_1 = binorm_spec_list[auroc_num, 'mu_y'], 
#                                            sigma_1 = binorm_spec_list[auroc_num, 'sigma_y'])
# }
# auc_binorm_true


#--------------------------------------------------------------------------------
# Binormal Model Generation
#--------------------------------------------------------------------------------

binorm_gen <- function(n_x, n_y, mu_x, mu_y, sigma_x, sigma_y) {
  
  binorm_df <- data.frame(outcomes = c(rep(FALSE, n_x),rep(TRUE, n_y)),
                          scores = c(rnorm(n = n_x, mean = mu_x, sd = sigma_x), 
                                     rnorm(n = n_y, mean = mu_y, sd = sigma_y)))
  # Add weights by outcome.
  binorm_df[binorm_df[, 'outcomes'], 'weights'] <- 1/n_y
  binorm_df[!binorm_df[, 'outcomes'], 'weights'] <- 1/n_x
  
  # Sort by scores.
  binorm_df <- binorm_df[order(binorm_df[, 'scores']), ]
  
  
}

#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# Biexponential Model
#--------------------------------------------------------------------------------

biexp_gen <- function(n_x, n_y, lambda_x, lambda_y) {
  
  biexp_df <- data.frame(outcomes = c(rep(0, n_x),rep(1, n_y)),
                         scores = c(rexp(n = n_x, rate = lambda_x), 
                                    rexp(n = n_y, rate = lambda_y)))
  
  # Add weights by outcome.
  biexp_df[outcomes, 'weights'] <- 1/n_y
  biexp_df[!outcomes, 'weights'] <- 1/n_x
  
  # Sort by scores.
  biexp_df <- biexp_df[order(biexp_df[, 'scores']), ]
  
  
}

#--------------------------------------------------------------------------------


#--------------------------------------------------------------------------------
# BiPower Model
#--------------------------------------------------------------------------------

bipower_gen <- function(n_x, n_y, 
                        x_min, alpha_x, y_min, gamma_y) {
  
  bipower_df <- data.frame(outcomes = c(rep(FALSE, n_x),rep(TRUE, n_y)),
                         scores = c(rpldis(n = n_x, xmin = x_min, alpha = alpha_x), 
                                    rpldis(n = n_y, xmin = y_min, alpha = gamma_y)))
  
  # print(summary(bipower_df))
  # print(summary(bipower_df[bipower_df[, 'outcomes'] == 0, ]))
  # print(summary(bipower_df[bipower_df[, 'outcomes'] == 1, ]))
  
  # Add weights by outcome.
  bipower_df[bipower_df[, 'outcomes'], 'weights'] <- 1/n_y
  bipower_df[!bipower_df[, 'outcomes'], 'weights'] <- 1/n_x
  
  # Sort by scores.
  bipower_df <- bipower_df[order(bipower_df[, 'scores']), ]
  
  
}

#--------------------------------------------------------------------------------




################################################################################
# AUROC Calculation
################################################################################

#--------------------------------------------------------------------------------
# AUC calculated - arbitrary weights.
#--------------------------------------------------------------------------------

# Calculate AUROC estimate under new weights.
# A_hat_current <- auc_weighted_loops_dt(dt_in)

auc_weighted_loops_dt <- function(dt_in) {
  
  # Initialize.
  auc_loops <- 0
  
  # Test rank-ordering of each pair of observations.
  for (i in which(binorm_dt[, !outcomes])) {
    
    for (j in which(binorm_dt[, outcomes])) {
      
      auc_loops <- auc_loops + 
        binorm_dt[j, weights]*binorm_dt[i, weights] * 
        (binorm_dt[j, scores] > binorm_dt[i, scores])
      
    }
    print(auc_loops)
    
  }
  
  # Calculate probability from number of permutations.
  auc_loops
  
}


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
  # auc_loops <- auc_loops/n_0/n_1 # Not required when weights are specified.
  auc_loops
  
}

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

# Closed-form expression for the special case of the binormal classification model.

auc_binorm <- function(mu_0, sigma_0, mu_1, sigma_1) {
  
  auc <- pnorm((mu_1 - mu_0)/sqrt(sigma_1^2 + sigma_0^2))
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 4 - Bi-exponential distribution
#--------------------------------------------------------------------------------

# Closed-form expression for the special case of the biexponential classification model.

auc_bi_exp <- function(lambda_0, lambda_1) {
  
  auc <- 7
  
}



#--------------------------------------------------------------------------------
# AUC Calculation 5 - Bi-power-law distribution
#--------------------------------------------------------------------------------

# Closed-form expression for the special case of the bipower-law classification model.

auc_bi_power <- function(x_min, alpha_x, 
                         y_min, gamma_y) {
  
  # Depends on lower bound.
  if (y_min > x_min) {
    
    auc <- 1 - (gamma_y - 1) / ((gamma_y - 1) + (alpha_x - 1) ) * 
      (x_min / y_min) ^ (alpha_x - 1)
    
  } else {
    
    auc <- (alpha_x - 1) / ((gamma_y - 1) + (alpha_x - 1) ) * 
      (y_min / x_min) ^ (gamma_y - 1)
    
  }
  
  
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 6 - Mean AUC from distribution-free, with fixed error rate
#--------------------------------------------------------------------------------

# Closed-form expression for the distribution-free case, with fixed error rate.

mean_auc_fixed_error <- function(n_0, n_1, k, max_n = 1000) {
  
  # Calculate magic ratio of counts of combinations.
  n <- n_0 + n_1
  count_ratio <- Z_ratio(i = 1, n_in = n, k_in = k, max_n = max_n)
  
  mean_auc <- 1 - k/(n_0 + n_1) - (n_0 - n_1)^2*(n_0 + n_1 + 1)/4/n_0/n_1 * 
    (k/(n_0 + n_1) - count_ratio)
  
}


#--------------------------------------------------------------------------------
# AUC Calculation 6.1 - Count ratio for AUC from distribution-free, with fixed error rate
#--------------------------------------------------------------------------------

# Literal application of formula gives numerical overflow for reasonable sample sizes.

Z_ratio <- function(i, n_in, k_in, max_n = 1000) {
  
  # Calculate requred ratios.
  error_rate <- k_in/n_in
  
  # Cap sample size to avoid overflow.
  n <- min(n_in, max_n)
  k <- round(error_rate*n, 0)
  
  
  count_ratio <- sum(sapply(seq(0, k - i), function(x) {choose(n + 1 - i, x)})) / 
    sum(sapply(seq(0, k), function(x) {choose(n + 1, x)}))
  
  return(count_ratio)
}



################################################################################
# AUROC Higher Order Terms
################################################################################

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



################################################################################
# AUROC Variance Calculation
################################################################################

# For most, these are calculated within the confidence interval calculations below.

# Sun-Xu fast calculation of Hanley-McNeil

# var_auroc_fast <- var_auroc_fast_calc(dt_in)

var_auroc_fast_calc <- function(dt_in) {
  
  
  # Get parameters from inputs.
  n_y <- dt_in[outcomes == TRUE, .N]
  n_x <- dt_in[outcomes == FALSE, .N]
  
  # Calculate ranks.
  dt_in[, rank_z := rank(scores)]
  # print(head(dt_in[, rank_z], 100))
  dt_in[outcomes == TRUE, rank_y_x := rank(scores)]
  # print(head(dt_in[outcomes == TRUE, rank_y_x], 100))
  dt_in[outcomes == FALSE, rank_y_x := rank(scores)]
  # print(head(dt_in[outcomes == FALSE, rank_y_x], 100))
  # print(head(dt_in[, rank_y_x], 100))
  
  # Estimate AUROC.
  auroc_hat_est <- dt_in[outcomes == TRUE, ( sum(rank_z) - 0.5*n_y*(n_y + 1) )/(n_y*n_x)]
  # print('A_hat = ')
  # print(auroc_hat_est)
  
  # Estimate variance.
  dt_in[outcomes == TRUE, var_y_x := (rank_z - rank_y_x)/n_x]
  dt_in[outcomes == FALSE, var_y_x := 1 - (rank_z - rank_y_x)/n_y]
  
  # Put it together.
  # var_auroc_fast <- var(dt_in[outcomes == TRUE, var_y_x] - auroc_hat_est) +
  #   var(dt_in[outcomes == FALSE, var_y_x] - auroc_hat_est)
  var_auroc_fast <- dt_in[outcomes == TRUE, sum((var_y_x - auroc_hat_est)^2) / (n_y - 1)] +
    dt_in[outcomes == FALSE, sum((var_y_x - auroc_hat_est)^2) / (n_x - 1)]
  
  return(var_auroc_fast)
  
}


################################################################################
# AUROC Confidence Interval Calculation
################################################################################



#--------------------------------------------------------------------------------
# Standard Error of the AUC Calculation 0 - From Variance
#--------------------------------------------------------------------------------

# Usual symmetric confidence interval.


# Calculate the bounds on the confidence interval.
confidence_interval <- function(mean, var, alpha) {
  
  z_critical <- qnorm(1 - alpha/2)
  # auc_ci_binorm <- c(pnorm(mean - z_critical*var), 
  #                    pnorm(mean + z_critical*var))
  auc_ci_binorm <- c(mean - z_critical*sqrt(var), 
                     mean + z_critical*sqrt(var))
  
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


auc_ci_fixed_error <- function(n_0, n_1, k, max_n = 1000, mean_auc, alpha) {
  
  
  # Calculate magic ratio of counts of combinations.
  n <- n_0 + n_1
  Z_1 <- Z_ratio(i = 1, n_in = n, k_in = k, max_n = max_n)
  Z_2 <- Z_ratio(i = 2, n_in = n, k_in = k, max_n = max_n)
  Z_3 <- Z_ratio(i = 3, n_in = n, k_in = k, max_n = max_n)
  Z_4 <- Z_ratio(i = 4, n_in = n, k_in = k, max_n = max_n)
  
  # Calculate other required polynomials.
  T_0 <- 3*((n_1 - n_0)^2 + n_1 + n_0) + 2
  Q_0 <- (n_0 + n_1 + 1)*T_0*k^2 + 
    ( ( - 3*n_0^2 + 3*n_1*n_0 + 3*n_1 + 1)*T_0 - 12*(3*n_1*n_0 + n_1 + n_0) - 8 )*k +
    ( - 3*n_1^2 + 7*n_1 + 10*n_0 + 3*n_1*n_0 + 10)*T_0 - 
    4*(3*n_1*n_0 + n_1 + n_0 + 1)
  Q_1 <- T_0*k^3 + 3*(n_1 - 1)*T_0*k^2 + 
    ( ( - 3*n_0^2 + 3*n_1*n_0 - 3*n_1 + 8)*T_0 - 6*(6*n_1*n_0 + n_1 + n_0) )*k +
    ( - 3*n_1^2 + 7*(n_1 + n_0) + 3*n_1*n_0)*T_0 - 
    4*(3*n_1*n_0 + n_1 + n_0 + 1)
  
  
  # Calculate the variance.
  sigma2_auc <- (n_0 + n_1 + 1)*(n_0 + n_1)*(n_0 + n_1 - 1)* T * 
    ((n_0 + n_1 - 2)*Z_4 - (2*n_1 - n_0 + 3*k - 10)*Z_3) /72/n_0^2/n_1^2 + 
    ((n_0 + n_1 + 1)*(n_0 + n_1) * T * 
       ( n_1^2 - n_1*n_0 + 3*k*n_1 - 5*n_1 + 2*k^2 - n_0*k + 12 - 9*k )*Z_2) /48/n_0^2/n_1^2 - 
    (n_0 + n_1 + 1)^2(n_0 - n_1)^4*Z_1^2 /16/n_0^2/n_1^2 - 
    (n_0 + n_1 + 1)*Q_1*Z_1 /72/n_0^2/n_1^2 + 
    k*Q_0 /144/n_0^2/n_1^2
  
  # Confidence interval.
  auc_ci_distn_free <- confidence_interval(mean_auc, sigma2_auc, alpha)
  
  
}





################################################################################
# End
################################################################################
