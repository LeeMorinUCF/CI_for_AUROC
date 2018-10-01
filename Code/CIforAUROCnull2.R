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
# This program calculates confidence intervals for the area under
# the ROC curve, under an assumed null hypothesis.
# 
# TODO: Consider optimal weighting of adjustments to positive and negative distributions
#       (or modify for distance between joint distributions). Done.
# TODO: Modify distance metrics for joint distributions. Done.
# TODO: Reconsider calculation of distribution shift functions. Done.
# TODO: Function to plot distributions within stages. Done.
# TODO: Function to check for fixed points to check for convergence.Done.
# TODO: Graphical demonstration of shifting distributions. Done.
#   May need literal mapping of observations to bins and summing of weights. Done.
# TODO: Simulate distance from optimized distribution. Done.
# TODO: Solving for AUROC such that distance is appropriate quantile. Done.
#       Problem: 95% quantile of distance implies huge range of AUROC.
# TODO: Find method to estimate appropriate quantile of distance.
# TODO: Estimate fixed error rate.
# 
################################################################################



################################################################################
# Simulation: Distance Metrics and AUROC for Confidence Intervals
################################################################################

#--------------------------------------------------------------------------------
# Power curves for a series of confidence intervals
#--------------------------------------------------------------------------------

# binorm_power <- binorm_power_sim(n_x, n_y, auroc_vec, auroc_grid, num_sims,
#                                  metric, eta, tol, max_iter, 
#                                  num_quantiles, num_boots, alpha)


binorm_power_sim <- function(n_x, n_y, auroc_vec, auroc_grid, num_sims,
                             metric, eta, tol, max_iter, 
                             num_quantiles, num_boots, alpha) {
  
  # Calculate parameters for desired list of A_0 values.
  binorm_null_spec_list <- binorm_spec(auroc_vec)
  num_null_auroc <- length(auroc_vec)
  num_true_auroc <- length(auroc_grid)
  
  # Initialize storage matrix for coverage rates.
  stats_list <- c('binorm_true', 'binorm_est', 'upper', # 'distn_free', 
                  'auroc_delong', 'auroc_boot', 'auroc_fixed_error', 'auroc_root',
                  'dist_sim', 'auroc_sim', 'dist_sim_adj', 'auroc_sim_adj')
  dist_stats_list <- stats_list[substring(stats_list, 1,5) == 'dist_']
  auroc_stats_list <- stats_list[!(stats_list %in% dist_stats_list) & 
                                   (stats_list != 'auroc_sim_adj')]
  
  binorm_power <- expand.grid(auroc_grid, auroc_vec, stats_list)
  binorm_power <- data.frame(A_0 = binorm_power[, 'Var2'],
                             A_true = binorm_power[, 'Var2'] + binorm_power[, 'Var1'],
                             stat = binorm_power[, 'Var3'],
                             power = rep(NA,nrow(binorm_power)))
  
  for (null_auroc_num in 1:num_null_auroc) {
    
    
    A_0 <- auroc_vec[null_auroc_num]
    
    print(sprintf('Running simulation for null AUROC value %f.', A_0))
    
    
    # mu_x_null = binorm_null_spec_list[null_auroc_num, 'mu_x']
    # sigma_x_null = binorm_null_spec_list[null_auroc_num, 'sigma_x']
    # mu_y_null = binorm_null_spec_list[null_auroc_num, 'mu_y']
    # sigma_y_null = binorm_null_spec_list[null_auroc_num, 'sigma_y']
    
    
    
    binorm_true_spec_list <- binorm_spec(A_0 + auroc_grid)
    
    for (true_auroc_num in 1:num_true_auroc) {
      
      A_true <- auroc_grid[true_auroc_num]
      
      mu_x_true = binorm_true_spec_list[true_auroc_num, 'mu_x']
      sigma_x_true = binorm_true_spec_list[true_auroc_num, 'sigma_x']
      mu_y_true = binorm_true_spec_list[true_auroc_num, 'mu_y']
      sigma_y_true = binorm_true_spec_list[true_auroc_num, 'sigma_y']
      
      print(sprintf('Running simulation for null AUROC value %f.', A_true))
      
      # Run simulations for each of these values.
      for (sim_num in 1:num_sims) {
        
        print(sprintf('Running simulation realization number %d.', sim_num))
        
        # Generate sample of data.
        binorm_dt <- data.table(binorm_gen(n_x, n_y, 
                                           mu_x_true, mu_y_true, sigma_x_true, sigma_y_true))
        
        # For this realization, generate the menu of confidence intervals.
        binorm_ci_df <- binorm_ci_calc(binorm_dt, 
                                       mu_x_true, mu_y_true, sigma_x_true, sigma_y_true,
                                       metric, eta, A_0, tol, max_iter, 
                                       num_quantiles, num_boots, alpha)
        
        
        #--------------------------------------------------------------------------------
        # Sepearate operations for AUROC vs distance.
        #--------------------------------------------------------------------------------
        
        # AUROC first, with knowm A_0.
        binorm_power[binorm_power[, 'A_0'] == A_0 & 
                       binorm_power[, 'A_true'] == A_true & 
                       binorm_power[, 'stat'] %in% auroc_stats_list, 'power'] <- 
          binorm_power[binorm_power[, 'A_0'] == A_0 & 
                         binorm_power[, 'A_true'] == A_true & 
                         binorm_power[, 'stat'] %in% auroc_stats_list, 'power'] + 
          (A_0 < binorm_ci_df[binorm_ci_df[, 'stat'] %in% auroc_stats_list, 'cu_l'] | 
             A_0 > binorm_ci_df[binorm_ci_df[, 'stat'] %in% auroc_stats_list, 'cu_u'])
        
        #--------------------------------------------------------------------------------
        # For distance metric, need to calculate optimized distance and quantiles.
        #--------------------------------------------------------------------------------
        
        
        # Calculate distribution for actual observations (uniform).
        quantile_df <- distn_quantiles(binorm_dt, num_quantiles)
        distn_table <- quantile_assign(binorm_dt, quantile_cuts = quantile_df, type = 'observed')
        
        
        # Find null distribution.
        adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
                       display = FALSE, ylim = NA)
        
        # Calculate distribution for optimized observations.
        distn_table[distn_table[, 'outcomes'] == TRUE, 'distn_sim'] <- 
          binorm_dt[outcomes == TRUE, sum(opt_weights), by = quantile][, V1]
        distn_table[distn_table[, 'outcomes'] == FALSE, 'distn_sim'] <- 
          binorm_dt[outcomes == FALSE, sum(opt_weights), by = quantile][, V1]
        
        # Calculate distance to "true" distribution, as in the closest distribution
        # for which the null is satisfied.
        auroc_dist_sim <- auroc_dist(outcomes = distn_table[, 'outcomes'], 
                                     weights_1 = distn_table[, 'distn_obs'], 
                                     weights_2 = distn_table[, 'distn_sim'], 
                                     metric = metric)
        
        # Check for power and add the count to the total.
        
        binorm_power[binorm_power[, 'A_0'] == A_0 & 
                       binorm_power[, 'A_true'] == A_true & 
                       binorm_power[, 'stat'] %in% dist_stats_list, 'power'] <- 
          binorm_power[binorm_power[, 'A_0'] == A_0 & 
                         binorm_power[, 'A_true'] == A_true & 
                         binorm_power[, 'stat'] %in% dist_stats_list, 'power'] + 
          (auroc_dist_sim > binorm_ci_df[binorm_ci_df[, 'stat'] %in% dist_stats_list, 'cu_u'])
        
        
        #--------------------------------------------------------------------------------
        # Check power for CI optimized to fit null hypothesis.
        # Test whether Estimated value is within confidence interval.
        #--------------------------------------------------------------------------------
        
        # Calculate AUROC from sample.
        binorm_dt[, rank_z := rank(scores)]
        A_hat <- binorm_dt[outcomes == TRUE, ( sum(rank_z) - 0.5*n_y*(n_y + 1) )/(n_y*n_x)]
        
        # AUROC first, with knowm A_0.
        binorm_power[binorm_power[, 'A_0'] == A_0 & 
                       binorm_power[, 'A_true'] == A_true & 
                       binorm_power[, 'stat'] == 'auroc_sim_adj', 'power'] <- 
          binorm_power[binorm_power[, 'A_0'] == A_0 & 
                         binorm_power[, 'A_true'] == A_true & 
                         binorm_power[, 'stat'] == 'auroc_sim_adj', 'power'] + 
          (A_hat < binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim_adj', 'cu_l'] | 
             A_hat > binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim_adj', 'cu_u'])
        
        
      }
      
      
    }
    
  }
  
  # Normalize counts by number of realizations.
  binorm_power[, 'power'] <- binorm_power[, 'power']/num_sims
  
  return(binorm_power)
  
}


#--------------------------------------------------------------------------------
# Coverage rates for a series of confidence intervals
#--------------------------------------------------------------------------------

# binorm_coverage <- binorm_coverage_sim(n_x, n_y, auroc_vec, num_sims,
#                             metric, eta, tol, max_iter, 
#                             num_quantiles, num_boots, alpha)

binorm_coverage_sim <- function(n_x, n_y, auroc_vec, num_sims,
                           metric, eta, tol, max_iter, 
                           num_quantiles, num_boots, alpha) {
  
  
  # Calculate parameters for desired list of A_0 values.
  binorm_spec_list <- binorm_spec(auroc_vec)
  num_auroc <- length(auroc_vec)
  
  # Initialize storage matrix for coverage rates.
  stats_list <- c('binorm_true', 'binorm_est', 'upper', # 'distn_free', 
                  'auroc_delong', 'auroc_boot', 'auroc_fixed_error', 'auroc_root',
                  'dist_sim', 'auroc_sim', 'dist_sim_adj', 'auroc_sim_adj')
  dist_stats_list <- stats_list[substring(stats_list, 1,5) == 'dist_']
  auroc_stats_list <- stats_list[!(stats_list %in% dist_stats_list)]
  binorm_coverage <- data.frame(A_0 = rep(0, num_auroc*length(stats_list)),
                                stat = rep(stats_list, num_auroc),
                                coverage = rep(0, num_auroc*length(stats_list)))
  
  for (auroc_num in 1:num_auroc) {
    
    
    A_0 <- auroc_vec[auroc_num]
    
    print(sprintf('Running simulation for AUROC value %f.', A_0))
    
    auroc_index <- (length(stats_list)*(auroc_num - 1) + 1) : (length(stats_list)*auroc_num)
    binorm_coverage[auroc_index, 'A_0'] <- A_0
    
    mu_x = binorm_spec_list[auroc_num, 'mu_x']
    sigma_x = binorm_spec_list[auroc_num, 'sigma_x']
    mu_y = binorm_spec_list[auroc_num, 'mu_y']
    sigma_y = binorm_spec_list[auroc_num, 'sigma_y']
    
    # Run simulations for each of these values.
    for (sim_num in 1:num_sims) {
      
      print(sprintf('Running simulation realization number %d.', sim_num))
      
      # Generate sample of data.
      binorm_dt <- data.table(binorm_gen(n_x, n_y, mu_x, mu_y, sigma_x, sigma_y))
      
      # For this realization, generate the menu of confidence intervals.
      binorm_ci_df <- binorm_ci_calc(binorm_dt, 
                                     mu_x, sigma_x, mu_y, sigma_y,
                                     metric, eta, A_0, tol, max_iter, 
                                     num_quantiles, num_boots, alpha)
      print('binorm_ci_df = ')
      print(binorm_ci_df)
      
      # Check for coverage and add the count to the total.
      # binorm_coverage[auroc_index, 'coverage'] <- 
      #   binorm_coverage[auroc_index, 'coverage'] + 
      #   binorm_ci_df[, 'cu_l'] < A_0 & A_0 < binorm_ci_df[, 'cu_u']
      
      
      #--------------------------------------------------------------------------------
      # Sepearate operations for AUROC vs distance.
      #--------------------------------------------------------------------------------
      
      # AUROC first, with knowm A_0.
      binorm_coverage[binorm_coverage[, 'A_0'] == A_0 & 
                        binorm_coverage[, 'stat'] %in% auroc_stats_list, 'coverage'] <- 
        binorm_coverage[binorm_coverage[, 'A_0'] == A_0 & 
                          binorm_coverage[, 'stat'] %in% auroc_stats_list, 'coverage'] + 
        (binorm_ci_df[binorm_ci_df[, 'stat'] %in% auroc_stats_list, 'cu_l'] < A_0 & 
        A_0 < binorm_ci_df[binorm_ci_df[, 'stat'] %in% auroc_stats_list, 'cu_u'])
      
      #--------------------------------------------------------------------------------
      # For distance metric, need to calculate optimized distance and quantiles.
      #--------------------------------------------------------------------------------
      
      
      # Calculate distribution for actual observations (uniform).
      quantile_df <- distn_quantiles(binorm_dt, num_quantiles)
      distn_table <- quantile_assign(binorm_dt, quantile_cuts = quantile_df, type = 'observed')
      
      
      # Find null distribution.
      adj_weights_dt(binorm_dt, eta, A_0, tol, max_iter, 
                     display = FALSE, ylim = NA)
      
      # Calculate distribution for optimized observations.
      distn_table[distn_table[, 'outcomes'] == TRUE, 'distn_sim'] <- 
        binorm_dt[outcomes == TRUE, sum(opt_weights), by = quantile][, V1]
      distn_table[distn_table[, 'outcomes'] == FALSE, 'distn_sim'] <- 
        binorm_dt[outcomes == FALSE, sum(opt_weights), by = quantile][, V1]
      
      # Calculate distance to "true" distribution, as in the closest distribution
      # for which the null is satisfied.
      auroc_dist_sim <- auroc_dist(outcomes = distn_table[, 'outcomes'], 
                                   weights_1 = distn_table[, 'distn_obs'], 
                                   weights_2 = distn_table[, 'distn_sim'], 
                                   metric = metric)
      print('dist = ')
      print(auroc_dist_sim)
      
      # Check for coverage and add the count to the total.
      
      binorm_coverage[binorm_coverage[, 'A_0'] == A_0 & 
                        binorm_coverage[, 'stat'] %in% dist_stats_list, 'coverage'] <- 
        binorm_coverage[binorm_coverage[, 'A_0'] == A_0 & 
                          binorm_coverage[, 'stat'] %in% dist_stats_list, 'coverage'] + 
        (auroc_dist_sim < binorm_ci_df[binorm_ci_df[, 'stat'] %in% dist_stats_list, 'cu_u'])
      
      # Print progress report.
      print('Coverage rate, so far:')
      print(binorm_coverage[, 'coverage']/sim_num)
      
      
    }
    
  }
  
  # Convert coverage frequency to coverage rate.
  binorm_coverage[, 'coverage'] <- binorm_coverage[, 'coverage']/num_sims
  
  return(binorm_coverage)
  
  
}

#--------------------------------------------------------------------------------
# Calculation of a series of confidence intervals
#--------------------------------------------------------------------------------

# binorm_ci_df <- binorm_ci_calc(dt_in, 
#                            mu_x, sigma_x, mu_y, sigma_y,
#                            metric, eta, A_0, tol, max_iter, 
#                            num_quantiles, num_boots, alpha)

binorm_ci_calc <- function(dt_in, 
                           mu_x, sigma_x, mu_y, sigma_y,
                           metric, eta, A_0, tol, max_iter, 
                           num_quantiles, num_boots, alpha) {
  
  
  # Initialize.
  stats_list <- c('binorm_true', 'binorm_est', 'upper', # 'distn_free', 
                  'auroc_delong', 'auroc_boot', 'auroc_fixed_error', 'auroc_root',
                  'dist_sim', 'auroc_sim', 'dist_sim_adj', 'auroc_sim_adj')
  binorm_ci_df <- data.frame(stat = stats_list,
                             cu_u = numeric(length(stats_list)),
                             cu_l = numeric(length(stats_list)))
  
  # Get parameters from inputs.
  n_y <- dt_in[outcomes == TRUE, .N]
  n_x <- dt_in[outcomes == FALSE, .N]
  
  # Need the estimated AUROC for midpoint of z-statistics.
  dt_in[, rank_sim := rank(scores)]
  auroc_hat_est <- dt_in[outcomes == TRUE, ( sum(rank_sim) - 0.5*n_y*(n_y + 1) )/(n_y*n_x)]
  
  # Standard confidence interval: True values of parameters.
  auc_ci_binorm_true <- pnorm(confidence_interval(mean = (mu_y - mu_x) / 
                                                    sqrt(sigma_x^2 + sigma_y^2), 
                                                  var = (sigma_x^2/n_x + sigma_y^2/n_y) / 
                                                    (sigma_x^2 + sigma_y^2), 
                                                  alpha))
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'binorm_true', 'cu_u'] <- auc_ci_binorm_true[2]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'binorm_true', 'cu_l'] <- auc_ci_binorm_true[1]

  # Standard confidence interval: Estimates of parameters.
  auc_ci_binorm <- auc_ci_binorm(scores = dt_in[, scores], 
                                 outcomes = dt_in[, outcomes], 
                                 alpha = alpha)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'binorm_est', 'cu_u'] <- auc_ci_binorm[2]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'binorm_est', 'cu_l'] <- auc_ci_binorm[1]
  
  
  # Upper bound on confidence interval.
  auc_ci_upper_bound <- auc_ci_upper_bound(auc = auroc_hat_est, n_0 = n_x, n_1 = n_y, alpha)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'upper', 'cu_u'] <- auc_ci_upper_bound[2]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'upper', 'cu_l'] <- auc_ci_upper_bound[1]
  
  
  # Standard confidence interval: Distribution-free.
  # Eliminated for computational efficiency.
  
  # From package pROC.
  rocobj <- roc(dt_in[, outcomes], dt_in[, scores])
  auc_ci_delong <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "delong", 
                          boot.n = num_boots, boot.stratified = TRUE, progress = 'none')
  auc_ci_boot <- ci.auc(roc = rocobj, conf.level = 1 - alpha, method = "bootstrap", 
                        boot.n = num_boots, boot.stratified = TRUE, progress = 'none')
  
  
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_delong', 'cu_u'] <- auc_ci_delong[3]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_delong', 'cu_l'] <- auc_ci_delong[1]
  
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_boot', 'cu_u'] <- auc_ci_boot[3]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_boot', 'cu_l'] <- auc_ci_boot[1]
  
  
  # Confidence interval with fixed error rate.
  k <- (1 - auroc_hat_est)*(n_y + n_x)
  max_n <- 1000
  mean_auc <- mean_auc_fixed_error(n_x, n_y, k, max_n = max_n)
  auc_ci_fixed_error <- auc_ci_fixed_error(n_x, n_y, k, max_n = max_n, mean_auc, alpha)
  
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_fixed_error', 'cu_u'] <- auc_ci_fixed_error[2]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_fixed_error', 'cu_l'] <- auc_ci_fixed_error[1]
  
  
  #--------------------------------------------------------------------------------
  # First, run simulation with actual weights.
  #--------------------------------------------------------------------------------
  
  # Set optimal weights to actual weights.
  dt_in[, opt_weights := weights]
  
  # Simulation confidence intervals.
  auroc_sim_df <- auroc_sim(dt_in, metric = metric, num_quantiles, num_boots)
  # print(summary(auroc_sim_df))
  
  # auroc_sim_df <- auroc_sim_df[order(auroc_sim_df[, 'distance'])]
  
  # Store the corresponding confidence intervals.
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'dist_sim', 'cu_u'] <- quantile(auroc_sim_df[, 'distance'], 1 - alpha)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'dist_sim', 'cu_l'] <- 0
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim', 'cu_u'] <- quantile(auroc_sim_df[, 'auroc'], 1 - alpha/2)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim', 'cu_l'] <- quantile(auroc_sim_df[, 'auroc'], alpha/2)

  
  #--------------------------------------------------------------------------------
  # Next, determine the AUROC levels within critical distance of the original sample.
  #--------------------------------------------------------------------------------
  
  # Set the rate to the critical value of distance.
  D_tail <- quantile(auroc_sim_df[, 'distance'], 1 - alpha)
  auc_ci_root <- c(min(auroc_sim_df[auroc_sim_df[, 'distance'] < D_tail, 'auroc']), 
                   max(auroc_sim_df[auroc_sim_df[, 'distance'] < D_tail, 'auroc']))
  
  # First, test the distance function for good behaviour.
  # Set the rate to the average distance.
  D_bar <- quantile(auroc_sim_df[, 'distance'], 0.50)
  dist_diff_u <- dist_A_lim_diff(A_lim = auc_ci_root[2], 
                               D_bar, dt_in = dt_in, metric, num_quantiles, min_eps = 10^(-4), 
                               eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                               display = FALSE, ylim = NA)
  dist_diff_l <- dist_A_lim_diff(A_lim = auc_ci_root[1], 
                               D_bar, dt_in = dt_in, metric, num_quantiles, min_eps = 10^(-4), 
                               eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                               display = FALSE, ylim = NA)
  print(dist_diff_u)
  print(dist_diff_l)
  
  # Solve for the upper root, if possible, else leave the starting value.
  if (dist_diff_u > 0) {
    
    # Solve for the upper root.
    ci_u_solve <- uniroot(f = dist_A_lim_diff, 
                          # interval = c(auroc_hat_est, max(auroc_sim_df[, 'auroc'])), 
                          # interval = c(auroc_hat_est, 1), 
                          interval = c(auroc_hat_est, auc_ci_root[2]), 
                          tol = .Machine$double.eps^0.25, maxiter = 100, 
                          D_bar, dt_in = dt_in, metric, num_quantiles, min_eps = 10^(-4), 
                          eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                          display = FALSE, ylim = NA)
    
    auc_ci_root[2] <- ci_u_solve$root
  }
  
  
  # Solve for the lower root, if possible, else leave the starting value.
  if (dist_diff_l > 0) {
    
    # Solve for the lower root.
    ci_l_solve <- uniroot(f = dist_A_lim_diff, 
                          # interval = c(min(auroc_sim_df[, 'auroc']), auroc_hat_est), 
                          # interval = c(0, auroc_hat_est), 
                          interval = c(auc_ci_root[1], auroc_hat_est), 
                          tol = .Machine$double.eps^0.25, maxiter = 100, 
                          D_bar, dt_in = dt_in, metric, num_quantiles, min_eps = 10^(-4), 
                          eta, tol_A_0 = tol, max_iter_A_0 = max_iter, 
                          display = FALSE, ylim = NA)
    
    auc_ci_root[1] <- ci_l_solve$root
  }
  
  
  
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_root', 'cu_u'] <- auc_ci_root[2]
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_root', 'cu_l'] <- auc_ci_root[1]
  
  
  
  #--------------------------------------------------------------------------------
  # Next, run simulation with optimal weights, to match null hypothesis.
  #--------------------------------------------------------------------------------
  
  # Find null distribution.
  adj_weights_dt(dt_in, eta, A_0, tol, max_iter, 
                 display = FALSE, ylim = NA)
  
  
  # Simulation confidence intervals.
  auroc_sim_df <- auroc_sim(dt_in, metric = metric, num_quantiles, num_boots)
  
  print(summary(auroc_sim_df))
  # auroc_sim_df <- auroc_sim_df[order(auroc_sim_df[, 'distance'])]
  
  # Store the corresponding confidence intervals.
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'dist_sim_adj', 'cu_u'] <- quantile(auroc_sim_df[, 'distance'], 1 - alpha)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'dist_sim_adj', 'cu_l'] <- 0
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim_adj', 'cu_u'] <- quantile(auroc_sim_df[, 'auroc'], 1 - alpha/2)
  binorm_ci_df[binorm_ci_df[, 'stat'] == 'auroc_sim_adj', 'cu_l'] <- quantile(auroc_sim_df[, 'auroc'], alpha/2)
  
  
  return(binorm_ci_df)
  
}


#--------------------------------------------------------------------------------
# Simulation of distance metrics and AUROC
#--------------------------------------------------------------------------------

# auroc_sim_df <- auroc_sim(dt_in, metric, num_quantiles, num_boots)


auroc_sim <- function(dt_in, metric, num_quantiles, num_boots) {
  
  # Initialize the data frame of simulation results.
  auroc_sim_df <- data.frame(distance = numeric(num_boots),
                             auroc = numeric(num_boots))
  
  # Get counts of outcomes.
  n_1 <- dt_in[outcomes == TRUE, .N]
  n_0 <- dt_in[outcomes == FALSE, .N]
  
  # Obtain the quantiles of the observed scores by outcome category.
  quantile_df <- distn_quantiles(dt_in, num_quantiles)
  # Obtain a table of relative frequencies for the observed scores.
  distn_table <- quantile_assign(dt_in, quantile_cuts = quantile_df, type = 'observed')
  
  
  
  # print(distn_table)
  
  for (boot_num in 1:num_boots) {
    
    
    # Generate a sample.
    dt_in[outcomes == TRUE, 'score_sim'] <- sample(x = dt_in[outcomes == TRUE, scores], 
                                           size = n_1, replace = TRUE, 
                                           prob = dt_in[outcomes == TRUE, opt_weights])
    dt_in[outcomes == FALSE, 'score_sim'] <- sample(x = dt_in[outcomes == FALSE, scores], 
                                            size = n_0, replace = TRUE, 
                                            prob = dt_in[outcomes == FALSE, opt_weights])
    
    
    # Calculate the AUROC (using the rank approach).
    dt_in[, rank_sim := rank(score_sim)]
    auroc_hat_sim <- dt_in[outcomes == TRUE, ( sum(rank_sim) - 0.5*n_1*(n_1 + 1) )/(n_1*n_0)]
    
    
    # Calculate the distance from the observed distribution.
    
    # Obtain a table of relative frequencies from the simulated scores.
    distn_table[, 'distn_sim'] <- quantile_assign(dt_in, quantile_cuts = quantile_df, 
                                                  type = 'simulted')[, 'distn_sim']
    
    # Pass the tabulated results to calculate distance.
    auroc_dist_sim <- auroc_dist(outcomes = distn_table[, 'outcomes'], 
                                 weights_1 = distn_table[, 'distn_obs'], 
                                 weights_2 = distn_table[, 'distn_sim'], 
                                 metric = metric)
    
    # print(distn_table)
    
    # Store the results in the data frame.
    auroc_sim_df[boot_num, 'auroc'] <- auroc_hat_sim
    auroc_sim_df[boot_num, 'distance'] <- auroc_dist_sim
    
  }
  
  return(auroc_sim_df)
  
}


################################################################################
# Output and Display
################################################################################

#--------------------------------------------------------------------------------
# Joint distance between pairs of positive and negative distributions.
#--------------------------------------------------------------------------------

# Plots two histograms of scores by outcome. 

plot_distn_shift <- function(scores, outcomes, weights, ylim) {
  
  colTransparency <- 1/4
  col_pos <- rgb(0,0,1,colTransparency)
  col_neg <- rgb(1,0,0,colTransparency)
  
  
  weighted.hist(main = 'Shifting of Positive and Negative Distributions',
                x = scores[outcomes], 
                w = weights[outcomes], 
                col = col_pos,
                ylim = ylim,
                xaxis = FALSE) #,
  # add = TRUE)
  weighted.hist(x = scores[!outcomes], 
                w = weights[!outcomes],
                col = col_neg,
                ylim = ylim,
                xaxis = FALSE,
                add = TRUE)
  
}



# Unused variables.
# col_pos_orig <- rgb(0,0,1, 1)
# col_neg_orig <- rgb(1,0,0, 1)

# x_lim <- c(min(binorm_df[, 'scores']), max(binorm_df[, 'scores']))
# y_lim <- c(0, 100)




#--------------------------------------------------------------------------------


################################################################################
# Distance Metrics to the null hypothesized AUROC
################################################################################

#--------------------------------------------------------------------------------
# Calculate relative frequencies of scores for a given set of quantile cuts.
#--------------------------------------------------------------------------------


# distn_table <- quantile_assign(dt_in, quantile_cuts = quantile_df, type = c('observed', 'simulated'))

quantile_assign <- function(dt_in, quantile_cuts, type = c('observed', 'simulated')) {
  
  # Get parameters.
  num_quantiles <- nrow(quantile_cuts) - 1
  n_1 <- dt_in[outcomes == TRUE, .N]
  n_0 <- dt_in[outcomes == FALSE, .N]
  
  # Initialize data frame.
  distn_table <- data.frame(quantiles = c(sprintf('Q_1_%d', 1:num_quantiles), 
                                          sprintf('Q_0_%d', 1:num_quantiles)), 
                            outcomes = c(rep(TRUE, num_quantiles), 
                                         rep(FALSE, num_quantiles)))
  
  if (type == 'observed') {
    
    dt_in[outcomes == TRUE , 
          quantile := cut(scores,
                          breaks = quantile_cuts[, 'pos_quantile'],
                          include.lowest = TRUE, labels = 1:num_quantiles) ]
    dt_in[outcomes == FALSE , 
          quantile := cut(scores,
                          breaks = quantile_cuts[, 'neg_quantile'],
                          include.lowest = TRUE, labels = 1:num_quantiles) ]
    
    distn_table[, 'distn_obs'] <- c(as.numeric(table(dt_in[outcomes == TRUE, 
                                                           quantile]) / n_1)[1:num_quantiles],
                                    as.numeric(table(dt_in[outcomes == FALSE, 
                                                           quantile]) / n_0)[1:num_quantiles])
    
  } else if (type == 'simulated') {
    
    dt_in[outcomes == TRUE , 
          quantile_sim := cut(score_sim,
                              breaks = quantile_cuts[, 'pos_quantile'],
                              include.lowest = TRUE, labels = 1:num_quantiles) ]
    dt_in[outcomes == FALSE , 
          quantile_sim := cut(score_sim,
                              breaks = quantile_cuts[, 'neg_quantile'],
                              include.lowest = TRUE, labels = 1:num_quantiles) ]
    
    distn_table[, 'distn_sim'] <- c(as.numeric(table(dt_in[outcomes == TRUE, 
                                                           quantile_sim]) / n_1)[1:num_quantiles],
                                    as.numeric(table(dt_in[outcomes == FALSE, 
                                                           quantile_sim]) / n_0)[1:num_quantiles])
    
    
  } # Else no third column.
  
  return(distn_table)
  
}






#--------------------------------------------------------------------------------
# Joint distance between pairs of positive and negative distributions.
# Binned version - sorted into quantiles.
#--------------------------------------------------------------------------------

# quantile_df <- distn_quantiles(dt_in, num_quantiles)
# num_quantiles <- 10

distn_quantiles <- function(dt_in, num_quantiles) {
  
  # Initialize by outcome category.
  quantile_df <- data.frame(pos_quantile = numeric(num_quantiles + 1),
                             neg_quantile = numeric(num_quantiles + 1))
  
  # Cut scores in each outcome class into quantiles.
  pos_quantile <- quantile(dt_in[outcomes == TRUE, scores], 
                           seq(0, 1, by = 1/num_quantiles))
  neg_quantile <- quantile(dt_in[outcomes == FALSE, scores], 
                           seq(0, 1, by = 1/num_quantiles))
  
  # Store quantiles for output.
  quantile_df[, 'pos_quantile'] <- pos_quantile
  quantile_df[, 'neg_quantile'] <- neg_quantile
  
  
  return(quantile_df)
  
}


#--------------------------------------------------------------------------------
# Joint distance between pairs of positive and negative distributions.
#--------------------------------------------------------------------------------

# dist <- auroc_dist(outcomes, weights_1, weights_2, metric, min_eps = 10^(-4))

# f_1 and f_2 are probability vectors of equal length for positive observations.
# g_1 and g_2 are probability vectors of equal length for negative observations.

auroc_dist <- function(outcomes, weights_1, weights_2, metric, min_eps = 10^(-4)) {
  
  # Parse out pairs of distributions.
  f_1 <- weights_1[outcomes]
  g_1 <- weights_1[!outcomes]
  f_2 <- weights_2[outcomes]
  g_2 <- weights_2[!outcomes]
  
  # Modify for minimum distance when necessary.
  if (metric %in% c('kld0_joint', 'kld0', 'kld1')) {
    
    f_1 <- pmax(f_1, min_eps)
    g_1 <- pmax(g_1, min_eps)
    f_2 <- pmax(f_2, min_eps)
    g_2 <- pmax(g_2, min_eps)
    
  }
  
  # Calculate distance corresponding to desired metric.
  if (tolower(metric) == 'chi2_joint') {
    
    dist <- 0
    for (i in 1:length(f_1)) {
      dist <- dist + chi2_distance(f_1[i]*g_1, f_2[i]*g_2)
    }
    
  } else if (tolower(metric) == 'kld0_joint') {
    
    dist <- 0
    for (i in 1:length(f_1)) {
      dist <- dist + kld_0_distance(f_1[i]*g_1, f_2[i]*g_2)
    }
    
  } else if (tolower(metric) == 'chi2') {
    
    dist <- chi2_distance(f_1, f_2) + chi2_distance(g_1, g_2)
    
  } else if (tolower(metric) == 'kld0') {
    
    dist <- kld_0_distance(f_1, f_2) + kld_0_distance(g_1, g_2)
    
  } else if (tolower(metric) == 'kld1') {
    
    dist <- kld_1_distance(f_1, f_2) + kld_1_distance(g_1, g_2)
    
    
  } else {
    
    print('Warning: Metric not recognized.')
    dist <- NA
  }
  
  return(dist)
  
}



#--------------------------------------------------------------------------------
# Chi-squared Distance
#--------------------------------------------------------------------------------

# distn_1 and distn_2 are probability vectors of equal length.
# dist is the scalar distance between these two distributions.

chi2_distance <- function(distn_1, distn_2) {
  
  n_obs <- length(distn_1)
  
  dist <- n_obs * sum((distn_1 - distn_2)^2)
  
}

#--------------------------------------------------------------------------------
# Kullback Leibler Divergence - Standard Definition
#--------------------------------------------------------------------------------

# distn_1 and distn_2 are probability vectors of equal length.
# distn_1 is the benchmark distribution.
# dist is the scalar distance between these two distributions.

kld_0_distance <- function(distn_1, distn_2) {
  
  n_obs <- length(distn_1)
  
  dist <- n_obs * sum(distn_1 * (log(distn_1) - log(distn_2)))
  
}

#--------------------------------------------------------------------------------
# Kullback Leibler Divergence - Symmetric Distance Metric
#--------------------------------------------------------------------------------

# distn_1 and distn_2 are probability vectors of equal length.
# dist is the scalar distance between these two distributions.

kld_1_distance <- function(distn_1, distn_2) {
  
  n_obs <- length(distn_1)
  
  dist <- n_obs * sum((distn_1 - distn_2) * (log(distn_1) - log(distn_2)))
  
}


################################################################################
# Imposing the null hypothesized AUROC
################################################################################

#--------------------------------------------------------------------------------
# Distance as a function of a suggested A_0
#--------------------------------------------------------------------------------

# dist_diff <- dist_A_lim_diff(A_lim, D_bar, dt_in, metric, min_eps = 10^(-4), 
#                              eta, tol_A_0, max_iter_A_0, 
#                              display = FALSE, ylim = NA)

dist_A_lim_diff <- function(A_lim, D_bar, dt_in, metric, num_quantiles = FALSE, min_eps = 10^(-4), 
                            eta, tol_A_0, max_iter_A_0, 
                            display = FALSE, ylim = NA) {
  
  # Calculate optimal weights to achieve the suggested A_lim.
  adj_weights_dt(dt_in, eta, A_0 = A_lim, tol = tol_A_0, max_iter = max_iter_A_0, 
                 display = FALSE, ylim = NA)
  
  # Calculate distance.
  if (num_quantiles == FALSE) {
    
    # Empirical from theoretical distance.
    distance_A_lim <- auroc_dist(outcomes = dt_in[, outcomes], 
                                 weights_1 = dt_in[, weights], 
                                 weights_2 = dt_in[, opt_weights], 
                                 metric, min_eps)
    
  } else {
    
    # Binned distribution for comparison with empirical distributions (apples to apples).
    
    # Calculate distribution for actual observations (uniform).
    quantile_df <- distn_quantiles(dt_in, num_quantiles)
    distn_table <- quantile_assign(dt_in, quantile_cuts = quantile_df, type = 'observed')
    
    # Calculate distribution for optimized observations.
    distn_table[distn_table[, 'outcomes'] == TRUE, 'distn_sim'] <- 
      binorm_dt[outcomes == TRUE, sum(opt_weights), by = quantile][, V1]
    distn_table[distn_table[, 'outcomes'] == FALSE, 'distn_sim'] <- 
      binorm_dt[outcomes == FALSE, sum(opt_weights), by = quantile][, V1]
    
    # Calculate distance to "true" distribution, as in the closest distribution
    # for which the null is satisfied.
    distance_A_lim <- auroc_dist(outcomes = distn_table[, 'outcomes'], 
                                 weights_1 = distn_table[, 'distn_obs'], 
                                 weights_2 = distn_table[, 'distn_sim'], 
                                 metric, min_eps)
  }
  
  # Calculate difference of current distance from that desired. 
  dist_diff <- distance_A_lim - D_bar
  
  # print(distance_A_lim)
  # print(A_lim)
  
  return(dist_diff)
  
}




#--------------------------------------------------------------------------------
# Solving for the limiting A_0 a set distance away.
# See above, using a root-solving algorithm.
# Need only define the distance objective function.
#--------------------------------------------------------------------------------

# A_lim <- find_A_lim(dt_in, eta, D_bar, tol, max_iter)

find_A_lim <- function(dt_in, eta, D_bar, tol, max_iter) {
  
  # Initialize current distance and weights.
  D_current <- 0
  A_hat_current <- auc_weighted_loops(scores = dt_in[, scores], 
                                      outcomes = dt_in[, outcomes], 
                                      weights = dt_in[, weights])
  dt_in[, opt_weights := weights]
  num_iter <- 0
  
  # Search toward further distance.
  while (abs(D_current - D_bar) > tol & num_iter <= max_iter) {
    
    num_iter <- num_iter + 1
    
    # Adjust weights.
    iter_weights_dt(dt_in, eta, A_0, A_hat_current, 
                    display = display, ylim = ylim)
    
  }
  
  return(A_lim)
  
}

#--------------------------------------------------------------------------------
# AUC Shift - Move to nearest distribution consistent with the null hypothesis.
# Data table edition
#--------------------------------------------------------------------------------

# Iteratively adjusts weights untill distribution is the nearest distribution 
# that is consistent with the null hypothesis.

# void <- adj_weights_dt(dt_in, eta, A_0, tol, max_iter, 
#                        display = FALSE, ylim = NA)

adj_weights_dt <- function(dt_in, eta, A_0, tol, max_iter, 
                        display = FALSE, ylim = NA) {
  
  # print(summary(dt_in))
  
  # Initialize current state with estimated AUROC and associated weights.
  # A_hat_current <- auc_weighted_loops_dt(dt_in) # Stoopidly inefficient.
  A_hat_current <- auc_weighted_loops(scores = dt_in[, scores], 
                                      outcomes = dt_in[, outcomes], 
                                      weights = dt_in[, weights])
  
  # Initialize optimized weights and iterate.
  dt_in[, opt_weights := weights]
  num_iter <- 0
  while (abs(A_hat_current - A_0) > tol & num_iter <= max_iter) {
    
    num_iter <- num_iter + 1
    
    # Adjust weights.
    iter_weights_dt(dt_in, eta, A_0, A_hat_current, 
                    display = display, ylim = ylim)
    
    # Recalculate AUROC estimate under new weights.
    # A_hat_current <- auc_weighted_loops_dt(dt_in) # Stoopidly inefficient.
    A_hat_current <- auc_weighted_loops(scores = dt_in[, scores], 
                                        outcomes = dt_in[, outcomes], 
                                        weights = dt_in[, opt_weights])
    
    # print(A_hat_current)
    
  }
  
  # print(summary(dt_in))
  
  # return(dt_in[, opt_weights])
  
}


#--------------------------------------------------------------------------------
# AUC Shift - Move to nearest distribution consistent with the null hypothesis.
#--------------------------------------------------------------------------------

# Iteratively adjusts weights untill distribution is the nearest distribution 
# that is consistent with the null hypothesis.

# new_weights <- adj_weights(scores, outcomes, weights, eta, A_0, A_hat_current)

adj_weights <- function(scores, outcomes, weights, eta, A_0, A_hat, tol, max_iter, 
                        display = FALSE, ylim = NA) {
  
  # Initialize current state with estimated AUROC and associated weights.
  A_hat_current <- A_hat # Could remove A_hat and calculate inside.
  new_weights <- weights
  num_iter <- 0
  while (abs(A_hat_current - A_0) > tol & num_iter <= max_iter) {
    
    num_iter <- num_iter + 1
    
    # Adjust weights.
    new_weights <- iter_weights(scores, outcomes, new_weights, 
                                eta, A_0, A_hat_current, 
                                display = display, ylim = ylim)
    
    # Recalculate AUROC estimate under new weights.
    # A_hat_current <- auc_weighted(scores, outcomes, weights)
    A_hat_current <- auc_weighted_loops(scores, outcomes, new_weights)
    
    print(A_hat_current)
    
  }
  
  return(new_weights)
  
}


#--------------------------------------------------------------------------------
# AUC Shift - Iteration on distributions (Data Table version)
#--------------------------------------------------------------------------------

# void <- iter_weights_dt(dt_in, eta, A_0, A_hat, 
#                         display = FALSE, ylim = NA)

iter_weights_dt <- function(dt_in, eta, A_0, A_hat, 
                            display = FALSE, ylim = NA) {
  
  # Calculate step size.
  lambda <- eta*(A_hat - A_0) # The last two are equivalent.
  
  # Calculate gradient in weights.
  # Note that this assumes that data are sorted by scores (increasing order).
  foc_gradient_dt(dt_in)
  
  # Calculate information criteria (note opposit outcomes).
  pos_exp_kld <- exp(sum(dt_in[outcomes == FALSE, opt_weights*log(opt_weights)]))
  neg_exp_kld <- exp(sum(dt_in[outcomes == TRUE, opt_weights*log(opt_weights)]))
  
  # Calculate updated weights.
  dt_in[outcomes == TRUE, opt_weights := opt_weights * pos_exp_kld * exp(lambda * foc_gradient)]
  dt_in[outcomes == FALSE, opt_weights := opt_weights * neg_exp_kld * exp(lambda * foc_gradient)]
  
  # Normalize for unit probability mass.
  dt_in[, opt_weights := opt_weights/sum(opt_weights), by = outcomes]
  
  if (display) {
    
    plot_distn_shift(scores = dt_in[, scores], 
                     outcomes = dt_in[, outcomes], 
                     weights = dt_in[, opt_weights], 
                     ylim = ylim)
    
    # print(summary(dt_in[outcomes, opt_weights]))
    # print(summary(dt_in[!outcomes, opt_weights]))
    
  }
  
  # # Calculate updated weights.
  # dt_in[, 'opt_weights'] <- NA_real_
  # # Positive weights.
  # # pos_wts_next <- pos_wts*exp(sum(neg_wts*log(neg_wts)) + lambda*pos_gradient)
  # set(dt_in, 
  #     i = which(dt_in[, get('outcomes')]), 
  #     j = 'opt_weights', 
  #     value = dt_in[outcomes == TRUE, get('weights')] * 
  #       exp(sum(dt_in[outcomes == FALSE, weights*log(weights)])) * 
  #       exp(lambda * dt_in[outcomes == TRUE, get('foc_gradient')]))
  # 
  # 
  # # Negative weights.
  # # neg_wts_next <- neg_wts*exp(sum(pos_wts*log(pos_wts)) + lambda*neg_gradient)
  # set(dt_in, 
  #     i = which(!dt_in[, get('outcomes')]), 
  #     j = 'opt_weights', 
  #     value = dt_in[outcomes == FALSE, get('weights')] * 
  #       exp(sum(dt_in[outcomes == TRUE, weights*log(weights)])) * 
  #       exp(lambda * dt_in[outcomes == FALSE, get('foc_gradient')]))
  
  
  # # Normalize for unit probability mass.
  # # pos_wts_next <- pos_wts_next/sum(pos_wts_next)
  # set(dt_in, 
  #     i = which(dt_in[, get('outcomes')]), 
  #     j = 'opt_weights', 
  #     value = dt_in[outcomes == TRUE, get('opt_weights')] / 
  #       sum(dt_in[outcomes == TRUE, get('opt_weights')]))
  # # neg_wts_next <- neg_wts_next/sum(neg_wts_next)
  # set(dt_in, 
  #     i = which(!dt_in[, get('outcomes')]), 
  #     j = 'opt_weights', 
  #     value = dt_in[outcomes == FALSE, get('opt_weights')] / 
  #       sum(dt_in[outcomes == FALSE, get('opt_weights')]))
  
  
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
iter_weights <- function(scores, outcomes, weights, eta, A_0, A_hat, 
                         display = FALSE, ylim = NA) {
  
  # Calculate step size.
  # lambda <- eta*(A_0 - A_hat)
  # lambda <- - eta*(A_0 - A_hat)
  # A_hat <- A_hat_curr_binom
  lambda <- eta*(A_hat - A_0) # The last two are equivalent.
  
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
  # Depend on particular distance metric.
  # pos_wts_next <- pos_wts*exp(lambda*pos_gradient)
  # pos_wts_next <- exp(sum(neg_wts*log(neg_wts)) + lambda*pos_gradient)
  pos_wts_next <- pos_wts*exp(sum(neg_wts*log(neg_wts)) + lambda*pos_gradient)
  # pos_wts_next <- pos_wts_iter(neg_wts, lambda, pos_gradient) # Make function later.
  # neg_wts_next <- neg_wts*exp(lambda*neg_gradient)
  # neg_wts_next <- exp(sum(pos_wts*log(pos_wts)) + lambda*neg_gradient)
  neg_wts_next <- neg_wts*exp(sum(pos_wts*log(pos_wts)) + lambda*neg_gradient)
  # neg_wts_next <- neg_wts_iter(pos_wts, lambda, neg_gradient, metric = 'kld') # Make function later.
  
  # Normalize for unit probability mass.
  pos_wts_next <- pos_wts_next/sum(pos_wts_next)
  neg_wts_next <- neg_wts_next/sum(neg_wts_next)
  
  # Place these weights in the data frame.
  next_weights <- 0*weights
  next_weights[outcomes] <- pos_wts_next
  next_weights[!outcomes] <- neg_wts_next
  
  # Check fixed ponits and progress of distribution shifting.
  if (display == TRUE) {
    pos_kld_0_dist <- kld_0_distance(distn_1 = pos_wts, distn_2 = pos_wts_next)
    neg_kld_0_dist <- kld_0_distance(distn_1 = neg_wts, distn_2 = neg_wts_next)
    
    print(sprintf('Positive weight change: %f', pos_kld_0_dist))
    print(sprintf('Negative weight change: %f', neg_kld_0_dist))
    
    
    plot_distn_shift(scores, outcomes, next_weights, 
                     ylim = ylim)
    
    # Compare FOC to check of conditions satisfied.
    m <- sum(!outcomes)
    n <- sum(outcomes)
    # Set up regression model. 
    curr_gradient <- NA*weights
    curr_gradient[outcomes] <- pos_gradient
    curr_gradient[!outcomes] <- neg_gradient
    curr_gradient <- curr_gradient
    foc_lm <- lm(formula = log_prob ~ curr_gradient, 
                 data = data.frame(log_prob = log(next_weights),
                                   curr_gradient = curr_gradient))
    lambda_hat <- coef(foc_lm)['curr_gradient']
    print(sprintf('Current Lagrange multiplier: %f', lambda_hat))
    
    pos_wt_FOC <- exp(sum(neg_wts*(log(neg_wts) - log(1/m)))  + 
                             lambda_hat*pos_gradient) * rep(1, n) 
    # summary(pos_wt_FOC)
    # summary(pos_wt_FOC/sum(pos_wt_FOC))
    # summary(pos_wts_next)
    # summary(pos_wt_FOC/sum(pos_wt_FOC) - 1/n)
    avg_pos_wt_FOC_diff <- sum(abs(pos_wts_next - pos_wt_FOC/sum(pos_wt_FOC)))/n
    
    neg_wt_FOC <- exp(sum(pos_wts*(log(pos_wts) - log(1/n)))  + 
                             lambda_hat*neg_gradient) * rep(1, m) 
    # summary(neg_wt_FOC)
    # summary(neg_wt_FOC/sum(neg_wt_FOC))
    # summary(neg_wts_next)
    avg_neg_wt_FOC_diff <- sum(abs(neg_wts_next - neg_wt_FOC/sum(neg_wt_FOC)))/m
    
    print(sprintf('Average positive weight FOC difference: %f', avg_pos_wt_FOC_diff))
    print(sprintf('Average negative weight FOC difference: %f', avg_neg_wt_FOC_diff))
    
    print(summary(pos_wts_next))
    print(summary(neg_wts_next))
    
    # Sys.sleep(0)
    
  }
  
  return(next_weights)
  
}


#--------------------------------------------------------------------------------
# AUC Shift - Positive and Negative conditional distribution
#--------------------------------------------------------------------------------

# Computes a reweighting gradient to adjust the distribution of the scores 
# corresponding to both positive and negative observations.
# Data table version.

# void <- foc_gradient_dt(dt_in)

foc_gradient_dt <- function(dt_in) {
  
  # Calculate vectors of selected weights.
  dt_in[, z_x_w := outcomes * opt_weights]
  dt_in[, zz_x_w := (1 - outcomes) * opt_weights]
  
  # Calculate partial sums of selected weights.
  dt_in[, sum_z_x_w := cumsum(z_x_w)]
  dt_in[, sum_zz_x_w := 1 - cumsum(zz_x_w)]
  
  # Assign values according to outcomes.
  dt_in[outcomes == FALSE, foc_gradient := sum_z_x_w]
  dt_in[outcomes == TRUE, foc_gradient := sum_zz_x_w]
  
  
  # Original data frame option.
  # neg_gradient <- cumsum(outcomes * weights)[!outcomes]
  # Data table version.
  # dt_in[, 'foc_gradient'] <- NA_real_
  # set(dt_in, 
  #     i = which(!dt_in[, get('outcomes')]), 
  #     j = 'foc_gradient', 
  #     value = dt_in[, sum_z_x_w][!dt_in[, get('outcomes')]])
  # value = dt_in[, cumsum(z_x_w)][!dt_in[, get('outcomes')]])
  # dt_in[!outcomes, 'foc_gradient'] <- dt_in[, cumsum(z_x_w)][!dt_in[, get('outcomes')]]
  
  # Original data frame option.
  # pos_gradient <- (1 - cumsum((1 - outcomes) * weights))[outcomes]
  # Data table version.
  # set(dt_in, 
  #     i = which(dt_in[, get('outcomes')]), 
  #     j = 'foc_gradient', 
  #     value = (1 - dt_in[, cumsum(zz_x_w)])[dt_in[, get('outcomes')]])
  # value = (1 - dt_in[, sum_zz_x_w])[dt_in[, get('outcomes')]])
  # dt_in[outcomes, 'foc_gradient'] <- dt_in[, cumsum(zz_x_w)][dt_in[, get('outcomes')]]
  
}



#--------------------------------------------------------------------------------
# AUC Shift - Positive conditional distribution
#--------------------------------------------------------------------------------

# Computes a reweighting gradient to adjust the distribution of the scores 
# corresponding to positive observations.

# pos_gradient <- neg_cond_shift(scores, outcomes, weights)

# Vector version first.
neg_cond_shift <- function(scores, outcomes, weights) {
  
  # pos_gradient <- cumsum(!outcomes[order(-scores)] * weights[order(-scores)])[outcomes[order(-scores)]]
  
  # Need to reverse order, since sorting by score is in reverse.
  # pos_gradient <- pos_gradient[order(seq(length(pos_gradient), 1, by = -1))]
  
  
  # Compare to opposite gradient.
  # neg_gradient <- cumsum(outcomes[order(scores)] * weights[order(scores)])[!outcomes[order(scores)]]
  # Assume already sorted.
  # neg_gradient <- cumsum(outcomes * weights)
  # neg_gradient <- cumsum(outcomes * weights)[!outcomes]
  # neg_gradient <- neg_gradient/sum(outcomes)/sum(!outcomes)
  
  # Current gradient.
  # pos_gradient <- cumsum(!outcomes[order(scores)] * weights[order(scores)])[outcomes[order(scores)]]
  # Assume already sorted.
  # pos_gradient <- cumsum(!outcomes * weights)
  # Got it! !outcomes coerces to logical (need brackets or 1 - outcomes: ! operates last).
  # pos_gradient <- 1 - cumsum((1 - outcomes) * weights)
  pos_gradient <- (1 - cumsum((1 - outcomes) * weights))[outcomes]
  
  # Normalize.
  # pos_gradient <- pos_gradient/sum(outcomes)
  # pos_gradient <- pos_gradient/sum(outcomes)/sum(!outcomes)
  
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
  
  # Complicated way without previously sorted.
  # neg_gradient <- cumsum(outcomes[order(scores)] * weights[order(scores)])[!outcomes[order(scores)]]
  
  # Assume already sorted.
  neg_gradient <- cumsum(outcomes * weights)[!outcomes]
  
  # Normalize.
  # neg_gradient <- neg_gradient/sum(!outcomes)
  # neg_gradient <- neg_gradient/sum(outcomes)/sum(!outcomes)
  
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








################################################################################
# End
################################################################################
