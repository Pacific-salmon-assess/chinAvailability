### Simulation Test
## March 14, 2024
## Example simulation whereby recent sampling events are compared to background 
## predictions to evaluate evidence for selectivity


# Define parameters
num_trials <- 2  # Sample size for each multinomial distribution
proportions1 <- c(0.25, 0.4, 0.3, 0.05)  # Proportions for first vector
# proportions2 <- c(0.3, 0.4, 0.3)  # Proportions for first vector
proportions2 <- c(0.9, 0.05, 0, 0.5)  # Proportions for second vector
num_simulations <- 1000  # Number of simulations

# Function to calculate test statistic (e.g., difference in proportions)
calculate_test_statistic <- function(sample1, sample2) {
  # Calculate proportions
  prop1 <- sample1 / sum(sample1)
  prop2 <- sample2 / sum(sample2)
  
  # Calculate difference in proportions
  diff_prop <- abs(prop1 - prop2)
  
  # Sum of differences
  sum_diff <- sum(diff_prop)
  
  return(sum_diff)
}

# Run simulations
obs_statistics <- test_statistics <- numeric(num_simulations)
diff_null <- diff_obs <- diff_obs2 <- matrix(NA, nrow = num_simulations, 
                                             ncol = length(proportions1))
for (i in 1:num_simulations) {
  # Generate multinomial samples under null hypothesis (i.e. same distribution)
  sample1_null <- rmultinom(1, num_trials, proportions1)
  sample2_null <- rmultinom(1, num_trials, proportions1)
  
  # Calculate test statistic for null samples
  test_statistics[i] <- calculate_test_statistic(sample1_null, sample2_null)
  
  
  sample_obs <- rmultinom(1, num_trials, proportions2)
  
  obs_statistics[i] <- calculate_test_statistic(sample1_null, sample_obs)
  
  diff_null[i, ] <- abs(sample2_null - sample1_null)
  diff_obs[i, ] <- abs(sample_obs - sample1_null)
  diff_obs2[i, ] <- sample_obs > sample1_null
}

sum(test_statistics >= obs_statistics) / num_simulations

diff_mat <- diff_null > diff_obs
apply(diff_mat, 2, mean)
apply(diff_mat, 2, sd)
