### Simulation Test
## March 14, 2024
## Example simulation whereby recent sampling events are compared to background 
## predictions to evaluate evidence for selectivity
## Specifically what is the probability that the most common stock in the diet
## is the most common stock in the environment

n_stocks <- 3
obs_ppn <- c(0.6, 0.3, 0.1)
true_ppn <- c(0.4, 0.5, 0.1)
sample_size <- 10
ntrials <- 500


selectivity_foo <- function (
    n_stocks = 3,
    obs_ppn = c(0.6, 0.3, 0.1),
    true_ppn = c(0.4, 0.5, 0.1),
    sample_size = 10,
    ntrials = 500
) {
  success1 <- rep(NA, ntrials)
  success2 <- rep(NA, ntrials)
  
  if (length(obs_ppn) != length(true_ppn)) {
    stop("Supply vectors of equal length")
  }
  
  for (i in 1:ntrials) {
    true_count <- rmultinom(sample_size, 1, true_ppn)
    sample_ppn_true <- apply(true_count, 1, sum) / sample_size
    
    # obs_count <- rmultinom(sample_size, 1, obs_ppn)
    # sample_ppn_obs <- apply(obs_count, 1, sum) / sample_size
    
    # is ppn of dominant stock >/= to observed
    success1[i] <- if (
      sample_ppn_true[which.max(obs_ppn)] >= sample_ppn_obs[which.max(obs_ppn)]
    ) {
      1
    } else {
      0
    }
    
    # is dominant stock the same 
    success2[i] <- if (which.max(sample_ppn_true) == which.max(obs_ppn)) {
      1
    } else {
      0
    }
  }
  list(
    prob_success1 = sum(success1) / ntrials, 
    prob_success2 = sum(success2) / ntrials
  )
}

dd <- selectivity_foo(
  n_stocks = 3,
  obs_ppn = c(0.6, 0.3, 0.1),
  true_ppn = c(0.6, 0.3, 0.1),
  sample_size = 10,
  ntrials = 500
)

samps_list <- purrr::map(
  list(5, 10, 50),
  ~ selectivity_foo(
    n_stocks = 3,
    obs_ppn = c(0.6, 0.3, 0.1),
    true_ppn = c(0.5, 0.2, 0.1),
    sample_size = .x
  )
)

purrr::map(samps_list, ~ .x$prob_success1)
purrr::map(samps_list, ~ .x$prob_success2)
