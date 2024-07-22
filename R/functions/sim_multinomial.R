## Function to simulate composition estimates from multinomial based on DF of 
# predictions

sim_foo <- function(pred_dat_in, nsim = 50) {
  sim_dat_out <- vector(mode = "list", length = nsim)
  
  for (i in 1:nsim) {
    # generate vector of random variables based on predicted proportions and SEs
    r_vec <- rep(NA, length = nrow(pred_dat_in))
    for (j in seq_along(r_vec)) {
      r_vec[j] <- rnorm(1, mean = pred_dat_in$fit[j], sd = pred_dat_in$se[j])
    }
    # remove negative random draws
    r_vec2 <- pmax(r_vec, 0)
    pred_dat_in$r_vec <- r_vec2
    
    # rescale so that proportions sum to one then split into a list
    scaled_pred_list <- pred_dat_in %>% 
      group_by(sample_id) %>% 
      mutate(
        new_ppn = r_vec / sum(r_vec)
      ) %>% 
      ungroup() %>% 
      split(., .$sample_id)
    
    # generate null multinomial distribution for each sampling event
    sim_dat <- purrr::map(
      scaled_pred_list,
      function (x) {
        # Generate multinomial samples under null hypothesis (i.e. same dist)
        x %>% 
          mutate(
            sample_size = unique(x$sample_size),
            sample1_null = rmultinom(1, sample_size, x$new_ppn) %>% 
              as.numeric(),
            sample2_null = rmultinom(1, sample_size, x$new_ppn) %>% 
              as.numeric(),
            sample_obs = rmultinom(1, sample_size, x$obs_ppn) %>% 
              as.numeric()
          )
      } 
    ) %>% 
      bind_rows()
    
    # for simulation calculate total counts then test whether obs differences
    # greater than null differences
    sim_dat_agg <- sim_dat %>% 
      group_by(
        stock_group
      ) %>% 
      summarize(
        sample1_count = sum(sample1_null),
        sample2_count = sum(sample2_null),
        obs_count = sum(sample_obs),
        null_diff = abs(sample1_count - sample2_count),
        obs_diff = abs(sample1_count - obs_count),
        test_stat = ifelse(null_diff >= obs_diff, 1, 0)
      ) %>% 
      mutate(
        sim_i = i %>% as.character()
      )
    sim_dat_out[[i]] <- sim_dat_agg
  }
  
  sim_dat_out %>% 
    bind_rows()
}