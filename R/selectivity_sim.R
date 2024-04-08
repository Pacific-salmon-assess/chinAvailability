## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)
library(mvtweedie)

# modified tweedie prediction that allows for exclude = list()
source(here::here("R", "functions", "pred_mvtweedie2.R"))


# import SRKW prey data
# rkw_dat <- readRDS(
#   here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
# ) %>% 
#   # downscale utm coordinates to match model
#   mutate(
#     utm_y = utm_y / 1000,
#     utm_x = utm_x / 1000
#     ) %>% 
#   rename(
#     week_n = week,
#     year_n = year
#   ) %>% 
#   filter(
#     era == "current"
#   ) 

rkw_dat_pooled <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat_pooled.rds")
) %>% 
  # downscale utm coordinates to match model
  mutate(
    utm_y = utm_y / 1000,
    utm_x = utm_x / 1000
  ) %>% 
  rename(
    week_n = week
  ) %>% 
  filter(
    era == "current"
  ) %>% 
  arrange(-n_samples)


# import original fishery data used to fit model
agg_dat <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
)


# import fitted model
model_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
  )
)


# split diet data by sample ID and simulate
rkw_dat_tbl <- rkw_dat_pooled %>% 
  group_by(sample_id_pooled) %>% 
  group_nest() 


# make new data frame consistent with SRKW diet to generate predictions from
new_dat_list <- purrr::map(
  rkw_dat_tbl$data,
  function (x) {
    newdata <- expand.grid(
      sample_size = unique(x$n_samples),
      week_n = unique(x$week_n),
      stock_group = levels(model_fit$model$stock_group),
      utm_x = unique(x$utm_x),
      utm_y = unique(x$utm_y),
      year_n = "2020" # dummy year so value doesn't matter
    ) %>% 
      left_join(
        ., 
        x %>% 
          select(stock_group, obs_prob = agg_prob) %>% 
          distinct(),
        by = c("stock_group")
      ) %>% 
      mutate(
        obs_prob = ifelse(is.na(obs_prob), 0, obs_prob)
      )
    
    # exclude year smoother and generate predictions based on fishery model
    excl <- grepl("year_n", gratia::smooths(model_fit))
    yr_coefs <- gratia::smooths(model_fit)[excl]
    preds <- pred_dummy(
      model_fit,
      se.fit = TRUE,
      category_name = "stock_group",
      origdata = agg_dat,
      newdata = newdata,
      exclude = yr_coefs
    )
    
    cbind(
      newdata, pred_prob = preds$fit, se_pred = preds$se.fit 
    ) %>% 
      distinct() 
  }
)


# Function to calculate test statistic (e.g., difference in proportions)
calculate_test_statistic <- function (sample1, sample2) {
  # Calculate proportions
  prop1 <- sample1 / sum(sample1)
  prop2 <- sample2 / sum(sample2)
  
  # Calculate difference in proportions
  diff_prop <- abs(prop1 - prop2)
  
  # Sum of differences
  sum_diff <- sum(diff_prop)
  
  return(sum_diff)
}


# simulate based on each observation dataset 
sim_list <- purrr::map(
  new_dat_list,
  function (x, num_simulations = 2000) {
    num_trials <- x$sample_size  # Sample size for each multinomial draw
    proportions1 <- x$pred_prob  # Proportions for first vector
    proportions2 <- x$obs_prob  # Proportions for second vector
    
    # vector to store results of simulation
    obs_statistics <- test_statistics <- numeric(num_simulations)
    # diff_null <- diff_obs <- diff_obs2 <- matrix(NA, nrow = num_simulations, 
    #                                              ncol = length(proportions1))
    for (i in 1:num_simulations) {
      # Generate multinomial samples under null hypothesis 
      # (i.e. same distribution)
      sample1_null <- rmultinom(1, num_trials, proportions1)
      sample2_null <- rmultinom(1, num_trials, proportions1)
      
      # Calculate test statistic for null samples
      test_statistics[i] <- calculate_test_statistic(sample1_null, sample2_null)
      
      # Generate multinomial samples from observations and calc test statistic
      sample_obs <- rmultinom(1, num_trials, proportions2)
      obs_statistics[i] <- calculate_test_statistic(sample1_null, sample_obs)
      # diff_null[i, ] <- abs(sample2_null - sample1_null)
      # diff_obs[i, ] <- abs(sample_obs - sample1_null)
    }
    
    # calculate p-value (i.e. proportion of sims where null differnece greater 
    # than or equal to observed difference)
    sum(test_statistics >= obs_statistics) / num_simulations
  }
)


p_dat <- rkw_dat_tbl %>%
  mutate(
    p_value = sim_list %>% unlist()
  ) %>% 
  unnest(cols = data) %>%
  select(sample_id_pooled, month, strata, n_samples, p_value) %>% 
  distinct()

# color pal used in sampling maps
strata_colour_pal <- RColorBrewer::brewer.pal(
  length(levels(agg_dat$strata)),
  "Paired"
)
names(strata_colour_pal) <- levels(agg_dat$strata) 

p_plot <- ggplot(p_dat) +
  geom_point(
    aes(x = n_samples, y = p_value, fill = strata), shape = 21
  ) +
  geom_hline(
    aes(yintercept = 0.05), lty = 2 
  ) +
  labs(
    y = "Selectivity p-value", x = "Number of Prey Samples"
  ) +
  scale_fill_manual(values = strata_colour_pal) +
  ggsidekick::theme_sleek()

png(
  here::here("figs", "ms_figs", "selectivity_pvalue.png"),
  height = 4.5, width = 5.5, units = "in", res = 250
)
p_plot
dev.off()