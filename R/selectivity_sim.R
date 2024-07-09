## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)
library(mvtweedie)

# modified tweedie prediction that allows for exclude = list()
source(here::here("R", "functions", "pred_mvtweedie2.R"))


# Function to calculate test statistic (e.g., difference in proportions)
# calculate_test_statistic <- function (sample1, sample2) {
#   # Calculate proportions
#   prop1 <- sample1 / sum(sample1)
#   prop2 <- sample2 / sum(sample2)
#   
#   # Calculate difference in proportions
#   diff_prop <- abs(prop1 - prop2)
#   
#   # Sum of differences
#   sum_diff <- sum(diff_prop)
#   
#   return(sum_diff)
# }


# import SRKW prey data
rkw_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) %>%
  rename(
    year_n = year
  ) %>%
  filter(
    era == "current"
  )

# rkw_dat_pooled <- readRDS(
#   here::here("data", "rkw_diet", "cleaned_ppn_dat_pooled.rds")
# ) %>% 
#   # downscale utm coordinates to match model
#   mutate(
#     utm_y = utm_y / 1000,
#     utm_x = utm_x / 1000
#   ) %>% 
#   filter(
#     era == "current"
#   ) 


# import original fishery data used to fit model
agg_dat <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
)


# import fitted models
model_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_mvtw.rds"
  )
)
slot_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_slot.rds"
  )
)
large_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_large.rds"
  )
)
# exclude slot model for now because no 2023 data
fit_list <- list(model_fit, slot_fit, large_fit)


# infill rkw data to ensure each stock group present for each sampling event
dat_list <- split(rkw_dat, rkw_dat$sample_id)

new_dat <- purrr::map(
  dat_list,
  function (x) {
    expand.grid(
      sample_id = unique(x$sample_id),
      sample_size = unique(x$n_samples),
      week_n = unique(x$week_n),
      stock_group = levels(model_fit$model$stock_group),
      utm_x = unique(x$utm_x),
      utm_y = unique(x$utm_y),
      year_n = unique(x$year_n) 
    ) %>% 
      # add observed proportion for each stock and sampling event
      left_join(
        ., 
        rkw_dat %>% 
          select(sample_id, stock_group, obs_ppn = agg_prob), 
        by = c("sample_id", "stock_group")
      ) %>% 
      mutate(
        obs_ppn = ifelse(is.na(obs_ppn), 0, obs_ppn),
        sg_year = paste(stock_group, year_n, sep = "_") %>% as.factor(),
        slot_limit = ifelse(
          year_n > 2018, "yes", "no"
        )
      )
  }
) %>% 
  bind_rows()


# generate predictions based on each fitted model and store as tibble
pred_tbl <- tibble(
  dataset = c("standard", "slot", "large")
) %>% 
  mutate(
    pred_dat = purrr::map(
      fit_list, function (x) {
        preds <- predict(
          x, se.fit = TRUE, category_name = "stock_group", origdata = x$model,
          newdata = new_dat
        )
        new_dat %>% 
          mutate(
            fit = preds$fit %>% as.numeric(),
            se = preds$se.fit %>% as.numeric()
          )
      }
    )
  )


# function to generate simulations for a dataframe of predictions
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


# generate simulations for each model fit
future::plan(future::multisession, workers = 6)
pred_tbl$sim_dat <- furrr::future_map(
  pred_tbl$pred_dat, ~ sim_foo(.x, nsim = 500),
  .options = furrr::furrr_options(seed = TRUE)
)


# calculate p-values for each and print
purrr::map(
  pred_tbl$sim_dat, 
  ~ .x %>% 
    group_by(stock_group) %>% 
    summarize(p_value = sum(test_stat) / length(unique(sim_i)))
  )


# calculate simulated proportion in each simulation to compare to observed
sim_ppn_dat <- pred_tbl %>% 
  mutate(
    dataset = factor(dataset, levels = c("standard", "slot", 
                                         "large")),
    sim_ppn = purrr::map(
      sim_dat,
      ~ .x %>% 
        group_by(sim_i) %>% 
        mutate(
          pred_sim_ppn = sample1_count / sum(sample1_count)
        )
    )
  ) %>% 
  select(dataset, sim_ppn) %>% 
  unnest(
    cols = sim_ppn
  )


# calculate observed proportions
obs_ppn_dat <- rkw_dat %>% 
  group_by(stock_group) %>% 
  summarize(
    sum_obs = sum(agg_prob)
  ) %>% 
  ungroup() %>% 
  mutate(
    obs_ppn = sum_obs / sum(sum_obs)
  )

dataset_pal <- c("#e0f3db", "#a8ddb5", 
                 "#43a2ca")
names(dataset_pal) <- levels(sim_ppn_dat$dataset)

sel_boxplot <- ggplot() +
  geom_boxplot(
    data = sim_ppn_dat,
    aes(x = stock_group, y = pred_sim_ppn, fill = dataset)
    ) +
  geom_point(
    data = obs_ppn_dat, 
    aes(x = stock_group, y = obs_ppn), 
    colour = "red"
  ) +
  labs(y = "Simulated Sample Composition") +
  facet_wrap(~dataset, ncol = 1) +
  scale_fill_manual(values = dataset_pal, name = "Model") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

png(
  here::here("figs", "selectivity", "selectivity_boxplot.png"),
  height = 7.5, width = 7.5, units = "in", res = 250
)
sel_boxplot
dev.off()
