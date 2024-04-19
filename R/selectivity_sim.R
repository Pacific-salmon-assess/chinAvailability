## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)
library(mvtweedie)

# modified tweedie prediction that allows for exclude = list()
source(here::here("R", "functions", "pred_mvtweedie2.R"))


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


# import SRKW prey data
rkw_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) %>%
  rename(
    week_n = week,
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


# import fitted model
model_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
  )
)

#### DEFUNCT IF GENERATING PREDICTIONS ALL AT ONCE #####
# split diet data by sample ID and simulate
# rkw_dat_tbl <- rkw_dat_pooled %>%
#   group_by(sample_id_pooled) %>%
#   group_nest()

# make new data frame consistent with SRKW diet to generate predictions from
# new_dat_list <- purrr::map(
#   rkw_dat_tbl$data,
#   function (x) {
#     newdata <- expand.grid(
#       sample_size = unique(x$n_samples),
#       week_n = unique(x$week_n),
#       stock_group = levels(model_fit$model$stock_group),
#       utm_x = unique(x$utm_x),
#       utm_y = unique(x$utm_y),
#       year_n = "2020" # dummy year so value doesn't matter
#     ) %>%
#       left_join(
#         .,
#         x %>%
#           select(stock_group, obs_prob = agg_prob) %>%
#           distinct(),
#         by = c("stock_group")
#       ) %>%
#       mutate(
#         obs_prob = ifelse(is.na(obs_prob), 0, obs_prob)
#       )
# 
#     # exclude year smoother and generate predictions based on fishery model
#     excl <- grepl("year_n", gratia::smooths(model_fit))
#     yr_coefs <- gratia::smooths(model_fit)[excl]
#     preds <- pred_dummy(
#       model_fit,
#       se.fit = TRUE,
#       category_name = "stock_group",
#       origdata = agg_dat,
#       newdata = newdata,
#       exclude = yr_coefs
#     )
# 
#     cbind(
#       newdata, pred_prob = preds$fit, se_pred = preds$se.fit
#     ) %>%
#       distinct()
#   }
# )

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
        obs_ppn = ifelse(is.na(obs_ppn), 0, obs_ppn)
      )
  }
) %>% 
  bind_rows()


# calculate mean ppn for each sampling event
# excl <- grepl("year_n", gratia::smooths(model_fit))
# yr_coefs <- gratia::smooths(model_fit)[excl]
preds <- pred_dummy(
  model_fit,
  se.fit = FALSE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = new_dat
)
new_dat$pred_ppn <- preds

# split by new sampling event
new_dat_list <- split(new_dat, new_dat$sample_id)

# simulate based on each observation dataset 
nsim <- 1000
sim_dat_out <- vector(mode = "list", length = nsim)

for (i in 1:nsim) {
  # generate null multinomial distribution for each sampling event
  sim_dat <- purrr::map(
    new_dat_list,
    function (x) {
      # Generate multinomial samples under null hypothesis (i.e. same dist)
      x %>% 
        mutate(
          nsim = i,
          sample_size = unique(x$sample_size),
          sample1_null = rmultinom(1, sample_size, x$pred_ppn) %>% as.numeric(),
          sample2_null = rmultinom(1, sample_size, x$pred_ppn) %>% as.numeric(),
          sample_obs = rmultinom(1, sample_size, x$obs_ppn) %>% as.numeric()
        )
    } 
  ) %>% 
    bind_rows()
  
  # for simulation calculate total counts then test whether obs greater than
  # null
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
    ) 
  sim_dat_out[[i]] <- sim_dat_agg
}

sim_dat_out %>% 
  bind_rows() %>% 
  group_by(stock_group) %>% 
  summarize(p_value = sum(test_stat) / nsim)

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