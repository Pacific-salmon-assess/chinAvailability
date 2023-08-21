## Model Fitting
# Fit data to estimate seasonal changes in stock composition and determine 
# whether comp similar inside/outside core habitat areas
# cap_region defines whether or not samples were collected in coarse region
# corresponding to habitat polygo
# rkw_habitat defines whether or not samples were collected within foraging 
# habitat polygon

library(tidyverse)
library(stockseasonr)

source(here::here("R", "utils.R"))

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  # redefine region based on analysis
  mutate(
    cap_region = case_when(
      lat < 48.8 & lon > -125.25 & lon < -124.25 ~ "swiftsure",
      lat < 48.45 & lon < -123.4 & lon > -124.25 ~ "sooke",
      TRUE ~ "outside"
    ),
    whale_samples_time = ifelse(
      (year < 2011 | year > 2017) & month_n %in% c("6", "7", "8") & 
        rkw_habitat == "yes",
      "yes",
      "no"
    ),
    year = as.factor(year),
    yday = lubridate::yday(date),
    day_samp = paste(yday, year)
  )

comp_in <- rec_raw %>% 
  filter(
    !legal == "sublegal",
    !cap_region == "outside"
  ) %>% 
  mutate(
    sample_id = paste(cap_region, week_n, rkw_habitat,
                      year, sep = "_")
  ) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)) %>% as.numeric,
    # add abbreviated key for model fitting
    stock_group2 = tolower(stock_group) %>% 
      str_replace_all(., " ", "_") %>% 
      paste("stock", ., sep = "-")
  ) %>% 
  ungroup() %>% 
  group_by(sample_id, 
           rkw_habitat,
           cap_region, week_n, month_n, year, nn, 
           stock_group2) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  filter(
    month_n >= 5 & month_n <= 9
  )

# subset predicted composition dataset
pred_dat_comp <- expand.grid(
  cap_region = unique(comp_in$cap_region),
  rkw_habitat = unique(comp_in$rkw_habitat),
  month_n = seq(min(comp_in$month_n),
                max(comp_in$month_n),
                by = 0.1)
  )


## stockseasonr fits -----------------------------------------------------------

fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + cap_region + rkw_habitat +
    s(month_n, bs = "tp", k = 4) + (1 | year),
  comp_dat = comp_in,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  # nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)

fit$ssdr


# make predictions
pred_fit <- clean_pred_foo(fit = fit, preds = pred_dat_comp)
preds <- pred_fit$preds %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    )
  )
obs <- pred_fit$obs_dat %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    )
  )
  

# seasonal patterns
p <- ggplot(data = pred_fit$preds, aes(x = month_n, colour = rkw_habitat)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) 

p_ribbon <- p +
  geom_ribbon(data = pred_fit$preds,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = rkw_habitat),
              alpha = 0.2)

p_obs <- p +
  geom_jitter(data = obs,
              aes(x = month_n, y = obs_ppn, size = samp_nn), alpha = 0.1) +
  scale_size_continuous() +
  theme(legend.position = "top")


# pull parameter estimates for rkw habitat effects by stock
habitat_effs <- pred_fit$pars %>% 
  filter(grepl("rkw_", parameter)) %>% 
  mutate(
    lo = estimate + qnorm(0.025) * std_error,
    up = estimate + qnorm(0.975) * std_error,
    stock = str_remove(parameter, "rkw_habitatyes_") %>% 
      factor(
        .,
        levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer_4.1",
                   "fraser_fall", "ecvi_somn", "psd", "other")
      )
  )

dot_plot <- ggplot(habitat_effs) +
  geom_pointrange(aes(x = stock, y = estimate, ymin = lo, ymax = up)) + 
  ggsidekick::theme_sleek() +
  geom_hline(aes(yintercept = 0), lty = 2) +
  labs(y = "Habitat Effect") +
  ggsidekick::theme_sleek()


pdf(here::here("figs", "rkw_habitat", "stratified_habitat_effects.pdf"))
p_ribbon
p_obs
dot_plot
dev.off()

# WCVI and Summer 4_1 fish more abundant inside zones than out