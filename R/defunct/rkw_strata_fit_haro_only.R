## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Extension of rkw_strata_fit including smaller number of stocks for Haro


library(tidyverse)
library(stockseasonr)

source(here::here("R", "utils.R"))

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) 

comp_in <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    !legal == "sublegal",
    #exclude samples collected outside areas in relatively close proximity to 
    # SRKW foraging areas
    !rkw_habitat == "outside",
    !is.na(strata)
  ) %>% 
  mutate(
    sample_id = paste(strata, week_n, year, sep = "_"),
    strata = as.factor(strata),
    #consolidate stocks
    stock_group = ifelse(
      stock_group %in% c("other", "NBC_SEAK"), "other", stock_group
      )
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
           strata,
           week_n, month_n, year, nn, 
           stock_group2) %>% 
  summarize(prob = sum(prob), .groups = "drop") 


# make distinct geographic datasets for model fitting
haro_dat <- comp_in %>% 
  mutate(
    strata2 = ifelse(grepl("n_haro", strata), "n_haro", as.character(strata)),
    strata = as.factor(strata2)
  ) %>% 
  filter(
    # grepl("vic", strata) |
    grepl("n_haro", strata) |# grepl("saanich", strata) |
      grepl("s_haro", strata),
    month_n > 3 & month_n < 10
  ) 


# make tibble
pred_dat <- expand.grid(
    strata = unique(haro_dat$strata),
    month_n = seq(min(haro_dat$month_n), max(haro_dat$month_n), by = 0.1)
  ) 


## stockseasonr fits -----------------------------------------------------------

haro_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + strata +
    s(month_n, bs = "tp", k = 3) +
    # s(month_n, by = strata, bs = "tp", k = 3, m = 1) +
    (1 | year),
  comp_dat = haro_dat,
  pred_dat = pred_dat,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)


# make predictions
# pred_swift <- clean_pred_foo(fit = swift_fit, preds = dat_tbl$pred_dat[[1]])
# pred_sooke <- clean_pred_foo(fit = sooke_fit, preds = dat_tbl$pred_dat[[2]])
pred_haro <- clean_pred_foo(fit = haro_fit, preds = pred_dat)


preds <- pred_haro$preds %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    ),
    stock = fct_recode(
      stock, "fr_yearling" = "fraser_yearling", "puget" = "psd",
      "fr_sum_4.1" = "fraser_summer_4.1", "fr_fall" = "fraser_fall"
    )
  )
obs <- pred_haro$obs_dat %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    ),
    stock = fct_recode(
      stock, "fr_yearling" = "fraser_yearling", "puget" = "psd",
      "fr_sum_4.1" = "fraser_summer_4.1", "fr_fall" = "fraser_fall"
    )
  )



# seasonal patterns inside/outside rkw_habitat
p_hab <- ggplot(data = preds, 
                aes(x = month_n, colour = strata)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) 

p_ribbon_hab <- p_hab +
  geom_ribbon(data = preds,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = strata),
              alpha = 0.2) 

p_obs_hab <- ggplot() +
  geom_jitter(data = obs,
              aes(x = month_n, y = obs_ppn, colour = strata, 
                  size = samp_nn),
              alpha = 0.1) +
  geom_line(data = preds,
            aes(x = month_n, colour = strata, y = pred_prob_est),
            linewidth = 1.25) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  scale_size_continuous() +
  theme(legend.position = "top")


# seasonal patterns across sizes
p_size <- ggplot(data = preds %>% filter(rkw_habitat == "yes"), 
                 aes(x = month_n, colour = size)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) +
  scale_color_brewer(type = "qual", palette = 2)

p_ribbon_size <- p_size +
  geom_ribbon(data = preds %>% filter(rkw_habitat == "yes"),
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = size),
              alpha = 0.2) +
  scale_fill_brewer(type = "qual", palette = 2)

p_obs_size <- ggplot() +
  geom_jitter(data = obs,
              aes(x = month_n, y = obs_ppn, colour = size, size = samp_nn),
              alpha = 0.1) +
  geom_line(data = preds %>% filter(rkw_habitat == "yes"),
            aes(x = month_n, colour = size, y = pred_prob_est), 
            linewidth = 1.25) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  scale_color_brewer(type = "qual", palette = 2) +
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

plot_list <- list(p_ribbon_hab, p_obs_hab, p_ribbon_size, p_obs_size, dot_plot)
saveRDS(
  plot_list, 
  here::here("figs", "rkw_habitat", "stratified_habitat_size_effects.rds")
)

pdf(here::here("figs", "rkw_habitat", "stratified_habitat_size_effects.pdf"),
    height =  4.5, width = 8)
p_ribbon_hab
p_obs_hab
p_ribbon_size
p_obs_size
dot_plot
dev.off()


# WCVI and Summer 4_1 fish more abundant inside zones than out