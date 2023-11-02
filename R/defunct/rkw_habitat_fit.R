## Model Fitting
# Fit data to estimate seasonal changes in stock composition and determine 
# whether comp similar inside/outside core habitat areas
# cap_region defines whether or not samples were collected in coarse region
# corresponding to habitat polygo
# rkw_habitat defines whether or not samples were collected within foraging 
# habitat polygon
# fit two models: all legal sized fish and only samples greater than 75 cm

library(tidyverse)
library(stockseasonr)

source(here::here("R", "utils.R"))

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) 

comp_in <- rec_raw %>% 
  filter(
    !legal == "sublegal",
    !cap_region == "outside"
  ) %>% 
  mutate(
    sample_id = paste(cap_region, week_n, rkw_habitat,
                      year, sep = "_"),
    cap_region = as.factor(cap_region)
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

comp_in_large <- rec_raw %>% 
  filter(
    !legal == "sublegal",
    !cap_region == "outside",
    fl >= 750
  ) %>% 
  mutate(
    sample_id = paste(cap_region, week_n, rkw_habitat,
                      year, sep = "_"),
    cap_region = as.factor(cap_region)
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
pred_dat_comp_week <- expand.grid(
  cap_region = unique(comp_in$cap_region),
  rkw_habitat = unique(comp_in$rkw_habitat),
  week_n = seq(min(comp_in$week_n),
                max(comp_in$week_n),
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
  newton_loops = 1,
  silent = FALSE
)

# week as covariate had minor effect on predictions 
# fit2 <- fit_stockseasonr(
#   comp_formula = stock_group2 ~ 1 + cap_region + rkw_habitat +
#     s(week_n, bs = "tp", k = 4) + (1 | year),
#   comp_dat = comp_in,
#   pred_dat = pred_dat_comp_week,
#   model = "dirichlet",
#   random_walk = FALSE,
#   fit = TRUE,
#   newton_loops = 1,
#   silent = FALSE
# )

fit_large <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + cap_region + rkw_habitat +
    s(month_n, bs = "tp", k = 4) + (1 | year),
  comp_dat = comp_in_large,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)
head(fit$ssdr)

# global smooth + regional smooth doesn't converge and regional smooth only 
# results in unrealistic predictions outside sampling domains
# fit2 <- fit_stockseasonr(
#   comp_formula = stock_group2 ~ 1 + cap_region + rkw_habitat +
#     s(month_n, bs = "tp", k = 4, m = 2) +
#     s(month_n, by = cap_region, bs = "tp", k = 4, m = 1) + (1 | year),
#   comp_dat = comp_in,
#   pred_dat = pred_dat_comp,
#   model = "dirichlet",
#   random_walk = FALSE,
#   fit = TRUE,
#   # nlminb_loops = 2,
#   newton_loops = 1,
#   silent = FALSE
# )
# 
# head(fit2$ssdr)


# make predictions
pred_fit <- clean_pred_foo(fit = fit2, preds = pred_dat_comp)
pred_fit_large <- clean_pred_foo(fit = fit_large, preds = pred_dat_comp)

preds <- rbind(
  pred_fit$preds %>% 
    mutate(size = "all"),
  pred_fit_large$preds %>% 
    mutate(size = "large")
  ) %>% 
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
obs <- rbind(
  pred_fit$obs_dat %>% 
    mutate(size = "all"),
  pred_fit_large$obs_dat %>% 
    mutate(size = "large")
) %>% 
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
p_hab <- ggplot(data = preds %>% filter(size == "all"), 
                 aes(x = month_n, colour = rkw_habitat)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) 

p_ribbon_hab <- p_hab +
  geom_ribbon(data = preds %>% filter(size == "all"),
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = rkw_habitat),
              alpha = 0.2) 

p_obs_hab <- ggplot() +
  geom_jitter(data = obs,
              aes(x = month_n, y = obs_ppn, colour = rkw_habitat, 
                  size = samp_nn),
              alpha = 0.1) +
  geom_line(data = preds %>% filter(size == "all"),
            aes(x = month_n, colour = rkw_habitat, y = pred_prob_est), 
            linewidth = 1.25) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
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