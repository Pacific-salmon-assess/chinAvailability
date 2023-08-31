## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Extension of rkw_habitat_fit including more granular spatial strata
# 1) Fit 3 models that include different spatial strata 
# 2) Within each spatial strata also fit alternates that stratify by size
# Collapse strata in subsequent analyses if no significant differences
# Also determine how fitting overlapping spatial models influences predictions
# (e.g linking S Haro to Sooke vs. to N Haro)


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
    !strata == "saanich",
    !is.na(strata)
  ) %>% 
  mutate(
    sample_id = paste(strata, week_n, year, sep = "_"),
    strata = as.factor(strata),
    #consolidate stocks
    stock_group = ifelse(
      stock_group %in% c("other", "NBC_SEAK"), "other", stock_group
    ),
    strata_region = ifelse(
      grepl("vic", strata) | grepl("sooke", strata) |
        grepl("n_haro", strata) | grepl("s_haro", strata),
      "east",
      "west"
    ) %>% 
      factor()
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
           strata, strata_region,
           week_n, month_n, year, nn, 
           stock_group2) %>% 
  summarize(prob = sum(prob), .groups = "drop") 


# make distinct geographic datasets for model fitting
swift_dat <- comp_in %>% 
  filter(
    grepl("renfrew", strata) | grepl("swiftsure", strata) | 
      grepl("bark", strata) | grepl("nitinat", strata)
  )
sooke_dat <- comp_in %>% 
  filter(
    grepl("sooke", strata) | grepl("vic", strata) | grepl("s_haro", strata),
    month_n > 3 & month_n < 10
  )
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
  ) %>% 
  select(-strata2)
east_dat <- comp_in %>% 
  mutate(
    strata2 = ifelse(grepl("n_haro", strata), "n_haro", as.character(strata)),
    strata = as.factor(strata2)
  ) %>% 
  filter(
    grepl("vic", strata) | grepl("sooke", strata) |
    grepl("n_haro", strata) | grepl("s_haro", strata),
    month_n > 3 & month_n < 10
  ) %>% 
  select(-strata2)
all_dat <- comp_in %>% 
  filter(
    month_n > 3 & month_n < 10
  )

# make tibble
dat_tbl <- tibble(
  group = c("swiftsure", "sooke", "haro", "east", "full"),
  data = list(swift_dat, sooke_dat, haro_dat, east_dat, all_dat)
) %>% 
  mutate(
    pred_dat = purrr::map(
      data,
      ~ expand.grid(
        strata = unique(.x$strata),
        month_n = seq(min(.x$month_n), max(.x$month_n), by = 0.1)
      ) %>% 
        left_join(., all_dat %>% select(strata, strata_region) %>% distinct(),
                  by = "strata")
    )#,
    # formula = ifelse(
    #   group == "swif"
    # )
  )


## stockseasonr fits -----------------------------------------------------------

swift_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + strata + s(month_n, bs = "tp", k = 3) + 
    (1 | year),
  comp_dat = dat_tbl$data[[1]],
  pred_dat = dat_tbl$pred_dat[[1]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)
# sooke_fit <- fit_stockseasonr(
#   comp_formula = stock_group2 ~ 1 + strata + s(month_n, bs = "cc", k = 4, m = 2) + 
#     # s(month_n, by = strata, bs = "cc", k = 4, m = 1) +
#     (1 | year),
#   comp_dat = dat_tbl$data[[2]],
#   pred_dat = dat_tbl$pred_dat[[2]],
#   model = "dirichlet",
#   random_walk = FALSE,
#   fit = TRUE,
#   newton_loops = 1,
#   silent = FALSE
# )
sooke_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 +
    strata + s(month_n, bs = "tp", k = 3) + 
    (1 | year),
  comp_dat = dat_tbl$data[[2]],
  pred_dat = dat_tbl$pred_dat[[2]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)
# haro_fit <- fit_stockseasonr(
#   comp_formula = stock_group2 ~ 1 + strata +
#     s(month_n, bs = "tp", k = 3) +
#     (1 | year),
#   comp_dat = dat_tbl$data[[3]],
#   pred_dat = dat_tbl$pred_dat[[3]],
#   model = "dirichlet",
#   random_walk = FALSE,
#   fit = TRUE,
#   newton_loops = 1,
#   silent = FALSE
# )
east_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + strata +
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = dat_tbl$data[[4]],
  pred_dat = dat_tbl$pred_dat[[4]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)
all_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + strata +
    s(month_n, bs = "tp", k = 4, m = 2, by = strata_region) +
    # s(month_n, bs = "tp", k = 4, m = 1, by = strata_region) +
    (1 | year),
  comp_dat = dat_tbl$data[[5]],
  pred_dat = dat_tbl$pred_dat[[5]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)


# make predictions
pred_swift <- clean_pred_foo(fit = swift_fit, preds = dat_tbl$pred_dat[[1]])
pred_sooke <- clean_pred_foo(fit = sooke_fit, preds = dat_tbl$pred_dat[[2]])
pred_all <- clean_pred_foo(fit = all_fit, preds = dat_tbl$pred_dat[[5]])
pred_east <- clean_pred_foo(fit = east_fit, preds = dat_tbl$pred_dat[[4]])

preds <- rbind(
  pred_swift$preds %>% 
    mutate(region = "swift"),
  pred_sooke$preds %>% 
    mutate(region = "sooke"),
  pred_all$preds %>%
    mutate(region = "all"),
  pred_east$preds %>%
    mutate(region = "east")
  ) %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    ),
    stock = fct_recode(
      stock, "fr_yearling" = "fraser_yearling", "puget" = "psd",
      "fr_sum_4.1" = "fraser_summer_4.1", "fr_fall" = "fraser_fall"
    )
  ) %>% 
  filter(strata_region == "west")
obs <- rbind(
  pred_swift$obs_dat %>% 
    mutate(region = "swift"),
  pred_sooke$obs_dat %>% 
    mutate(region = "sooke"),
  pred_all$obs_dat %>%
    mutate(region = "all"),
  pred_east$obs_dat %>%
    mutate(region = "east")
) %>% 
  mutate(
    stock = factor(
      stock,
      levels = c("wcvi", "fraser_yearling", "fraser_summer_4.1",
                 "fraser_fall", "ecvi_somn", "psd", "other")
    ),
    stock = fct_recode(
      stock, "fr_yearling" = "fraser_yearling", "puget" = "psd",
      "fr_sum_4.1" = "fraser_summer_4.1", "fr_fall" = "fraser_fall"
    )
  ) %>% 
  filter(strata_region == "west")
  

# seasonal patterns inside/outside rkw_habitat
p_hab <- ggplot(data = preds, 
                 aes(x = month_n, colour = strata)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(region~stock) +
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
  facet_grid(region~stock) +
  ggsidekick::theme_sleek() +
  scale_size_continuous() +
  theme(legend.position = "top")



pdf(here::here("figs", "rkw_habitat", "strata_effects_pooling_east.pdf"),
    height =  4.5, width = 8)
p_ribbon_hab
p_obs_hab
dev.off()

