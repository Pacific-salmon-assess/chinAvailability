## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Extension of rkw_strata_fit, but constrained to only large (>75 cm FL)


library(tidyverse)
library(stockseasonr)

source(here::here("R", "utils.R"))

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region)

rec_trim <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    !legal == "sublegal",
    #exclude samples collected outside areas in relatively close proximity to 
    # SRKW foraging areas
    !rkw_habitat == "outside",
    !strata == "saanich",
    !is.na(strata),
    fl >= (700),
    month_n >= 5 & month_n <= 10
  ) %>% 
  mutate(
    sample_id = paste(strata, week_n, year, sep = "_"),
    strata = ifelse(grepl("n_haro", strata), "n_haro", strata) %>% 
      factor(),
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
  )

comp_in <- rec_trim  %>% 
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



## MODEL FIT  ------------------------------------------------------------------


# make distinct geographic datasets for model fitting
# swift_dat <- comp_in %>% 
#   filter(
#     strata_region == "west"
#   ) %>% 
#   droplevels()
# east_dat <- comp_in %>% 
#   filter(strata_region == "east") %>% 
#   droplevels()
# 
# 
# # make tibble
# dat_tbl <- tibble(
#   group = c("swiftsure", "east"),
#   data = list(swift_dat, east_dat)
# ) %>% 
#   mutate(
#     pred_dat = purrr::map(
#       data,
#       ~ expand.grid(
#         strata = unique(.x$strata),
#         month_n = seq(min(.x$month_n), max(.x$month_n), by = 0.1)
#       ) 
#     )
#   )

pred_dat <- expand.grid(
    strata = unique(comp_in$strata),
    month_n = seq(min(comp_in$month_n), max(comp_in$month_n), by = 0.1)
  ) %>% 
    left_join(., comp_in %>% select(strata, strata_region) %>% distinct(),
              by = "strata")



## stockseasonr fits -----------------------------------------------------------

all_fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + strata +
    s(month_n, bs = "tp", k = 3, m = 2, by = strata_region) +
    (1 | year),
  comp_dat = comp_in,
  pred_dat = pred_dat,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)


# make predictions
preds_all <- clean_pred_foo(fit = all_fit, preds = pred_dat)

preds <- preds_all$preds %>% 
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
  ) 
obs <- preds_all$obs_dat %>% 
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
  ) 


p <- ggplot(data = preds, 
                 aes(x = month_n, colour = strata)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(strata_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) 

p_ribbon <- p +
  geom_ribbon(data = preds,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = strata),
              alpha = 0.2) 

p_obs <- ggplot() +
  geom_jitter(data = obs,
              aes(x = month_n, y = obs_ppn, colour = strata, 
                  size = samp_nn),
              alpha = 0.1) +
  geom_line(data = preds,
            aes(x = month_n, colour = strata, y = pred_prob_est),
            linewidth = 1.25) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(strata_region~stock) +
  ggsidekick::theme_sleek() +
  scale_size_continuous() +
  theme(legend.position = "top")


west_stacked <- ggplot(data = preds %>% filter(strata_region == "swift"), 
                       aes(x = month_n)) +
  geom_area(aes(y = pred_prob_est, colour = stock, fill = stock), 
            stat = "identity") +
  scale_fill_brewer(name = "Stock", palette = "Spectral") +
  scale_colour_brewer(name = "Stock", palette = "Spectral") +
  labs(y = "Predicted Composition", x = "Month") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size=9),
        plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points")
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA), xlim = c(1, 12)) +
  facet_wrap(~strata)

east_stacked <- ggplot(data = preds %>% filter(strata_region == "east"), 
       aes(x = month_n)) +
  geom_area(aes(y = pred_prob_est, colour = stock, fill = stock), 
            stat = "identity") +
  scale_fill_brewer(name = "Stock", palette = "Spectral") +
  scale_colour_brewer(name = "Stock", palette = "Spectral") +
  labs(y = "Predicted Composition", x = "Month") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size=9),
        plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points")
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA), xlim = c(1, 12)) +
  facet_wrap(~strata)

pdf(here::here("figs", "rkw_habitat", "large_preds.pdf"),
    height = 7, width = 8.5)
p
p_ribbon
p_obs
west_stacked
east_stacked
dev.off()
