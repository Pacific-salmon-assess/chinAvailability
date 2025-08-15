## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024

# Set French language option
FRENCH <- FALSE

# Create appropriate figure directories
if (FRENCH) {
  dir.create("figs-french", showWarnings = FALSE)
  dir.create("figs-french/selectivity", showWarnings = FALSE)
  fig_dir <- "figs-french"
} else {
  dir.create("figs/selectivity", showWarnings = FALSE)
  fig_dir <- "figs"
}

# Translation helper function
tr <- function(english, french) {
  if (FRENCH) french else english
}

# Helper function for figure paths
fig_path <- function(filename) {
  file.path(fig_dir, filename)
}

# Translation function for dataset categories
translate_dataset <- function(dataset_values) {
  if (FRENCH) {
    case_when(
      dataset_values == "standard" ~ "standard",
      dataset_values == "large" ~ "grand", 
      dataset_values == "filtered" ~ "filtré",
      TRUE ~ dataset_values
    )
  } else {
    dataset_values
  }
}

library(tidyverse)
library(mvtweedie)

# modified tweedie prediction that allows for exclude = list()
source(here::here("R", "functions", "pred_mvtweedie2.R"))
source(here::here("R", "functions", "sim_multinomial.R"))


dataset_pal <- c("#e0f3db", "#43a2ca")
# Set palette names based on language setting
if (FRENCH) {
  names(dataset_pal) <- translate_dataset(c("standard", "large"))
} else {
  names(dataset_pal) <- c("standard", "large")
}

# import SRKW prey data list
# NOTE: uses mean week and location for samples collected within a given strata;
# exploratory analysis using single observations provided qualitatively similar
# results with much greater uncertainty
rkw_dat1 <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) 

# clean both datasets
rkw_dat <- purrr::map(
  rkw_dat1, 
  ~ .x %>% 
    rename(
      year_n = year
    ) %>%
    filter(
      era == "current"
    ) 
)


# make keys representing proportion of samples in rec fishery by stock and
# size class
stock_sample_key <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
) %>% 
  filter(
    strata %in% c(rkw_dat$stock$strata),
    week_n >= min(rkw_dat$stock$week_n) & week_n <= max(rkw_dat$stock$week_n)
  ) %>% 
  group_by(
    stock_group
  ) %>% 
  summarize(
    nn = sum(agg_prob)
  ) %>% 
  ungroup() %>% 
  mutate(
    total_samp = sum(nn),
    ppn = nn /total_samp
  )

size_sample_key <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_size_xy.rds")
) %>% 
  filter(
    strata %in% c(rkw_dat$size$strata),
    week_n >= min(rkw_dat$size$week_n) & week_n <= max(rkw_dat$size$week_n)
  ) %>% 
  group_by(
    size_bin
  ) %>% 
  summarize(
    nn = sum(agg_prob)
  ) %>% 
  ungroup() %>% 
  mutate(
    total_samp = sum(nn),
    ppn = nn /total_samp
  )


## STOCK SELECTIVITY -----------------------------------------------------------

# import fitted models
model_fit <- readRDS(
  here::here(
    "data", "model_fits", "fit_spatial_fishery_ri_mvtw.rds"
  )
)
large_fit <- readRDS(
  here::here(
    "data", "model_fits", "fit_large.rds"
  )
)
fit_list <- list(model_fit, large_fit)


# infill rkw data to ensure each stock group present for each sampling event
dat_list <- split(rkw_dat$stock, rkw_dat$stock$sample_id)

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
        rkw_dat$stock %>% 
          select(sample_id, stock_group, obs_ppn = agg_prob), 
        by = c("sample_id", "stock_group")
      ) %>% 
      mutate(
        obs_ppn = ifelse(is.na(obs_ppn), 0, obs_ppn),
        sg_year = paste(stock_group, year_n, sep = "_") %>% as.factor()
      )
  }
) %>% 
  bind_rows()


# generate predictions based on each fitted model and store as tibble
pred_tbl <- tibble(
  dataset = c("standard", "large")
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


# generate simulations for each model fit
future::plan(future::multisession, workers = 6)
pred_tbl$sim_dat <- furrr::future_map(
  pred_tbl$pred_dat, ~ sim_foo(.x, nsim = 500),
  .options = furrr::furrr_options(seed = TRUE)
)


# calculate p-values for each and print
# p_val <- purrr::map2(
#   pred_tbl$sim_dat, pred_tbl$dataset,
#   ~ .x %>% 
#     group_by(stock_group) %>% 
#     summarize(p_value = sum(test_stat) / length(unique(sim_i))) %>% 
#     mutate(dataset = .y)
#   ) %>% 
#   bind_rows()
# p_val_sig <- p_val %>% 
#   filter(p_value < 0.05) %>% 
#   mutate(
#     dataset = factor(dataset, levels = c("standard", "large"))
#   )


# calculate simulated proportion in each simulation to compare to observed
sim_ppn_dat <- pred_tbl %>% 
  mutate(
    dataset = factor(dataset, levels = c("standard", "large")),
    sim_ppn = purrr::map(
      sim_dat,
      ~ .x %>% 
        group_by(sim_i) %>% 
        mutate(
          pred_sim_ppn = sample1_count / sum(sample1_count),
          diff_ppn = obs_diff_raw / sum(sample1_count)
        )
    )
  ) %>% 
  select(dataset, sim_ppn) %>% 
  unnest(
    cols = sim_ppn
  )


dd <- sim_ppn_dat %>%
  group_by(stock_group, dataset) %>% 
  summarize(
    med_dif = median(diff_ppn),
    up_dif = rethinking::HPDI(diff_ppn, 0.95)[2],
    lo_dif = rethinking::HPDI(diff_ppn, 0.95)[1]
  ) %>% 
  left_join(
    ., stock_sample_key %>% select(stock_group, ppn), by = "stock_group"
  )

# export for use in other sensitivity analysis
saveRDS(dd %>% filter(dataset == "standard"),
        here::here("data", "rec", "original_stock_selectivity.rds"))


sel_bean <- ggplot() +
  geom_pointrange(
    data = dd %>% filter(dataset == "standard"),
    aes(x = med_dif, xmin = lo_dif, xmax = up_dif, y = stock_group,
        fill = ppn),
    shape = 21
  ) +
  scale_fill_continuous(
    name = tr("Proportion of\nFishery Samples\nin Western Strata", "Proportion des\néchantillons de pêche\ndans les strates occidentales"),
    trans = "sqrt",
    breaks = c(0.05, 0.15, 0.25)
    ) +
  labs(x = tr("Difference Between Observed and Predicted Composition", "Différence entre la composition observée et prédite"),
       y = tr("Stock", "Stock"),
       fill = tr("Proportion of Fishery Samples in Western Strata", "Proportion des échantillons de pêche dans les strates occidentales")) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        legend.key.size = unit(0.75, "cm"),
        plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5)) +
  geom_vline(xintercept = 0, lty = 2)

sel_bean2 <- ggplot() +
  geom_pointrange(
    data = dd %>% mutate(dataset = translate_dataset(dataset)),
    aes(x = med_dif, xmin = lo_dif, xmax = up_dif, y = stock_group, 
        fill = ppn),
    shape = 21
  ) +
  scale_fill_continuous(
    name = tr("Proportion of\nFishery Samples\nin Western Strata", "Proportion des\néchantillons de pêche\ndans les strates occidentales"),
    trans = "sqrt",
    breaks = c(0.05, 0.15, 0.25)
  ) +
  labs(x = tr("Difference Between Observed and Predicted Composition", "Différence entre la composition observée et prédite"),
       y = tr("Stock", "Stock"),
       fill = tr("Proportion of Fishery Samples in Western Strata", "Proportion des échantillons de pêche dans les strates occidentales")) +
  facet_wrap(~dataset, ncol = 1) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        legend.key.size = unit(0.75, "cm"),
        # legend.text = element_text(size = 9),
        plot.margin = margin(t = 5.5, r = 10, b = 5.5, l = 5.5)) +
  geom_vline(xintercept = 0, lty = 2)

png(
  fig_path("selectivity/selectivity_bean_stock.png"),
  height = 6.5, width = 5, units = "in", res = 250
)
sel_bean
dev.off()

png(
  fig_path("selectivity/selectivity_bean_stock_comp.png"),
  height = 7.5, width = 5, units = "in", res = 250
)
sel_bean2
dev.off()


# calculate observed proportions
# obs_ppn_dat <- rkw_dat %>% 
#   group_by(stock_group) %>% 
#   summarize(
#     sum_obs = sum(agg_prob)
#   ) %>% 
#   ungroup() %>% 
#   mutate(
#     obs_ppn = sum_obs / sum(sum_obs),
#     se_obs_ppn = sqrt(obs_ppn * (1 - obs_ppn) / sum(sum_obs)),
#     lo = pmax(0, obs_ppn - (1.96 * se_obs_ppn)),
#     up = pmin(1, obs_ppn + (1.96 * se_obs_ppn))
#   )
# 
# sel_boxplot <- ggplot() +
#   geom_boxplot(
#     data = sim_ppn_dat %>% filter(dataset == "standard"),
#     aes(x = stock_group, y = pred_sim_ppn)
#     ) +
#   geom_pointrange(
#     data = obs_ppn_dat, 
#     aes(x = stock_group, y = obs_ppn, ymin = lo, ymax = up), 
#     colour = "red",
#     position = position_nudge(x = 0.1)
#     ) +
#   geom_text(
#     data = p_val_sig %>% filter(dataset == "standard"),
#     aes(x = stock_group, y = max(sim_ppn_dat$pred_sim_ppn + 0.1)), 
#     label = "*", size = 7.5
#   ) +
#   lims(y = c(0, max(sim_ppn_dat$pred_sim_ppn + 0.1))) +
#   labs(y = "Sample Composition") +
#   ggsidekick::theme_sleek() +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# sel_boxplot2 <- ggplot() +
#   geom_boxplot(
#     data = sim_ppn_dat,
#     aes(x = stock_group, y = pred_sim_ppn, fill = dataset)
#   ) +
#   geom_pointrange(
#     data = obs_ppn_dat, 
#     aes(x = stock_group, y = obs_ppn, ymin = lo, ymax = up), 
#     colour = "red",
#     position = position_nudge(x = 0.1)
#   ) +
#   geom_text(
#     data = p_val_sig,
#     aes(x = stock_group, y = max(sim_ppn_dat$pred_sim_ppn + 0.1)), 
#     label = "*", size = 7.5
#   ) +
#   lims(y = c(0, max(sim_ppn_dat$pred_sim_ppn + 0.13))) +
#   labs(y = "Sample Composition") +
#   facet_wrap(~dataset, ncol = 1) +
#   scale_fill_manual(values = dataset_pal, name = "Model") +
#   ggsidekick::theme_sleek() +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# png(
#   here::here("figs", "selectivity", "selectivity_boxplot_stock.png"),
#   height = 4.5, width = 7.5, units = "in", res = 250
# )
# sel_boxplot
# dev.off()
# 
# png(
#   here::here("figs", "selectivity", "selectivity_boxplot_stock_comp.png"),
#   height = 7.5, width = 7.5, units = "in", res = 250
# )
# sel_boxplot2
# dev.off()



## SIZE SELECTIVITY ------------------------------------------------------------

# import fitted models
fit_size <- readRDS(
  here::here(
    "data", "model_fits", "fit_size_mvtw.rds"
  )
)
fit_list <- list(fit_size)


size_bins <- levels(fit_size$model$size_bin)

# infill rkw data to ensure each stock group present for each sampling event
rkw_size_dat <- rkw_dat$size 
dat_list_size <- split(rkw_size_dat, rkw_size_dat$sample_id)

new_dat_size <- purrr::map(
  dat_list_size,
  function (x) {
    expand.grid(
      sample_id = unique(x$sample_id),
      sample_size = unique(x$n_samples),
      week_n = unique(x$week_n),
      size_bin = size_bins,
      utm_x = unique(x$utm_x),
      utm_y = unique(x$utm_y),
      year_n = unique(x$year_n) 
    ) %>% 
      # add observed proportion for each stock and sampling event
      left_join(
        ., 
        rkw_size_dat %>% 
          select(sample_id, size_bin, obs_ppn = agg_prob), 
        by = c("sample_id", "size_bin")
      ) %>% 
      mutate(
        obs_ppn = ifelse(is.na(obs_ppn), 0, obs_ppn),
        sg_year = paste(size_bin, year_n, sep = "_") %>% as.factor()
      )
  }
) %>% 
  bind_rows()


# generate predictions based on each fitted model and store as tibble
pred_tbl_size <- tibble(
  dataset = c("standard")
) %>% 
  mutate(
    pred_dat = purrr::map(
      fit_list, function (x) {
        preds <- pred_dummy(
          x, se.fit = TRUE, category_name = "size_bin", origdata = x$model,
          newdata = new_dat_size
        )
        new_dat_size %>% 
          mutate(
            fit = preds$fit %>% as.numeric(),
            se = preds$se.fit %>% as.numeric()
          )
      }
    )
  )


# generate simulations for each model fit
future::plan(future::multisession, workers = 6)
pred_tbl_size$sim_dat <- furrr::future_map(
  pred_tbl_size$pred_dat, ~ sim_foo(.x, category_name = "size_bin", nsim = 500),
  .options = furrr::furrr_options(seed = TRUE)
)


# calculate p-values for each and print
# p_val_size <- purrr::map2(
#   pred_tbl_size$sim_dat, pred_tbl_size$dataset,
#   ~ .x %>% 
#     group_by(size_bin) %>% 
#     summarize(p_value = sum(test_stat) / length(unique(sim_i))) %>% 
#     mutate(dataset = .y)
# ) %>% 
#   bind_rows()
# p_val_sig_size <- p_val_size %>% 
#   filter(p_value < 0.05) %>% 
#   mutate(
#     dataset = factor(dataset, levels = c("standard"))
#   )


# calculate simulated proportion in each simulation to compare to observed
sim_ppn_dat_size <- pred_tbl_size %>% 
  mutate(
    sim_ppn = purrr::map(
      sim_dat,
      ~ .x %>% 
        group_by(sim_i) %>% 
        mutate(
          pred_sim_ppn = sample1_count / sum(sample1_count),
          diff_ppn = obs_diff_raw / sum(sample1_count)
        )
    )
  ) %>% 
  select(sim_ppn) %>% 
  unnest(
    cols = sim_ppn
  )


diff_quantile <- sim_ppn_dat_size %>%
  group_by(size_bin) %>% 
  summarize(
    med_dif = median(diff_ppn),
    up_dif = quantile(diff_ppn, 0.975),
    lo_dif = quantile(diff_ppn, 0.025)
  ) %>% 
  left_join(
    ., size_sample_key %>% select(size_bin, ppn), by = "size_bin"
  )

saveRDS(diff_quantile %>% mutate(dataset = "standard"),
        here::here("data", "rec", "original_size_selectivity.rds"))


sel_bean_size <- ggplot() +
  geom_pointrange(
    data = diff_quantile,
    aes(x = med_dif, xmin = lo_dif, xmax = up_dif, y = size_bin,
        fill = ppn),
    shape = 21
  ) +
  scale_fill_continuous(
    name = tr("Proportion of\nFishery Samples\nin Western Strata", "Proportion des\néchantillons de pêche\ndans les strates occidentales"),
    trans = "sqrt"
  ) +
  labs(x = tr("Difference Between Observed and Predicted Composition", "Différence entre la composition observée et prédite"),
       y = tr("Size Bin (cm)", "Classe de taille (cm)"),
       fill = tr("Proportion of Fishery Samples in Western Strata", "Proportion des échantillons de pêche dans les strates occidentales")) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top"
        ) +
  geom_vline(xintercept = 0, lty = 2)


png(
  fig_path("selectivity/selectivity_bean_size.png"),
  height = 3.25, width = 5, units = "in", res = 250
)
sel_bean_size
dev.off()


# png(
#   here::here("figs", "selectivity", "selectivity_bean_size_comp.png"),
#   height = 6.5, width = 5, units = "in", res = 250
# )
# sel_bean2_size
# dev.off()


# 
# # # calculate observed proportions
# obs_ppn_dat_size <- rkw_dat$size %>%
#   group_by(size_bin) %>%
#   summarize(
#     sum_obs = sum(agg_prob)
#   ) %>%
#   ungroup() %>%
#   # # no smallest sizes observed so add manually as zero
#   # rbind(
#   #   .,
#   #   data.frame(
#   #     size_bin = "55-65", sum_obs = 0
#   #   )
#   # ) %>%
#   mutate(
#     obs_ppn = ifelse(sum_obs == "0", 0, sum_obs / sum(sum_obs)),
#     se_obs_ppn = sqrt(obs_ppn * (1 - obs_ppn) / sum(sum_obs)),
#     lo = pmax(0, obs_ppn - (1.96 * se_obs_ppn)),
#     up = pmin(1, obs_ppn + (1.96 * se_obs_ppn))
#   )
# 
# sel_boxplot_size <- ggplot() +
#   geom_boxplot(
#     data = sim_ppn_dat_size,
#     aes(x = size_bin, y = pred_sim_ppn)
#   ) +
#   geom_pointrange(
#     data = obs_ppn_dat_size,
#     aes(x = size_bin, y = obs_ppn, ymin = lo, ymax = up),
#     colour = "red",
#     position = position_nudge(x = 0.1)
#   ) +
#   # geom_text(
#   #   data = p_val_sig_size %>% filter(dataset == "standard"),
#   #   aes(x = size_bin, y = max(obs_ppn_dat_size$up + 0.1)),
#   #   label = "*", size = 7.5
#   # ) +
#   lims(y = c(0, max(obs_ppn_dat_size$up + 0.2))) +
#   labs(y = "Sample Composition") +
#   ggsidekick::theme_sleek() +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# sel_boxplot_size2 <- ggplot() +
#   geom_boxplot(
#     data = sim_ppn_dat_size,
#     aes(x = size_bin, y = pred_sim_ppn, fill = dataset)
#   ) +
#   geom_pointrange(
#     data = obs_ppn_dat_size,
#     aes(x = size_bin, y = obs_ppn, ymin = lo, ymax = up),
#     colour = "red",
#     position = position_nudge(x = 0.1)
#   ) +
#   geom_text(
#     data = p_val_sig_size,
#     aes(x = size_bin, y = max(obs_ppn_dat_size$up + 0.1)),
#     label = "*", size = 7.5
#   ) +
#   lims(y = c(0, max(obs_ppn_dat_size$up + 0.2))) +
#   labs(y = "Sample Composition") +
#   facet_wrap(~dataset, ncol = 1) +
#   scale_fill_manual(values = dataset_pal, name = "Model") +
#   ggsidekick::theme_sleek() +
#   theme(legend.position = "none",
#         axis.title.x = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1))
# 
# png(
#   here::here("figs", "selectivity", "selectivity_boxplot_size.png"),
#   height = 4.5, width = 5.5, units = "in", res = 250
# )
# sel_boxplot_size
# dev.off()
# 
# png(
#   here::here("figs", "selectivity", "selectivity_boxplot_size_comp.png"),
#   height = 6, width = 5.5, units = "in", res = 250
# )
# sel_boxplot_size2
# dev.off()
