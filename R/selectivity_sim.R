## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)
library(mvtweedie)

# modified tweedie prediction that allows for exclude = list()
source(here::here("R", "functions", "pred_mvtweedie2.R"))
source(here::here("R", "functions", "sim_multinomial.R"))



## STOCK SELECTIVITY -----------------------------------------------------------

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


# generate simulations for each model fit
future::plan(future::multisession, workers = 6)
pred_tbl$sim_dat <- furrr::future_map(
  pred_tbl$pred_dat, ~ sim_foo(.x, nsim = 500),
  .options = furrr::furrr_options(seed = TRUE)
)


# calculate p-values for each and print
p_val <- purrr::map2(
  pred_tbl$sim_dat, pred_tbl$dataset,
  ~ .x %>% 
    group_by(stock_group) %>% 
    summarize(p_value = sum(test_stat) / length(unique(sim_i))) %>% 
    mutate(dataset = .y)
  ) %>% 
  bind_rows()
p_val_sig <- p_val %>% 
  filter(p_value < 0.05) 


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
    obs_ppn = sum_obs / sum(sum_obs),
    se_obs_ppn = sqrt(obs_ppn * (1 - obs_ppn) / sum(sum_obs)),
    lo = pmax(0, obs_ppn - (1.96 * se_obs_ppn)),
    up = pmin(1, obs_ppn + (1.96 * se_obs_ppn))
  )

dataset_pal <- c("#e0f3db", "#a8ddb5", 
                 "#43a2ca")
names(dataset_pal) <- levels(sim_ppn_dat$dataset)

sel_boxplot <- ggplot() +
  geom_boxplot(
    data = sim_ppn_dat,
    aes(x = stock_group, y = pred_sim_ppn, fill = dataset)
    ) +
  geom_pointrange(
    data = obs_ppn_dat, 
    aes(x = stock_group, y = obs_ppn, ymin = lo, ymax = up), 
    colour = "red",
    position = position_nudge(x = 0.1)
    ) +
  geom_text(
    data = p_val_sig,
    aes(x = stock_group, y = max(sim_ppn_dat$pred_sim_ppn + 0.1)), 
    label = "*", size = 7.5
  ) +
  lims(y = c(0, max(sim_ppn_dat$pred_sim_ppn + 0.2))) +
  labs(y = "Simulated Sample Composition") +
  facet_wrap(~dataset, ncol = 1) +
  scale_fill_manual(values = dataset_pal, name = "Model") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

png(
  here::here("figs", "selectivity", "selectivity_boxplot_stock.png"),
  height = 7.5, width = 7.5, units = "in", res = 250
)
sel_boxplot
dev.off()


## SIZE SELECTIVITY ------------------------------------------------------------


rkw_dat_size <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat_size.rds")
) %>%
  rename(
    year_n = year
  ) %>%
  filter(
    era == "current"
  )


# import fitted models
fit_size <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_size_mvtw.rds"
  )
)
slot_fit_size <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_size_slot.rds"
  )
)

fit_list <- list(fit_size, slot_fit_size)


# ensure size bins the same
identical(
  levels(fit_size$model$size_bin), levels(slot_fit_size$model$size_bin))
size_bins <- levels(fit_size$model$size_bin)


# infill rkw data to ensure each stock group present for each sampling event
dat_list <- split(rkw_dat_size, rkw_dat_size$sample_id)

new_dat <- purrr::map(
  dat_list,
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
        rkw_dat_size %>% 
          select(sample_id, size_bin, obs_ppn = agg_prob), 
        by = c("sample_id", "size_bin")
      ) %>% 
      mutate(
        obs_ppn = ifelse(is.na(obs_ppn), 0, obs_ppn),
        sg_year = paste(size_bin, year_n, sep = "_") %>% as.factor(),
        slot_limit = ifelse(
          year_n > 2018, "yes", "no"
        )
      )
  }
) %>% 
  bind_rows()


# generate predictions based on each fitted model and store as tibble
pred_tbl <- tibble(
  dataset = c("standard", "management")
) %>% 
  mutate(
    pred_dat = purrr::map(
      fit_list, function (x) {
        preds <- pred_dummy(
          x, se.fit = TRUE, category_name = "size_bin", origdata = x$model,
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
p_val <- purrr::map2(
  pred_tbl$sim_dat, pred_tbl$dataset,
  ~ .x %>% 
    group_by(stock_group) %>% 
    summarize(p_value = sum(test_stat) / length(unique(sim_i))) %>% 
    mutate(dataset = .y)
) %>% 
  bind_rows()
p_val_sig <- p_val %>% 
  filter(p_value < 0.05) 


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
    obs_ppn = sum_obs / sum(sum_obs),
    se_obs_ppn = sqrt(obs_ppn * (1 - obs_ppn) / sum(sum_obs)),
    lo = pmax(0, obs_ppn - (1.96 * se_obs_ppn)),
    up = pmin(1, obs_ppn + (1.96 * se_obs_ppn))
  )

dataset_pal <- c("#e0f3db", "#a8ddb5", 
                 "#43a2ca")
names(dataset_pal) <- levels(sim_ppn_dat$dataset)

sel_boxplot <- ggplot() +
  geom_boxplot(
    data = sim_ppn_dat,
    aes(x = stock_group, y = pred_sim_ppn, fill = dataset)
  ) +
  geom_pointrange(
    data = obs_ppn_dat, 
    aes(x = stock_group, y = obs_ppn, ymin = lo, ymax = up), 
    colour = "red",
    position = position_nudge(x = 0.1)
  ) +
  geom_text(
    data = p_val_sig,
    aes(x = stock_group, y = max(sim_ppn_dat$pred_sim_ppn + 0.1)), 
    label = "*", size = 7.5
  ) +
  lims(y = c(0, max(sim_ppn_dat$pred_sim_ppn + 0.2))) +
  labs(y = "Simulated Sample Composition") +
  facet_wrap(~dataset, ncol = 1) +
  scale_fill_manual(values = dataset_pal, name = "Model") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1))

png(
  here::here("figs", "selectivity", "selectivity_boxplot.png"),
  height = 7.5, width = 7.5, units = "in", res = 250
)
sel_boxplot
dev.off()
