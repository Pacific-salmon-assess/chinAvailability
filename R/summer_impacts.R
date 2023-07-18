### Summer 5.2 Stock Comp
## Generate estimates of composition for SWVI and NWVI (pooling stat areas)
## with particular emphasis on contribution of summer 5.2 stocks.
## July 14, 2023
## To evaluate data quality and model performance generate:
## 1) Weekly (by year and across years) average comp, with sample sizes, for 
## each region and fishery
## 2) Model predicted mean stock comp, integrating both datasets
## 3) Estimates of relative impact in removals using fixed catches 

library(tidyverse)
library(brms)


## RECREATIONAL DATA -----------------------------------------------------------


comp_in_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  mutate(week = lubridate::week(date))

comp_in_raw %>% 
  filter(Region1Name == "Fraser_Summer_5.2") %>% 
  select(stock, Region1Name, Region1Name) %>% 
  distinct() %>% 
  arrange(stock, Region1Name) %>% 
  print(n = Inf)


rec_dat <- comp_in_raw %>% 
  filter(cap_region %in% c("NWVI", "SWVI"),
         area > 100,
         fl >= 550#,
         # month %in% c("June", "July", "August", "September")
         ) %>% 
  mutate(
    smu = ifelse(Region1Name == "Fraser_Summer_5.2", Region1Name, "other"),
    sample_id = paste(week, year, cap_region, "rec", sep = "_"),
    year = as.factor(year),
    dataset = "rec"
  ) %>% 
  select(
    dataset, fish_id = id, sample_id, week, year, region = cap_region, area, 
    smu, prob
  )



## COMMERCIAL DATA -------------------------------------------------------------

stock_key <- readRDS(here::here("data", "rec", "finalStockList_Jul2023.rds"))

comm_dat <- readRDS(here::here("data", "comm", "wcviIndProbsLong.rds")) %>% 
  select(-c(Region1Name:pst_agg)) %>% 
  #adjust stock IDs to make sure correct
  left_join(., stock_key, by = "stock") %>% 
  filter(as.numeric(as.character(area)) > 100#,
         # month %in% c("6", "7", "8", "9")
         ) %>% 
  mutate(
    smu = ifelse(Region1Name == "Fraser_Summer_5.2", Region1Name, "other"),
    sample_id = paste(week, year, region, "comm", sep = "_"),
    dataset = "comm"
  ) %>% 
  select(
    dataset, fish_id = id, sample_id, week, year, region, area, smu,
    prob = adj_prob
  )


## TAGGING DATA ----------------------------------------------------------------


set_dat <- readRDS(here::here("data", "tagging", "cleanSetData.RDS")) %>% 
  mutate(nearshore = ifelse(coast_dist_km < 9.26, TRUE, FALSE))

tag_dat <- readRDS(here::here("data", "tagging", "clean_catch.RDS")) %>% 
  left_join(., set_dat %>% select(event, nearshore), by = "event") %>% 
  mutate(
    mu = ifelse(
      cu %in% c("MFR-summer", "NTh-sum", "STh-1.3"), 
      "Fraser_Summer_5.2",
      "other"
    ),
    week = lubridate::week(deployment_time),
    region = "SWVI",
    prob = 1,
    area = 123
    ) %>% 
  filter(!is.na(mu),
         nearshore == FALSE#,
         # month %in% c("6", "7", "8", "9")
         ) %>% 
  mutate(
    sample_id = paste(week, year, region, "tag", sep = "_"),
    dataset = "tag"
  ) %>% 
  select(
    dataset, fish_id = fish, sample_id, week, year, region, area, smu = mu,
    prob
  )


## COMBINE AND EXPLORE ---------------------------------------------------------

all_dat <- rbind(comm_dat, rec_dat, tag_dat)

dumm <- all_dat %>% select(fish_id, year) %>% distinct() %>% group_by(year) %>% tally
sum(dumm$n)

all_dat %>% 
  filter(!smu == "other") %>% 
  pull(prob) %>% 
  sum()

samp_size <- all_dat %>% 
  group_by(sample_id) %>% 
  summarize(samp_nn = length(unique(fish_id))) 

expand_dat <- expand.grid(
  sample_id = unique(all_dat$sample_id),
  smu = unique(all_dat$smu)
) %>% 
  left_join(
    ., 
    all_dat %>% 
      select(dataset, sample_id, week, year, region) %>% 
      distinct(),
    by = "sample_id"
  ) %>% 
  left_join(., samp_size, by = "sample_id")

all_samps <- all_dat %>% 
  group_by(sample_id) %>% 
  mutate(samp_nn = length(unique(fish_id))) %>% 
  group_by(dataset, sample_id, samp_nn, week, year, region, smu) %>% 
  summarize(smu_prob = sum(prob),
            smu_ppn = smu_prob / samp_nn,
            .groups = "drop") %>% 
  distinct() %>% 
  rename(prob = smu_ppn) 

all_samps2 <- expand_dat %>% 
  left_join(., 
            all_samps %>%
              select(sample_id, smu, prob), 
            by = c("sample_id", "smu")) %>% 
  mutate(
    prob = replace_na(prob, 0)
  ) 

dd <- all_samps2 %>% filter(smu == "Fraser_Summer_5.2")
dd$dataset2 <- ifelse(dd$dataset == "tag", "comm", dd$dataset)

week_comp <- ggplot(dd, 
       aes(x = week, y = prob, size = samp_nn, fill = year)) +
  geom_jitter(alpha = 0.5, shape = 21) +
  facet_grid(dataset2 ~ region) +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Catch", x = "Week")


png(here::here("figs", "summer_impacts", "weekly_comp.png"), res = 250,
    units = "in", height = 4.5, width = 5.5)
week_comp
dev.off()

# pool years
ggplot(dd, aes(x = as.factor(week), y = prob)) +
  geom_boxplot(alpha = 0.5, shape = 21) +
  facet_grid(dataset2 ~ region) +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Sample", x = "Week")


## annual means (i.e. sum all probs within a year then divide)
# all_dat_yr <- rbind(comm_dat, rec_dat) %>% 
#   filter(week > 30 & week < 35) %>% 
#   mutate(sample_id = paste(year, region, dataset, sep = "_")) 
#   
# samp_size_yr <- all_dat %>% 
#   group_by(sample_id) %>% 
#   summarize(samp_nn = length(unique(fish_id))) 
# 
# expand_dat_yr <- expand.grid(
#   sample_id = unique(all_dat_yr$sample_id),
#   smu = unique(all_dat_yr$smu)
# ) %>% 
#   left_join(
#     ., 
#     all_dat_yr %>% 
#       select(dataset, sample_id, year, region) %>% 
#       distinct(),
#     by = "sample_id"
#   ) %>% 
#   left_join(., samp_size_yr, by = "sample_id")
# 
# all_samps_yr <- all_dat_yr %>% 
#   group_by(sample_id) %>% 
#   mutate(samp_nn = length(unique(fish_id))) %>% 
#   group_by(dataset, sample_id, samp_nn, year, region, smu) %>% 
#   reframe(smu_prob = sum(prob),
#           smu_ppn = smu_prob / samp_nn,
#           .groups = "drop") %>% 
#   distinct() %>% 
#   rename(prob = smu_ppn) 
# 
# all_samps2_yr <- expand_dat_yr %>% 
#   left_join(., 
#             all_samps_yr %>%
#               select(sample_id, smu, prob), 
#             by = c("sample_id", "smu")) %>% 
#   mutate(
#     prob = replace_na(prob, 0)
#   ) 
# 
# annual_mean <- all_samps2_yr %>% 
#   filter(smu == "Fraser_Summer_5.2") %>% 
#   group_by(dataset, region, year) %>% 
#   summarize(mean_prob = mean(prob))
# 
# annual_comp <- ggplot(annual_mean, aes(x = region, y = mean_prob)) +
#   geom_boxplot() +
#   facet_wrap(~dataset) +
#   ggsidekick::theme_sleek() +
#   labs(y = "Observed Proportion of\nSummer 5_2 in August Catch", 
#        x = "Area G Region")
# 
# png(here::here("figs", "summer_impacts", "annual_comp.png"), res = 250,
#     units = "in", height = 5.5, width = 5.5)
# annual_comp
# dev.off()



## MODEL FIT -------------------------------------------------------------------

library(brms)
library(tidybayes)

dd$prob2 <- dd$prob + 0.000001

beta_priors <- c(
  prior(normal(-4.8, 2.5), class = "Intercept"),
  prior(gamma(0.01, 0.01), class = "phi")
)

fit3 <- brm(
  bf(prob2 ~ region + week + I(week^2) + dataset2 + (1 | year),
     phi ~ dataset2
  ),
  prior = c(
    prior(normal(-4.8, 2.5), class = "Intercept")#,
    # prior(normal(0, 5), class = "b")
  ),
  data = dd,
  family = Beta(),
  init = "0",
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234
)


# if convergence issues found can remove week interaction and samp_nn pars
fit4 <- brm(
  bf(prob ~ region*week + region*I(week^2) + dataset2 + (1 | year),
     phi ~ dataset2 + samp_nn,
     zi ~ region*week + region*I(week^2) + dataset2 + (1 | year)
  ),
  prior = c(
    prior(normal(-4.8, 2.5), class = "Intercept", dpar = ""),
    prior(normal(0, 5), class = "b")
  ),
  # control = list(#adapt_delta = 0.91,
  #                max_treedepth = 11),
  data = dd,
  init = "0",
  family = zero_inflated_beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234
)

plot(fit4)


# check ppn zeros
nrow(dd[dd$prob == "0", ]) / nrow(dd)
obs_preds <- dd %>%
  add_predicted_draws(fit4) 
nrow(obs_preds[obs_preds$.prediction == 0, ]) / nrow(obs_preds)

pp_check(fit4, "dens_overlay_grouped", ndraws = 1000, group = "dataset2") +
  lims(y = c(0, 100))
pp_check(fit4, "dens_overlay", ndraws = 1000) +
  lims(y = c(0, 100))
pp_check(fit4, type = "scatter_avg", ndraws = 100)
pp_check(fit4, type = "stat_2d")

loo1 <- loo(fit4, save_psis = TRUE, cores = 4)
psis1 <- loo1$psis_object
lw <- weights(psis1)

bayesplot::ppc_loo_pit_overlay(dd$prob, posterior_predict(fit4), lw = lw)
bayesplot::ppc_loo_pit_qq(dd$prob, posterior_predict(fit4), lw = lw)


new_data <- expand.grid(
  region = unique(dd$region),
  dataset2 = unique(dd$dataset2),
  week = seq(min(dd$week), max(dd$week), by = 0.2),
  year = "2021",
  # minor imapct of sample size
  samp_nn = 50
) 

preds_zi <- fit4 %>% 
  predicted_draws(newdata = new_data) %>% 
  mutate(is_zero = .prediction == 0,
         .prediction = ifelse(is_zero, .prediction - 0.01, .prediction))

ggplot(preds_zi, aes(x = .prediction)) +
  geom_histogram(aes(fill = is_zero), binwidth = 0.025, 
                 boundary = 0, color = "white") +
  facet_grid(dataset2~region)


preds4 <- fit4 %>% 
  epred_draws(newdata = new_data, re_formula = NA, scale = "response")
# assume 21 year effects
preds4_21 <- fit4 %>% 
  epred_draws(newdata = new_data, re_formula = NULL, scale = "response")


preds4_mean <- preds4 %>% 
  group_by(region, dataset2, week) %>% 
  summarize(.epred = mean(.epred), .groups = "drop")

pred_line <- ggplot() +
  geom_line(data = preds4 %>% 
              filter(.draw %in% sample.int(200)), 
            aes(x = week, y = .epred, group = paste(dataset2, region, .draw)),
            alpha = 0.2
            ) +
  geom_line(data = preds4_21 %>% 
              filter(.draw %in% sample.int(200)), 
            aes(x = week, y = .epred, group = paste(dataset2, region, .draw)),
            alpha = 0.2, colour = "red"
  ) +
  facet_grid(dataset2~region) +
  labs(y = "Proportion Summer 5_2 in Catch", x = "Week") +
  ggsidekick::theme_sleek()

obs_comp <- ggplot() +
  geom_line(data = preds4_mean, aes(x = week, y = .epred), 
            color = "red", lwd = 2, group = 1) +
  geom_jitter(data = dd, aes(x = week, y = prob, size = samp_nn, fill = year), 
              alpha = 0.5, shape = 21) +
  facet_grid(dataset2~region) +
  labs(y = "Proportion Summer 5_2 in Catch", x = "Week") +
  ggsidekick::theme_sleek()

pred_ribbon <- ggplot(data = preds4, aes(x = week, y = .epred)) +
  stat_lineribbon() +
  facet_grid(dataset2~region) + 
  labs(y = "Predicted Proportion Summer 5_2 in Catch", 
       x = "Week") +
  ggsidekick::theme_sleek()

preds4_mean %>% 
  group_by(region, dataset2) %>% 
  mutate(max_pred = max(.epred)) %>% 
  filter(.epred == max_pred)


png(here::here("figs", "summer_impacts", "comm_preds.png"), res = 250,
    units = "in", height = 3.5, width = 5.5)
preds4 %>% 
  filter(dataset2 == "comm") %>% 
  group_by(week, region, dataset2) %>% 
  summarize(med_pred = median(.epred),
            lo_pred = quantile(.epred, 0.05),
            up_pred = quantile(.epred, 0.95)) %>% 
  ggplot(data = ,
       aes(x = week, y = med_pred)) +
  geom_line() +
  geom_ribbon(aes(ymin = lo_pred, ymax = up_pred), alpha = 0.2) +
  facet_wrap(~region) + 
  labs(y = "Predicted Proportion Summer 5_2 in Catch", 
       x = "Week") +
  ggsidekick::theme_sleek()
dev.off()


## convert to harvested fish in man scenarios
scen_tbl <- expand.grid(
  scenario = c("status_quo", "2023", "adj"),
  week = seq(31, 34, by = 1),
  region = unique(dd$region),
  dataset2 = "comm",
  samp_nn = 100,
  year = "2021"
) %>% 
  mutate(
    week_group = ifelse(week %in% c(31, 32), "early", "late"),
    catch = case_when(
      scenario == "2023" & week_group == "early" ~ 0,
      scenario == "adj" & region == "SWVI" & week_group == "early" ~ 0,
      TRUE ~ 10000
    )
  ) %>% 
  as_tibble() %>% 
  group_by(scenario) %>% 
  group_nest()

scen_tbl$preds <- purrr::map(
  scen_tbl$data,
  ~ epred_draws(fit4, newdata = .x, re_formula = NA, scale = "response") %>% 
    mutate(catch_52 = catch * .epred)
)
# scen_tbl$preds <- purrr::map(
#   scen_tbl$data,
#   ~ epred_draws(fit4, newdata = .x, re_formula = NULL, scale = "response") %>% 
#     mutate(catch_52 = catch * .epred)
# )

# exploitation rate preds
pred_dat <- scen_tbl$preds[[1]] %>% 
  group_by(week_group, region, .draw) %>% 
  summarize(
    exp = mean(.epred), .groups = "drop"
  ) %>%
  group_by(week_group, region) %>% 
  summarize(
    mean_exp = mean(exp),
    lo_exp = quantile(exp, 0.05),
    up_exp = quantile(exp, 0.95)
  ) 

mean_stock_comp <- ggplot(pred_dat, aes(x = week_group, y = mean_exp)) +
  geom_pointrange(aes(ymin = lo_exp, ymax = up_exp)) +
  facet_wrap(~region) +
  labs(y = "Predicted Commercial\nSummer 5_2 Composition", x = "August Period") +
  ggsidekick::theme_sleek()

# catch preds
pred_catch <- scen_tbl %>%
  # select(-data) %>%
  unnest(cols = c(preds)) %>%
  group_by(scenario, region, week_group, catch, .draw) %>% 
  summarize(
    exp = mean(.epred), .groups = "drop"
  ) %>% 
  mutate(
    catch_52 = catch * exp
  )

# catch by region and week group
pred_catch %>% 
  group_by(scenario, region, week_group, catch) %>%
  summarize(
    mean_catch = mean(catch_52),
    lo_catch = quantile(catch_52, 0.05),
    up_catch = quantile(catch_52, 0.95)
  ) %>% 
  ggplot(., aes(x = week_group, y = mean_catch)) +
  geom_pointrange(aes(ymin = lo_catch, ymax = up_catch)) +
  facet_grid(scenario ~ region) +
  ggsidekick::theme_sleek()

# catch by region
pred_catch_est_dat <- pred_catch %>% 
  group_by(scenario, .draw) %>% 
  summarize(
    month_catch = sum(catch_52)
  ) %>%
  group_by(scenario) %>% 
  summarize(
    mean_catch = mean(month_catch),
    lo_catch = quantile(month_catch, 0.05),
    up_catch = quantile(month_catch, 0.95),
    mean_er = mean_catch / 25000,
    lo_er = lo_catch / 25000,
    up_er = up_catch / 25000
  )

pred_catch_est <- ggplot(pred_catch_est_dat, aes(x = scenario, y = mean_catch)) +
  geom_pointrange(aes(ymin = lo_catch, ymax = up_catch)) +
  labs(y = "Predicted Summer 5_2 Commercial Catch (pieces)", x = "Area G Region") +
  ggsidekick::theme_sleek()


# difference in total catch
diff_catch <- pred_catch %>% 
  group_by(scenario, .draw) %>% 
  summarize(
    open_catch = sum(catch_52)
  ) %>%
  pivot_wider(names_from = scenario, values_from = open_catch, 
              names_prefix = "scen_") %>% 
  mutate(
    `New to Status Quo` = scen_status_quo - scen_2023,
    `New to Adjusted` = scen_adj - scen_2023
  ) %>% 
  select(.draw, `New to Adjusted`, `New to Status Quo`) %>% 
  pivot_longer(cols = c(`New to Adjusted`, `New to Status Quo`)) %>% 
  group_by(name) %>% 
  summarize(
    mean_diff = mean(value),
    lo_diff = quantile(value, 0.05),
    up_diff = quantile(value, 0.95)
  ) %>% 
  ggplot(., aes(y = mean_diff, x = name)) +
  geom_pointrange(aes(ymin = lo_diff, ymax = up_diff),
                  fill = "white", shape = 21) +
  labs(y = "Difference in Predicted\nSummer 5_2 Catch (pieces)", x = "Scenario") +
  ggsidekick::theme_sleek()


png(here::here("figs", "summer_impacts", "catch_effects.png"), res = 250,
    units = "in", height = 3.5, width = 5.5)
diff_catch
dev.off()

pdf(here::here("figs", "summer_impacts", "summary_figs.pdf"))
annual_comp
obs_comp
pred_ribbon
mean_stock_comp
pred_catch_est
diff_catch
dev.off()