### Summer 5.2 Stock Comp
## Generate estimates of composition for SWVI and NWVI (pooling stat areas)
## with particular emphasis on contribution of summer 5.2 stocks.
## July 14, 2023
## To evaluate data quality and model performance generate:
## 1) Weekly (by year and across years) average comp, with sample sizes, for 
## each region and fishery
## 2) Model predicted mean stock comp, integrating both datasets
## 3) Estimates of relative impact in removals using fixed catches derived from
## historical commercial catches (ppn of August catch caught in W1/2 vs W3/4 and
## in NWVI vs SWVI)

library(tidyverse)
library(brms)
library(tidybayes)


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
dd$week_z <- scale(dd$week) %>% as.numeric()


week_comp <- ggplot(dd, 
       aes(x = week, y = prob, size = samp_nn, fill = year)) +
  geom_jitter(alpha = 0.5, shape = 21) +
  scale_fill_discrete(guide = "none") +
  scale_size_continuous(name = "Weekly\nSample\nSize") +
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

dd_aug <- dd %>% 
  filter(week > 30 & week < 35)
mean(dd_aug$prob)
nrow(dd_aug[dd_aug$prob == 0, ]) / nrow(dd_aug)


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


## COMMERCIAL CATCH DATA -------------------------------------------------------

catch <- read.csv(here::here("data", "comm", "area_g_catch.csv")) %>% 
  janitor::clean_names() %>% 
  mutate(
    region = case_when(
      mgmt_area %in% c("123", "124") ~ "SWVI",
      mgmt_area %in% c("125", "126", "127") ~ "NWVI",
      TRUE ~ "other"
    ),
    date = as.POSIXct(fishing_date,
                         format = "%m/%d/%Y %H:%M"
    ),
    month = lubridate::month(date),
    week = lubridate::week(date),
    day = lubridate::day(date),
    year = lubridate::year(date),
    period = ifelse(day < 16, "early", "late")
  ) %>% 
  filter(!opng_cat == "Exploratory",
         targets_chinook == "YES",
         !region == "other",
         !comments == "ESTIMATES NOT AVAILABLE") 

annual_catch <- catch %>% 
  group_by(week, region, year) %>% 
  summarize(
    total_catch = sum(chinook_kept),
    total_effort = sum(vessels_op),
    total_cpue = total_catch / total_effort,
    .groups = "drop"
  ) %>% 
  distinct() 

ggplot(annual_catch) +
  geom_boxplot(aes(x = as.factor(week), y = total_cpue)) 


aug_catch <- catch %>% 
  filter(month == "8") %>% 
  group_by(year) %>%
  mutate(
    subset = paste(week, region, sep = "_"),
    total_aug_catch = sum(chinook_kept)
  ) %>%
  group_by(week, region, subset, year, total_aug_catch) %>% 
  summarize(
    total_sub_catch = sum(chinook_kept),
    ppn_catch = total_sub_catch / total_aug_catch,
    .groups = "drop"
  ) %>% 
  distinct()


# visualize ppn
pdf(here::here("figs", "summer_impacts", "catch_ppn.pdf"))
ggplot(aug_catch, aes(fill = subset, y = ppn_catch, x = year)) + 
  geom_bar(position = "fill", stat = "identity") +
  labs(y = "Proportion of Total August Catch", x = "Year")
dev.off()


# no late catch prior to 2019 so exclude those years and calculate proportion 
# over time period to sum to one
# NOTE: eventually should be modeled
week_ppn <- catch %>% 
  filter(month == "8",
         year > 2018) %>% 
  # group_by(year) %>%
  mutate(
    week = ifelse(week == 35, 34, week),
    subset = paste(week, region, sep = "_"),
    total_aug_catch = sum(chinook_kept),
    #consolidate week 34 and 35 to keep balanced
  ) %>%
  group_by(week, region, subset, total_aug_catch) %>% 
  summarize(
    total_sub_catch = sum(chinook_kept),
    ppn_catch = total_sub_catch / total_aug_catch,
    .groups = "drop"
  ) %>% 
  distinct()


## MODEL FIT -------------------------------------------------------------------

fit4 <- brm(
  bf(prob ~ region*week_z + region*I(week_z^2) + dataset2 + (1 | year),
     phi ~ dataset2 + samp_nn,
     zi ~ region*week_z + region*I(week_z^2) + dataset2 + (1 | year)
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
saveRDS(fit4, 
        here::here("data", "model_fits", "summer_impacts", "fit4.rds"))

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


# m1 <- glmmTMB::glmmTMB(
#   prob2 ~ region * week_z + region * I(week_z^2) + (1 | year), 
#   data = dd, 
#   family= glmmTMB::beta_family()
# )
# 
# library(DHARMa)
# res <- simulateResiduals(m1)
# plot(res, rank = T)



## PREDICTIONS -----------------------------------------------------------------

new_data <- expand.grid(
  region = unique(dd$region),
  dataset2 = unique(dd$dataset2),
  week = seq(min(dd$week), max(dd$week), by = 0.2),
  year = "2021",
  # minor imapct of sample size
  samp_nn = 50
)
new_data$week_z <- (new_data$week - mean(dd$week)) / sd(dd$week) 

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
    units = "in", height = 3, width = 5.5)
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
  labs(y = "Predicted Proportion\nSummer 5_2 in Catch", 
       x = "Week") +
  ggsidekick::theme_sleek()
dev.off()


## convert to harvested fish in management scenarios

# tac (max catch in recent years was 25k)
tac <- 30000
esc_forecast <- 25000


# ppns as estimated above
ppn_key <- expand.grid(
  week = c(31, 32, 33, 34),
  region = c("SWVI", "NWVI")
) %>% 
  mutate(
    week_group = ifelse(week < 33, "early", "late")
  ) %>% 
  left_join(
    ., 
    week_ppn %>% select(-c(subset:total_sub_catch)), 
    by = c("week", "region")
  )


scen_tbl <- expand.grid(
  scenario = c("status_quo", "2023", "adj_area", "adj_time"),
  week = seq(31, 34, by = 1),
  region = unique(dd$region),
  dataset2 = "comm",
  samp_nn = 100,
  year = "2021"
) %>% 
  left_join(., ppn_key, by = c("week", "region")) %>% 
  mutate(
    week_z = (week - mean(dd$week)) / sd(dd$week),
    # define management scenarios and calculate catch
    ppn_catch_adj = case_when(
      scenario == "2023" & week_group == "early" ~ 0,
      scenario == "adj_area" & week_group == "early" & region == "SWVI" ~ 0,
      scenario == "adj_time" & week == "31" ~ 0,
      TRUE ~ ppn_catch
    ),
    catch_adj = ppn_catch_adj * tac
  ) %>% 
  as_tibble() %>% 
  group_by(scenario) %>% 
  mutate(total_catch = sum(catch_adj)) %>% 
  group_by(scenario, total_catch) %>% 
  group_nest()

scen_tbl$preds <- purrr::map(
  scen_tbl$data,
  ~ epred_draws(fit4, newdata = .x, re_formula = NA, scale = "response") 
)
# scen_tbl$preds <- purrr::map(
#   scen_tbl$data,
#   ~ epred_draws(fit4, newdata = .x, re_formula = NULL, scale = "response") 
# )

# exploitation rate preds
pred_dat <- scen_tbl$preds[[1]] %>% 
  group_by(week, region, .draw) %>% 
  summarize(
    exp = mean(.epred), .groups = "drop"
  ) %>%
  group_by(week, region) %>% 
  summarize(
    mean_exp = mean(exp),
    lo_exp = quantile(exp, 0.05),
    up_exp = quantile(exp, 0.95)
  ) 

mean_stock_comp <- ggplot(pred_dat, aes(x = week, y = mean_exp)) +
  geom_pointrange(aes(ymin = lo_exp, ymax = up_exp)) +
  facet_wrap(~region) +
  labs(y = "Predicted Commercial\nSummer 5_2 Composition", x = "August Period") +
  ggsidekick::theme_sleek()

# catch preds
pred_catch <- scen_tbl %>%
  # select(-data) %>%
  unnest(cols = c(preds)) %>%
  group_by(scenario, region, week, catch_adj, .draw) %>% 
  summarize(
    exp = mean(.epred), .groups = "drop"
  ) %>% 
  mutate(
    catch_52 = catch_adj * exp
  )

# catch by region and week 
pred_catch %>% 
  group_by(scenario, region, week) %>%
  summarize(
    mean_catch = mean(catch_52),
    lo_catch = quantile(catch_52, 0.05),
    up_catch = quantile(catch_52, 0.95)
  ) %>% 
  ggplot(., aes(x = week, y = mean_catch)) +
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
    mean_er = mean_catch / esc_forecast,
    lo_er = lo_catch / esc_forecast,
    up_er = up_catch / esc_forecast
  )

pred_catch_est <- ggplot(pred_catch_est_dat, 
                         aes(x = scenario, y = mean_catch)) +
  geom_pointrange(aes(ymin = lo_catch, ymax = up_catch)) +
  labs(y = "Predicted Summer 5_2/nCommercial Catch (pieces)", 
       x = "Area G Region") +
  ggsidekick::theme_sleek()

total_catch_dat <- scen_tbl %>% 
  select(scenario, total_catch) %>% 
  mutate(scenario = fct_recode(
    scenario, "Current 2023" = "2023", "Adjusted Area" = "adj_area", 
    "Adjusted Time" = "adj_time"
  ))

# difference in total escapement
kobe_plot <- pred_catch %>% 
  group_by(scenario, .draw) %>% 
  summarize(
    total_catch52 = sum(catch_52)
  ) %>%
  pivot_wider(names_from = scenario, values_from = total_catch52, 
              names_prefix = "scen_") %>% 
  mutate(
    `Current 2023` = scen_status_quo - scen_2023,
    `Adjusted Area` = scen_status_quo - scen_adj_area,
    `Adjusted Time` = scen_status_quo - scen_adj_time
  ) %>% 
  select(.draw, `Current 2023`, `Adjusted Area`, `Adjusted Time`) %>% 
  pivot_longer(cols = c(`Current 2023`, `Adjusted Area`, `Adjusted Time`),
               names_to = "scenario") %>% 
  group_by(scenario) %>% 
  summarize(
    mean_diff = mean(value) / esc_forecast,
    lo_diff = quantile(value, 0.05) / esc_forecast,
    up_diff = quantile(value, 0.95) / esc_forecast
  ) %>% 
  left_join(., total_catch_dat, by = "scenario") %>% 
  mutate(diff_catch = total_catch / tac) %>% 
  ggplot(., aes(y = mean_diff, x = diff_catch)) +
  geom_pointrange(aes(ymin = lo_diff, ymax = up_diff, fill = scenario),
                  shape = 21) +
  labs(y = "Difference in Predicted\nEscapement", 
       x = "Predicted Proportion of TAC Caught") +
  ggsidekick::theme_sleek()


png(here::here("figs", "summer_impacts", "kobe_plot.png"), res = 250,
    units = "in", height = 3.5, width = 5.5)
kobe_plot
dev.off()

pdf(here::here("figs", "summer_impacts", "summary_figs.pdf"))
annual_comp
obs_comp
pred_ribbon
mean_stock_comp
pred_catch_est
diff_catch
dev.off()