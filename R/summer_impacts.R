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
library(stockseasonr)


## RECREATIONAL DATA -----------------------------------------------------------


comp_in_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  mutate(week = lubridate::week(date))

comp_in_raw %>% 
  filter(Region1Name == "Fraser_Summer_5.2") %>% 
  select(stock, Region1Name, Region1Name) %>% 
  distinct() %>% 
  arrange(stock, Region1Name) %>% 
  print(n = Inf)

comp_in_raw %>% 
  filter(cap_region %in% c("NWVI", "SWVI")) %>% 
  select(area, cap_region) %>% 
  distinct()

rec_dat <- comp_in_raw %>% 
  filter(cap_region %in% c("NWVI", "SWVI"),
         area > 100,
         fl >= 550,
         month %in% c("June", "July", "August", "September")) %>% 
  mutate(
    smu = ifelse(Region1Name == "Fraser_Fall", Region1Name, "other"),
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
  filter(as.numeric(as.character(area)) > 100,
         month %in% c("6", "7", "8", "9")) %>% 
  mutate(
    smu = ifelse(Region1Name == "Fraser_Fall", Region1Name, "other"),
    sample_id = paste(week, year, region, "comm", sep = "_"),
    dataset = "comm"
  ) %>% 
  select(
    dataset, fish_id = id, sample_id, week, year, region, area, smu,
    prob = adj_prob
  )


## COMBINE AND EXPLORE ---------------------------------------------------------

all_dat <- rbind(comm_dat, rec_dat)


samps <- all_dat %>% 
  group_by(sample_id) %>% 
  mutate(samp_nn = length(unique(fish_id))) 

all_samps <- samps %>% 
  group_by(dataset, sample_id, week, year, region, smu) %>% 
  summarize(smu_prob = sum(prob),
            smu_ppn = smu_prob / samp_nn,
            .groups = "drop") %>% 
  distinct() %>% 
  rename(prob = smu_ppn) %>% 
  select(-smu_prob)

focal_dat <- samps %>% 
  select(dataset, sample_id, week, year, region, samp_nn) %>% 
  distinct() %>% 
  left_join(., 
            all_samps %>%
              filter(smu == "Fraser_Fall") %>% 
              select(sample_id, prob), 
            by = c("sample_id")) %>% 
  mutate(
    prob = replace_na(prob, 0)
  ) 

ggplot(focal_dat, aes(x = week, y = prob, size = samp_nn)) +
  geom_jitter(alpha = 0.5) +
  facet_grid(dataset ~ region) +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Catch", x = "Week (July and August Only)")

ggplot(focal_dat, aes(x = as.factor(week), y = smu_ppn)) +
  geom_boxplot() +
  facet_grid(dataset ~ region) +
  ggsidekick::theme_sleek() +
  labs(y = "Proportion Catch", x = "Week (July and August Only)")

focal_dat %>% 
  filter(week > 30) %>% 
  group_by(cap_region, week) %>% 
  summarize(mean_comp = mean(smu_ppn))


## MODEL FIT -------------------------------------------------------------------

# subset predicted composition dataset
pred_dat_comp <- expand.grid(
  region = unique(focal_dat$region)#,
  # dataset = unique(focal_dat$dataset),
  # week = seq(min(focal_dat$week),
  #               max(focal_dat$week),
  #               by = 0.2
  # )
) 

dd <- all_samps %>% filter(smu == "Fraser_Fall")


model_beta_bayes_1 <- brm(
  bf(prob ~ 1 + region + s(week, k = 4, bs = "tp", by = region) + (1 | year)#,
     # phi ~ quota + polyarchy
     ),
  data = dd,
  family = Beta(),
  chains = 4, iter = 2000, warmup = 1000,
  cores = 4, seed = 1234
)




fit1 <- fit_stockseasonr(
  comp_formula = smu ~ 1 + region,
    # s(week, bs = "tp", k = 4, by = region) + (1 | year),
  comp_dat = all_samps,
  pred_dat = pred_dat_comp ,
  model = "dirichlet",
  # random_walk = TRUE,
  fit = TRUE,
  # nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)
