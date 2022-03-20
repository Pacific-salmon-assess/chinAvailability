### MVN Random Intercepts Dirichlet
## Adapt stockseasonr model with MVN intercepts
## NOTE: attempted to incorporate random smooths, but unclear how to proceed 
## given parameters are a matrix, not vector
## Nov. 29, 2021
## Updated March 19 to explore composition at area levels

library(tidyverse)
library(TMB)


# tmb models
compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))



# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal") %>% 
  mutate(reg = abbreviate(cap_region, minlength = 4),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         sample_id = paste(month_n, reg, yday, year, sep = "_")) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id))
  ) %>% 
  ungroup()
stock_comp <- comp1 %>%  
  group_by(sample_id, reg, reg_c = cap_region, month, month_n, year, nn, 
           pst_agg) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(reg %in% c("JdFS", "SSoG")) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg))


# prediction datasets 
pred_dat_comp1 <- group_split(stock_comp, reg) %>%
  map_dfr(., function(x) {
    expand.grid(
      reg = unique(x$reg),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  mutate(
    reg_month_year = paste(reg, month_n, year, sep = "_"),
    key_var = as.factor(reg_month_year)
  )