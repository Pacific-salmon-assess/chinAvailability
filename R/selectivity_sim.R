## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)
library(mvtweedie)


# import SRKW prey data
rkw_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) %>% 
  # downscale utm coordinates to match model
  mutate(
    utm_y = utm_y / 1000,
    utm_x = utm_x / 1000
    ) %>% 
  rename(
    week_n = week,
    year_n = year
  ) %>% 
  filter(
    era == "current"
  ) 

rkw_dat_pooled <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat_pooled.rds")
) %>% 
  # downscale utm coordinates to match model
  mutate(
    utm_y = utm_y / 1000,
    utm_x = utm_x / 1000
  ) %>% 
  rename(
    week_n = week
  ) %>% 
  filter(
    era == "current"
  ) 


# import original fishery data used to fit model
agg_dat <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
)


# import fitted model
model_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
  )
)

# identify unique stock groups within model (since some may be absent from
# prey dataset)
stks <- levels(model_fit$model$stock_group)


dum <- rkw_dat %>% 
  filter(sample_id == rkw_dat$sample_id[1])


# preds 

# make new data frame consistent with SRKW diet to generate predictions from
newdata <- expand.grid(
  week_n = unique(dum$week_n),
  stock_group = stks,
  year_n = unique(dum$year_n)
) %>% 
  left_join(
    ., 
    dum %>% 
      select(week_n, year_n, utm_x, utm_y) %>% 
      distinct(),
    by = c("week_n", "year_n")
  ) 

preds <- predict(
  model_fit,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata#,
  # exclude = yr_coefs
)

new_dat <- cbind(
  newdata, fit = preds$fit, se.fit = preds$se.fit 
) %>% 
  distinct()

# simulate from multinomial 
set.seed(123)
dd <- rmultinom(n = 1000, size = 1, prob = new_dat$fit) %>% 
  as.matrix() 
apply(dd, 1, sum)

set.seed(123)
dd <- rmultinom(n = 1, size = 1000, prob = new_dat$fit) %>% 
  as.matrix() 
apply(dd, 1, sum)



true_ppn <- c(0.3, 0.2, 0.5)