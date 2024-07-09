## Model Fitting
# Fit equivalent model as mvtweedie_fit but with spatial component estimated
# in sdmTMB
## Severe convergence issues

library(tidyverse)
library(sdmTMB)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() 


dat <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    !fl < 551,
    #exclude samples collected outside areas in relatively close proximity to 
    # SRKW foraging areas
    !strata == "other"
  ) %>% 
  mutate(
    # spatial location variable since some fishing sites have multiple lat/lon
    spatial_loc = paste(lat, lon) %>% 
      as.factor() %>% 
      as.numeric() %>% 
      paste("site", ., sep = ""),
    sample_id = paste(spatial_loc, week_n, year, sep = "_"),
    strata = factor(
      strata,
      levels = c("swiftsure", "swiftsure_nearshore", "renfrew", "vic",
                 "haro", "saanich"),
      labels = c("Swiftsure", "Nitinat", "Renfrew", "Sooke/\nVictoria",
                 "S. Gulf\nIslands", "Saanich")
    ),
    stock_group = fct_relevel(stock_group, "PSD", after = 3)
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "km",
    utm_names = c("utm_x", "utm_y")
  ) %>% 
  droplevels() %>% 
  group_by(sample_id) %>% 
  mutate(sample_id_n = sum(prob)) %>% 
  ungroup()

sample_key <- dat %>% 
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x, shore_dist,
         slot_limit) %>% 
  distinct()


# identify a single spatial location for each strata based on mean
loc_key <- readRDS(
  here::here("data", "spatial", "strata_key.rds")
) %>% 
  # match scale of fitted model
  mutate(
    utm_x = utm_x_m / 1000,
    utm_y = utm_y_m / 1000
  )


# add zero observations
agg_dat <- expand.grid(
  sample_id = unique(dat$sample_id),
  stock_group = unique(dat$stock_group)
) %>% 
  left_join(., sample_key, by = "sample_id", relationship = "many-to-many") %>% 
  left_join(
    ., 
    dat %>% 
      group_by(sample_id, stock_group) %>% 
      summarize(
        agg_prob = sum(prob)
      ) %>% 
      ungroup(),
    by = c("sample_id", "stock_group")
  ) %>% 
  mutate(
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    agg_ppn = agg_prob / sample_id_n,
    year_n = as.numeric(year),
    year = as.factor(year),
    stock_group = as.factor(stock_group),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000,
    sg_year = paste(stock_group, year, sep = "_") %>% 
      as.factor(),
    week_z = scale(week_n) %>% as.numeric(),
    shore_dist_z = scale(shore_dist) %>% as.numeric()
  ) 


## FIT MODEL -------------------------------------------------------------------

mesh_in <- make_mesh(
  agg_dat,
  c("utm_x", "utm_y"),
  n_knots = 100
)

fit_sdmTMB <- sdmTMB(
  agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 7, bs = "cc") #+
    # (1 | sg_year)
  ,
  data = agg_dat,
  groups = "stock_group",
  mesh = mesh_in,
  spatial = "on",
  spatiotemporal = "off",
  family = tweedie(link = "log"),
  knots = list(week_n = c(0, 52)),
  anisotropy = FALSE,
  silent = FALSE
)

