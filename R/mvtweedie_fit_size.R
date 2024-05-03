## Model Fitting
# Fit data to estimate seasonal changes in size composition 
# Adaptation of rkw_strata_fit to fit a spatially explicit model using mvtweedie
# and mgcv or glmmtmb packages
# Includes sensitivity analyses to evaluate effects of a) slot limits and b)
# excluding fish < 75 cm


library(tidyverse)
library(mvtweedie)
library(mgcv)

# modified prediction function
source(here::here("R", "functions", "pred_mvtweedie2.R"))


size_raw <- readRDS(here::here("data", "rec", "rec_size.rds")) %>% 
  janitor::clean_names() 


dat <- size_raw %>% 
  filter(
    !is.na(size_bin),
    !is.na(lat),
    !is.na(lon),
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
    )
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "km",
    utm_names = c("utm_x", "utm_y")
  ) %>% 
  droplevels() %>% 
  group_by(sample_id) %>% 
  mutate(sample_id_n = length(unique(id))) %>% 
  ungroup()


sample_key <- dat %>% 
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x, 
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
  size_bin = unique(dat$size_bin)
) %>% 
  left_join(., sample_key, by = "sample_id") %>% 
  left_join(
    ., 
    dat %>% 
      group_by(sample_id, size_bin) %>% 
      summarize(
        agg_prob = sum(length(unique(id)))
      ) %>% 
      ungroup(),
    by = c("sample_id", "size_bin")
  ) %>% 
  mutate(
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    agg_ppn = agg_prob / sample_id_n,
    year_n = as.numeric(year),
    year = as.factor(year),
    size_bin = as.factor(size_bin),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
  ) 


# SMU colour palette
size_colour_pal <- c("grey30", "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5", 
                    "#5ab4ac", "#01665e")
names(size_colour_pal) <- c(NA, levels(agg_dat$size_bin))


## DATA FIGURES ----------------------------------------------------------------

# sampling coverage 
size_samp_cov <- ggplot(sample_key) +
  geom_jitter(aes(x = week_n, y = year, size = sample_id_n),
              alpha = 0.4
  ) +
  facet_wrap(~strata) +
  scale_size_continuous(name = "Sample\nSize") +
  scale_x_continuous(
    breaks = c(2, 20, 36, 50),
    labels = c("Jan", "May", "Sep", "Dec")
  ) + 
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )


# stacked bar plot
rec_size_bar <- ggplot(dat) +
  geom_bar(aes(x = month_n, fill = size_bin)) +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = size_colour_pal, name = "Size\nBin") +
  labs(
    y = "Recreational Fishery\nComposition"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 9, 12),
    labels = c("Jan", "May", "Sep", "Dec")
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )

# subset of monthly samples that matches RKW diet
rec_size_bar_summer <- dat %>% 
  filter(month_n %in% c("6", "7", "8", "9", "10")) %>% 
  ggplot(.) +
  geom_bar(aes(x = month_n, fill = size_bin)) +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = size_colour_pal, name = "Size\nBin") +
  labs(
    y = "Recreational Fishery\nComposition"
  ) +
  scale_x_continuous(
    breaks = c(6, 7, 8, 9, 10),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )

png(
  here::here("figs", "ms_figs", "rec_monthly_size_bar.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_size_bar
dev.off()

png(
  here::here("figs", "ms_figs", "rec_monthly_size_bar_summer.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_size_bar_summer
dev.off()


## FIT MODEL -------------------------------------------------------------------

system.time(
  fit <- gam(
    agg_prob ~ 0 + size_bin + s(week_n, by = size_bin, k = 7, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
      s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 25) +
      s(year_n, by = size_bin, k = 4, bs = "tp"),
    data = agg_dat, family = "tw", method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
class(fit) = c( "mvtweedie", class(fit) )
saveRDS(
  fit,
  here::here(
    "data", "model_fits", "mvtweedie", "fit_size_yr_s_mvtw.rds"
  )
)


## DIAGNOSTICS -----------------------------------------------------------------

## CHECK -----------------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_prob == 0) / nrow(agg_dat)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
library(sdmTMB)
fit_sdmTMB <- sdmTMB(
  agg_prob ~ 0 + size_bin + s(week_n, by = size_bin, k = 7, bs = "cc") +
    s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
    s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 25) +
    s(year_n, by = size_bin, k = 4, bs = "tp"),
  data = agg_dat,
  spatial = "off",
  spatiotemporal = "off",
  family = tweedie(link = "log"),
  knots = list(week_n = c(0, 52))
)
saveRDS(
  fit_sdmTMB,
  here::here("data", "model_fits", "mvtweedie", "fit_size_yr_s_sdmTMB.rds")
)

# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
s_sdmTMB <- simulate(fit_sdmTMB, nsim = 500)
sum(s_sdmTMB == 0) / length(s_sdmTMB)

sdmTMBextra::dharma_residuals(s_sdmTMB, fit_sdmTMB)


# look at average stock comp 
avg_sim_comp <- s_sdmTMB %>% 
  as.data.frame() %>% 
  mutate(stock_group = agg_dat$stock_group,
         obs_prob = agg_dat$agg_prob) %>% 
  pivot_longer(
    cols = starts_with("V"),
    names_to = "sim_number",
    values_to = "prob"
  ) %>% 
  group_by(stock_group, sim_number) %>% 
  summarize(
    mean_sim_prob = mean(prob),
    mean_obs_prob = mean(obs_prob)
  )

ggplot(data = avg_sim_comp) +
  geom_boxplot(aes(x = stock_group, y = mean_sim_prob)) +
  geom_point(aes(x = stock_group, y = mean_obs_prob), col = "red") +
  ggsidekick::theme_sleek()
