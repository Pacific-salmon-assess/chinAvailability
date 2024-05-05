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
    ),
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 651, 751, 851, Inf), 
      labels = c("<65", "65-75", "75-85", ">85")
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
size_colour_pal <- c("grey30", "#8c510a", "#d8b365", "#f6e8c3", "#c7eae5"#, 
                    # "#5ab4ac", "#01665e"
                    )
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
      s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 25) #+
    #remove since didn't fit well with sdmTMB
#      s(year_n, by = size_bin, k = 4, bs = "tp")
,
    data = agg_dat, family = "tw", method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
class(fit) = c( "mvtweedie", class(fit) )
saveRDS(
  fit,
  here::here(
    "data", "model_fits", "mvtweedie", "fit_size_mvtw.rds"
  )
)
fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_size_mvtw.rds"
  )
)


## DIAGNOSTICS -----------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_prob == 0) / nrow(agg_dat)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
library(sdmTMB)
fit_sdmTMB <- sdmTMB(
  agg_prob ~ 0 + size_bin + s(week_n, by = size_bin, k = 7, bs = "cc") +
    s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
    s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 25)# +
    # (1 | year)
    ,
  data = agg_dat,
  spatial = "off",
  spatiotemporal = "off",
  family = tweedie(link = "log"),
  knots = list(week_n = c(0, 52))
)
saveRDS(
  fit_sdmTMB,
  here::here("data", "model_fits", "mvtweedie", "fit_size_sdmTMB.rds")
)

# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
s_sdmTMB <- simulate(fit_sdmTMB, nsim = 500)
sum(s_sdmTMB == 0) / length(s_sdmTMB)

sdmTMBextra::dharma_residuals(s_sdmTMB, fit_sdmTMB)


# look at average stock comp 
avg_sim_comp <- s_sdmTMB %>% 
  as.data.frame() %>% 
  mutate(size_bin = agg_dat$size_bin,
         obs_prob = agg_dat$agg_prob) %>% 
  pivot_longer(
    cols = starts_with("V"),
    names_to = "sim_number",
    values_to = "prob"
  ) %>% 
  group_by(size_bin, sim_number) %>% 
  summarize(
    mean_sim_prob = mean(prob),
    mean_obs_prob = mean(obs_prob)
  )

ggplot(data = avg_sim_comp) +
  geom_boxplot(aes(x = size_bin, y = mean_sim_prob)) +
  geom_point(aes(x = size_bin, y = mean_obs_prob), col = "red") +
  ggsidekick::theme_sleek()


## PREDICT ---------------------------------------------------------------------

newdata <- expand.grid(
  strata = unique(dat$strata),
  week_n = unique(agg_dat$week_n),
  size_bin = levels(agg_dat$size_bin)
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  mutate(
    strata = factor(strata, levels = levels(agg_dat$strata))
  ) %>% 
  filter(
    !strata == "Saanich"
  )

pred = predict(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = newdata
)

newdata2 <- cbind( newdata, fit=pred$fit, se.fit=pred$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) %>%
  filter(!strata == "Saanich")

summer_preds <- ggplot(newdata2, aes(week_n, fit)) +
  geom_point(
    data = agg_dat %>% 
      filter(week_n %in% newdata2$week_n,
             !strata == "Saanich"),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "red") +
  facet_grid(size_bin~strata) +
  coord_cartesian(xlim = c(24, 40), ylim = c(0, 1)) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top")

summer_preds_fullx <- summer_preds +
  coord_cartesian(xlim = c(0, 52)) 


## stacked ribbon predictions
summer_pred_stacked <- ggplot(
  data = newdata2 %>% 
    filter(week_n > 21 & week_n < 39), 
  aes(x = week_n)
) +
  geom_area(aes(y = fit, colour = size_bin, fill = size_bin), 
            stat = "identity") +
  scale_fill_manual(name = "Size Bin", values = size_colour_pal) +
  scale_colour_manual(name = "Size Bin", values = size_colour_pal) +
  labs(y = "Predicted Mean Composition of Fishery Sample", x = "Week") +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "top",
    axis.text = element_text(size=9),
    plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points"),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA)) +
  facet_wrap(~strata) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37),
    labels = c("Jun", "Jul", "Aug", "Sep")
  )


png(
  here::here("figs", "ms_figs", "size_smooth_preds_chinook.png"),
  height = 8.5, width = 6.5, units = "in", res = 250
)
summer_preds
dev.off()

png(
  here::here("figs", "ms_figs", "size_smooth_preds_chinook_xaxis.png"),
  height = 8.5, width = 6.5, units = "in", res = 250
)
summer_preds_fullx
dev.off()

png(
  here::here("figs", "ms_figs", "size_smooth_preds_chinook_stacked.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
summer_pred_stacked
dev.off()


## spatial predictions


# pfma sf dataframe to subset pred grid
pfma_sf <- readRDS(
  here::here(
    "data", "spatial", "pfma_subareas_sBC.rds"
  )
) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m"))

# predictive grid (1 km and 1 km)
pred_grid <- readRDS(
  here::here(
    "data", "spatial", "pred_bathy_grid_utm_no_bark.RDS"
  )
) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>%
  sf::st_crop(
    ., 
    xmin = min(agg_dat$utm_x_m) - 500, 
    ymin = min(agg_dat$utm_y_m) - 500,
    xmax = max(agg_dat$utm_x_m) + 500, 
    ymax = max(agg_dat$utm_y_m) + 500
  ) %>% 
  # constrain grid to Canadian waters
  sf::st_intersection(., pfma_sf)


new_dat_sp <- expand.grid(
  week_n = c(20, 25, 29, 33, 36),
  size_bin = levels(agg_dat$size_bin)
) %>%
  mutate(
    fac = paste(week_n, as.factor(size_bin), sep = "_")
  ) %>%
  distinct() %>%
  split(., .$fac) %>%
  purrr::map(., function (x) {
    data.frame(
      sf::st_coordinates(pred_grid[ , 1]),
      week_n = x$week_n,
      size_bin = x$size_bin
    )
  }) %>%
  bind_rows %>%
  mutate(
    utm_y = Y / 1000,
    utm_x = X / 1000
  )

# key for month labels
month_key <- data.frame(
  week_n = unique(new_dat_sp$week_n)
) %>% 
  mutate(
    month = c("May", "Jun", "Jul", "Aug", "Sep") %>% 
      as.factor() %>% 
      fct_reorder(., week_n)
  )

pred_sp <- predict(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = new_dat_sp
)

coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                         returnclass = "sf"), 
               rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>%
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>% 
  sf::st_crop(
    ., 
    xmin = min(new_dat_sp$X) - 1500, 
    ymin = min(new_dat_sp$Y) - 2000,
    xmax = max(new_dat_sp$X) + 1500, 
    ymax = max(new_dat_sp$Y) + 2000
  )


new_dat_sp_plot <- cbind(
  new_dat_sp, fit=pred_sp$fit, se.fit=pred_sp$se.fit 
) %>% 
  group_by(size_bin) %>% 
  mutate(
    scaled_fit = fit / max(fit)
  ) %>% 
  ungroup() %>% 
  left_join(
    ., month_key, by = "week_n"
  )

spatial_pred <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_grid(
    size_bin ~ month
  ) +
  scale_fill_viridis_c(
    name = "Predicted Proportion\nof Rec Catch"
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 5)
  ) 


spatial_pred_scaled <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = scaled_fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_grid(
    size_bin ~ month
  ) +
  scale_fill_viridis_c(
    option = "A",
    name = "Predicted Scaled\nProportion\nof Rec Catch"
  ) +
  ggsidekick::theme_sleek()  +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 5)
  ) 


spatial_pred_se <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = se.fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_grid(
    size_bin ~ month
  ) +
  scale_fill_gradient2(
    name = "Predicted SE\nof Proportion\nEstimate"
  ) +
  ggsidekick::theme_sleek()  +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "top",
    strip.text = element_text(size = 5)
  ) 


png(
  here::here("figs", "ms_figs", "spatial_size_preds.png"),
  height = 6, width = 6, units = "in", res = 250
)
spatial_pred
dev.off()

png(
  here::here("figs", "ms_figs", "spatial_preds_size_scaled.png"),
  height = 6, width = 6, units = "in", res = 250
)
spatial_pred_scaled
dev.off()

png(
  here::here("figs", "ms_figs", "spatial_preds_size_se.png"),
  height = 6, width = 6, units = "in", res = 250
)
spatial_pred_se
dev.off()