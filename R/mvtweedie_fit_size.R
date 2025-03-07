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
    # remove smallest sizes that are absent from diets entirely
    !fl < 551,
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
      labels = c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew",
                 "Sooke/\nVictoria", "S. Gulf\nIslands", "Saanich")
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
    utm_y_m = utm_y * 1000,
    sg_year = paste(size_bin, year, sep = "_") %>% 
      as.factor()
  ) 
# saveRDS(
#   agg_dat, here::here("data", "rec", "cleaned_ppn_data_rec_size_xy.rds")
# )


# size colour palette
size_colour_pal <- c("grey30", "#8c510a", #"#d8b365", 
                     "#f6e8c3", "#c7eae5", 
                    # "#5ab4ac",
                    "#01665e"
                    )
names(size_colour_pal) <- c(NA, levels(agg_dat$size_bin))

# sample sizes
full_samp_size <- dat %>% 
  group_by(strata) %>% 
  summarize(
    n = length(unique(id))
  )
summer_samp_size <- dat %>% 
  filter(month_n %in% c("5", "6", "7", "8", "9", "10")) %>% 
  group_by(strata) %>% 
  summarize(
    n = length(unique(id))
  )


## DATA FIGURES ----------------------------------------------------------------

# stacked bar plot
rec_size_bar <- ggplot(dat) +
  geom_bar(aes(x = month_n, fill = size_bin)) +
  facet_wrap(~strata, scales = "free_y") +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = size_colour_pal, name = "Size\nBin") +
  labs(
    y = "Recreational Fishery Composition\n(Individual Samples)"
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 9, 12),
    labels = c("Jan", "May", "Sep", "Dec")
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  geom_text(
    data = full_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )


# subset of monthly samples that matches RKW diet
rec_size_bar_summer <- dat %>% 
  filter(month_n %in% c("5", "6", "7", "8", "9", "10")) %>% 
  ggplot(.) +
  geom_bar(aes(x = month_n, fill = size_bin)) +
  facet_wrap(~strata, scales = "free_y") +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = size_colour_pal, name = "Size\nBin") +
  labs(
    y = "Recreational Fishery Composition\n(Individual Samples)"
  ) +
  scale_x_continuous(
    breaks = c(5, 6, 7, 8, 9, 10),
    labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  geom_text(
    data = summer_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )


png(
  here::here("figs", "size_comp_fishery", "rec_monthly_size_bar.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_size_bar
dev.off()

png(
  here::here("figs", "size_comp_fishery", "rec_monthly_size_bar_summer.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_size_bar_summer
dev.off()


## FIT MODEL -------------------------------------------------------------------

system.time(
  fit <- gam(
    agg_prob ~ 0 + size_bin + s(week_n, by = size_bin, k = 20, bs = "cc") +
      # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
      s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 35) +
      s(sg_year, bs = "re"),
    data = agg_dat, family = nb(), method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
saveRDS(
  fit,
  here::here(
    "data", "model_fits", "fit_size_mvtw.rds"
  )
)
fit <- readRDS(
  here::here(
    "data", "model_fits", "fit_size_mvtw.rds"
  )
)


## DIAGNOSTICS -----------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_prob == 0) / nrow(agg_dat)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
library(sdmTMB)
fit_sdmTMB <- sdmTMB(
  agg_prob ~ 0 + size_bin + s(week_n, by = size_bin, k = 10, bs = "tp") +
    # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
    s(utm_y, utm_x, by = size_bin, m = c(0.5, 1), bs = "ds", k = 25)  +
    (1 | sg_year)
    ,
  data = agg_dat,
  spatial = "off",
  spatiotemporal = "off",
  family = sdmTMB::nbinom2(),
  # family = tweedie(link = "log"),
  knots = list(week_n = c(0, 52))
)
saveRDS(
  fit_sdmTMB,
  here::here("data", "model_fits", "fit_size_sdmTMB.rds")
)
fit_sdmTMB <- readRDS(
  here::here("data", "model_fits", "fit_size_sdmTMB.rds")
)


# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
s_sdmTMB <- simulate(fit_sdmTMB, nsim = 500)
sum(s_sdmTMB == 0) / length(s_sdmTMB)


samp <- sdmTMBextra::predict_mle_mcmc(fit_sdmTMB, mcmc_iter = 101, 
                                      mcmc_warmup = 100)
mcmc_res <- residuals(fit_sdmTMB, type = "mle-mcmc", mcmc_samples = samp)

png(
  here::here("figs", "size_comp_fishery", "qq_plot_size.png"),
  height = 4, width = 4, units = "in", res = 250
)
qqnorm(mcmc_res); qqline(mcmc_res)
dev.off()


# check outlier
# qq_dat <- qqnorm(mcmc_res, plot.it = FALSE)
# qq_df <- data.frame(theoretical = qq_dat$x, residual = qq_dat$y) %>% 
#   mutate(
#     dev = abs(residual - qq_dat$x)
#   ) 
# 
# agg_dat[order(-qq_df$dev)[1], ]
# nothing overly unusual about the sample; large ppn of over 85 cm fish, but


# look at average stock comp 
sim_comp <- s_sdmTMB %>% 
  as.data.frame() %>% 
  cbind(
    ., 
    agg_dat %>% 
      select(sample_id, size_bin, strata, week_n)
  ) %>% 
  pivot_longer(
    cols = starts_with("V"),
    names_to = "sim_number",
    values_to = "conc"
  )

# focus on summer samples
sum_samps <- dat %>% 
  filter(month_n %in% c("6", "7", "8", "9", "10")) %>% 
  select(week_n, sample_id) %>% 
  distinct()

week_key <- data.frame(
  week_n = seq(min(sum_samps$week_n), max(sum_samps$week_n), by = 1)
) %>% 
  mutate(
    month = cut(
      week_n, breaks = 5, labels = c("Jun", "Jul", "Aug", "Sep", "Oct"))
  )

avg_sim_comp <- sim_comp %>% 
  group_by(sim_number, sample_id) %>% 
  mutate(
    # calculate number of fish simulated per simulation and sampling event
    sim_sample_conc = sum(conc),
    sim_ppn = ifelse(sim_sample_conc == "0", 0, conc / sim_sample_conc)
  ) %>% 
  filter(week_n %in% sum_samps$week_n) %>% 
  left_join(., week_key, by = "week_n") %>% 
  group_by(strata, month, size_bin, sim_number) %>% 
  summarize(
    mean_sim_ppn = mean(sim_ppn)
  )
avg_obs_comp1 <- agg_dat %>% 
  filter(sample_id %in% sum_samps$sample_id) %>%
  left_join(., week_key, by = "week_n") %>% 
  group_by(strata, month) %>% 
  mutate(
    strata_n = sum(agg_prob)
  ) %>%
  ungroup()
avg_obs_comp <- avg_obs_comp1 %>% 
  group_by(strata, month, strata_n, size_bin) %>% 
  summarize(
    mean_obs_ppn = mean(agg_ppn),
    se_obs_ppn = sqrt(mean_obs_ppn * (1 - mean_obs_ppn) / strata_n),
    lo = pmax(0, mean_obs_ppn - (1.96 * se_obs_ppn)),
    up = pmin(1, mean_obs_ppn + (1.96 * se_obs_ppn))
  ) %>% 
  distinct()

post_sim <- ggplot() +
  geom_boxplot(data = avg_sim_comp ,
               aes(x = month, y = mean_sim_ppn)) +
  geom_pointrange(
    data = avg_obs_comp,
    aes(x = month, y = mean_obs_ppn, ymin = lo, ymax = up), 
    col = "red", alpha = 0.6) +
  ggsidekick::theme_sleek() +
  facet_grid(strata~size_bin) +
  labs(y = "Mean Composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  geom_text(
    data = summer_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

png(
  here::here("figs", "size_comp_fishery", "posterior_sims_size.png"),
  height = 8.5, width = 7.5, units = "in", res = 250
)
post_sim
dev.off()



## PREDICT ---------------------------------------------------------------------

# modified prediction function
source(here::here("R", "functions", "pred_mvtweedie2.R"))


newdata1 <- expand.grid(
  strata = unique(dat$strata),
  week_n = unique(agg_dat$week_n),
  size_bin = levels(agg_dat$size_bin),
  year_n = unique(agg_dat$year_n)
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  mutate(
    # year = as.factor(year_n),
    year = as.factor(year_n),
    sg_year = paste(size_bin, year, sep = "_") %>% 
      as.factor(),
    strata = factor(strata, levels = levels(agg_dat$strata))
  ) %>% 
  filter(
    !strata == "Saanich"
  )

newdata <- newdata1 %>% 
  filter(
    year == agg_dat$year[1]
  )


# year-specific predictions (average strata)
pred3 = pred_dummy(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = newdata1
)

newdata_yr <- cbind( newdata1, fit=pred3$fit, se.fit=pred3$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) %>% 
  filter(
    !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & week_n < 22),
    !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & week_n > 40)
  )

year_preds <- ggplot(newdata_yr, aes(week_n, fit)) +
  geom_line(aes(colour = year)) +
  facet_grid(size_bin~strata, scales = "free_y") +
  scale_colour_discrete(name = "Year") +
  coord_cartesian(xlim = c(20, 43)#, ylim = c(0, 1)
  ) +
  labs(y="Predicted Proportion") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top",
        axis.title.x = element_blank()) +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
  )



## average seasonal predictions by size bin
newdata_season <- newdata %>% 
  filter(strata == agg_dat$strata[1])
excl2 <- grepl("utm", gratia::smooths(fit)) | 
  grepl("year", gratia::smooths(fit))
utm_yr_coefs <- gratia::smooths(fit)[excl2]
pred_season = pred_dummy(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = newdata_season,
  exclude = utm_yr_coefs
)

newdata_season2 <- cbind( 
  newdata_season, fit=pred_season$fit, se.fit=pred_season$se.fit 
  ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) 

season_preds <- ggplot(newdata_season2, aes(week_n, fit)) +
  geom_line(aes(colour = size_bin)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = size_bin), alpha = 0.5) +
  facet_wrap(~size_bin, scales = "free_y") +
  scale_fill_manual(name = "Size Bin", values = size_colour_pal) +
  scale_colour_manual(name = "Size Bin", values = size_colour_pal) +
  coord_cartesian(xlim = c(20, 43)#, ylim = c(0, 1)
                  ) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top")  +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
  )

season_preds_fullx <- season_preds +
  coord_cartesian(xlim = c(0, 52), expand = FALSE)  +
  scale_x_continuous()


# average across year predictions
excl <- grepl("year", gratia::smooths(fit))
yr_coefs <- gratia::smooths(fit)[excl]
full_pred = pred_dummy(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = newdata,
  exclude = yr_coefs
)

newdata_full <- cbind( newdata, fit=full_pred$fit, se.fit=full_pred$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) 

# integrates out yearly variation
summer_preds <- ggplot(newdata_full, aes(week_n, fit)) +
  geom_point(
    data = agg_dat %>% 
      filter(week_n %in% newdata$week_n,
             !strata == "Saanich"),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "red") +
  facet_grid(size_bin~strata) +
  coord_cartesian(xlim = c(20, 43), ylim = c(0, 1)) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
  )

summer_preds_fullx <- summer_preds +
  coord_cartesian(xlim = c(0, 52)) +
  scale_x_continuous()


## stacked ribbon predictions
summer_pred_stacked <- ggplot(
  data = newdata_full %>% 
    filter(
      week_n > 19 & week_n < 44,
      !(strata %in% c("Swiftsure", "Nitinat", "Renfrew") & week_n < 22),
      !(strata %in% c("Swiftsure", "Nitinat", "Renfrew") & week_n > 40)
    ), 
  aes(x = week_n)
) +
  geom_area(aes(y = fit, colour = size_bin, fill = size_bin), 
            stat = "identity") +
  scale_fill_manual(values = size_colour_pal) +
  scale_colour_manual(values = size_colour_pal) +
  labs(y = "Predicted Mean Composition of Fishery Sample", x = "Week") +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "off",
    axis.text = element_text(size=9),
    plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points"),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA)) +
  facet_wrap(~strata)  +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
  )
summer_pred_legend <- cowplot::get_legend(
  ggplot(
    data = newdata_full, aes(x = week_n)
  ) +
    geom_area(aes(y = fit, colour = size_bin, fill = size_bin), 
              stat = "identity") +
    scale_fill_manual(name = "Size Bin", values = size_colour_pal) +
    scale_colour_manual(name = "Size Bin", values = size_colour_pal) +
    theme(
      legend.title = element_text(size = 15),
      legend.key.size = unit(1.4, "lines"),  # Adjust the size of the legend keys
      legend.text = element_text(size = 14)   # Adjust the text size of the legend
    )
)


png(
  here::here("figs", "size_comp_fishery", "size_smooth_preds_chinook_year.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
year_preds
dev.off()

png(
  here::here("figs", "size_comp_fishery", "size_season_preds_chinook.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
season_preds
dev.off()

png(
  here::here("figs", "size_comp_fishery", "size_smooth_preds_chinook.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
summer_preds
dev.off()

png(
  here::here("figs", "size_comp_fishery", "size_smooth_preds_chinook_xaxis.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
summer_preds_fullx
dev.off()

png(
  here::here("figs", "size_comp_fishery", "size_smooth_preds_chinook_stacked.png"),
  height = 6.5, width = 8.5, units = "in", res = 250
)
cowplot::ggdraw(summer_pred_stacked) +
  cowplot::draw_plot(summer_pred_legend,
                     height = 0.33, x = 0.33, y = 0.08)
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
  week_n = c(20, 25, 29, 33, 36, 40),
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
    utm_x = X / 1000,
    year_n = unique(agg_dat$year_n)[1],
    sg_year = paste(size_bin, year_n, sep = "_") %>% 
      as.factor()
  ) 


# exclude seasonal and annual effects
excl3 <- #grepl("week_n", gratia::smooths(fit)) | 
  grepl("year", gratia::smooths(fit))
yr_coefs <- gratia::smooths(fit)[excl3]
pred_sp <- pred_dummy(
  fit,
  se.fit = TRUE,
  category_name = "size_bin",
  origdata = agg_dat,
  newdata = new_dat_sp,
  exclude = yr_coefs
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


# calculate the average summer distribution (SE dropped to avoid delta method):
# 1) calculate mean comp for each stock/cell
# 2) rescale each cell so that summed comp = 1
# 3) calculate scaled_fit for each size bin
new_dat_sp_plot <- cbind(
  new_dat_sp, fit=pred_sp$fit#, se.fit=pred_sp$se.fit 
) %>% 
  group_by(size_bin, X, Y, utm_y, utm_x) %>% 
  summarize(
    mean_fit = mean(fit)
  ) %>% 
  group_by(X, Y, utm_y, utm_x) %>% 
  mutate(fit = mean_fit / sum(mean_fit)) %>% 
  group_by(size_bin) %>% 
  mutate(
    scaled_fit = fit / max(fit)
  ) %>% 
  ungroup() 


spatial_pred <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_wrap(
    ~ size_bin
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
    strip.text = element_text(size = 9)
  ) 


spatial_pred_scaled <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = scaled_fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_wrap(
    ~ size_bin
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
    strip.text = element_text(size = 9)
  ) 


# spatial_pred_se <- ggplot() +
#   geom_raster(data = new_dat_sp_plot, 
#               aes(x = X, y = Y, fill = se.fit)) +
#   geom_sf(data = coast, color = "black", fill = "grey") +
#   facet_wrap(
#     ~ size_bin
#   ) +
#   scale_fill_gradient2(
#     name = "Predicted SE\nof Proportion\nEstimate"
#   ) +
#   ggsidekick::theme_sleek()  +
#   theme(
#     axis.title = element_blank(),
#     axis.text = element_blank(), 
#     axis.ticks = element_blank(),
#     legend.position = "top",
#     strip.text = element_text(size = 5)
#   ) 


png(
  here::here("figs", "size_comp_fishery", "size_spatial_preds.png"),
  height = 4, width = 6, units = "in", res = 250
)
spatial_pred
dev.off()

png(
  here::here("figs", "size_comp_fishery", "size_spatial_preds_scaled.png"),
  height = 4, width = 6, units = "in", res = 250
)
spatial_pred_scaled
dev.off()

# png(
#   here::here("figs", "size_comp_fishery", "spatial_preds_se.png"),
#   height = 4, width = 6, units = "in", res = 250
# )
# spatial_pred_se
# dev.off()
# 



