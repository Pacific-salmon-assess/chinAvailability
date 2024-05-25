## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Adaptation of rkw_strata_fit to fit a spatially explicit model using mvtweedie
# and mgcv or glmmtmb packages
# Includes sensitivity analyses to evaluate effects of a) slot limits and b)
# excluding fish < 75 cm


library(tidyverse)
library(mvtweedie)
library(mgcv)

# modified prediction function
source(here::here("R", "functions", "pred_mvtweedie2.R"))


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


# trim data to exclude swiftsure samples west of Nitinat
# NOTE: doesn't strongly impact preds
# w_border <- dat %>% 
#   filter(strata == "Nitinat") %>% 
#   pull(utm_x) %>% 
#   min()
# dat <- dat %>% 
#   filter(!utm_x < w_border - 1)


# add zero observations
agg_dat <- expand.grid(
  sample_id = unique(dat$sample_id),
  stock_group = unique(dat$stock_group)
) %>% 
  left_join(., sample_key, by = "sample_id") %>% 
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
      as.factor()
  ) 
# saveRDS(
#   agg_dat, here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
# )


# SMU colour palette
smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(dat$stock_group)

# hatchery origin colour palette
hatchery_colour_pal <- c("#31a354", "#e5f5e0", "grey30", "grey60", "#756bb1",
                         "#efedf5")
names(hatchery_colour_pal) <- levels(dat$origin2)


## DATA FIGURES ----------------------------------------------------------------

# sampling coverage 
rec_samp_cov <- ggplot(sample_key) +
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
rec_samp_bar <- ggplot(dat) +
  geom_bar(aes(x = month_n, y = prob, fill = stock_group), 
           stat = "identity") +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = "Stock\nGroup") +
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
rec_samp_bar_summer <- dat %>% 
  filter(month_n %in% c("6", "7", "8", "9", "10")) %>% 
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = stock_group), 
           stat = "identity") +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = "Stock\nGroup") +
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


# as above but for hatchery origin
rec_samp_bar_h <- ggplot(dat) +
  geom_bar(aes(x = month_n, y = prob, fill = origin2),
           stat = "identity") +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = hatchery_colour_pal, name = "Hatchery\nOrigin") +
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
rec_samp_bar_summer_h <- dat %>%
  filter(month_n %in% c("6", "7", "8", "9", "10"), year > 2017) %>%
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = origin2),
           stat = "identity") +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = hatchery_colour_pal, name = "Hatchery\nOrigin") +
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
  here::here("figs", "stock_comp_fishery", "rec_temporal_sample_coverage.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_cov
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "rec_monthly_comp_bar.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "rec_monthly_comp_bar_summer.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_summer
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "rec_monthly_hatchery_comp_bar.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_h
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "rec_monthly_hatchery_comp_bar_summer.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_summer_h
dev.off()


## FIT MODEL -------------------------------------------------------------------


# Does not include year effects
# system.time(
#   fit <- gam(
#     agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 7, bs = "cc") +
#       s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25) ,
#     data = agg_dat, family = "tw", method = "REML",
#     knots = list(week_n = c(0, 52))
#   )
# )
# # ~850 seconds to converge
# saveRDS(
#   fit,
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_tw.rds")
# )
# 
# class(fit) = c( "mvtweedie", class(fit) )
# saveRDS(
#   fit,
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_mvtw.rds")
# )


# Includes year/stock as RE; remove global smooth after sdmTMB v fails to 
# converge
system.time(
  fit2 <- gam(
    agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 7, bs = "cc") +
      # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25) +
      s(sg_year, bs = "re"),
    data = agg_dat, family = "tw", method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
class(fit2) = c( "mvtweedie", class(fit2) )
saveRDS(
  fit2,
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_mvtw.rds")
)


# Includes smooth for year by stock group
# system.time(
#   fit3 <- gam(
#     agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 7, bs = "cc") +
#       s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25) +
#       s(year_n, by = stock_group, k = 4, bs = "tp"),
#     data = agg_dat, family = "tw", method = "REML",
#     knots = list(week_n = c(0, 52))
#   )
# )
# class(fit3) = c( "mvtweedie", class(fit3) )
# saveRDS(
#   fit3,
#   here::here(
#     "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
#   )
# )



## favor third model since it generates year-specific estimates and appears
# to converge well; UPDATE -- does not converge with sdmTMB switch to RIs
# fit_raw <- readRDS(
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_tw.rds"))
# fit <- readRDS(
#     here::here(
#       "data", "model_fits", "mvtweedie", "fit_spatial_fishery_mvtw.rds")
#     )
fit2 <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_mvtw.rds")
)
# fit3 <- readRDS(
#   here::here(
#       "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
#     )
# )
  

## CHECK -----------------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_prob == 0) / nrow(agg_dat)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
# library(sdmTMB)
# fit_sdmTMB <- sdmTMB(
# agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 7, bs = "cc") +
#   # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 15) +
#   s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 15) +
#   # s(year_n, by = stock_group, bs = "tp", k = 7)
#   # (1 | year)
#   (1 | sg_year)
# ,
#   data = agg_dat,
#   spatial = "off",
#   spatiotemporal = "off",
#   family = tweedie(link = "log"),
#   knots = list(week_n = c(0, 52))
# )
# saveRDS(
#   fit_sdmTMB,
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_sdmTMB.rds")
# )
fit_sdmTMB <- readRDS(
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_sdmTMB.rds")
)



# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
s_sdmTMB <- simulate(fit_sdmTMB, nsim = 500)
sum(s_sdmTMB == 0) / length(s_sdmTMB)

png(
  here::here("figs", "stock_comp_fishery", "qq_plot.png"),
  height = 4, width = 4, units = "in", res = 250
)
sdmTMBextra::dharma_residuals(s_sdmTMB, fit_sdmTMB)
dev.off()


# look at average stock comp 
sim_comp <- s_sdmTMB %>% 
  as.data.frame() %>% 
  cbind(
    ., 
    agg_dat %>% 
      select(sample_id, stock_group, strata, week_n      )
    ) %>% 
  # mutate(stock_group = agg_dat$stock_group,
  #        obs_prob = agg_dat$agg_prob) %>% 
  pivot_longer(
    cols = starts_with("V"),
    names_to = "sim_number",
    values_to = "conc"
  )

avg_sim_comp <- sim_comp %>% 
  group_by(sim_number, sample_id) %>% 
  mutate(
    # calculate number of fish simulated per simulation and sampling event
    sim_sample_conc = sum(conc),
    sim_ppn = conc / sim_sample_conc
  ) %>% 
  filter(!is.na(sim_ppn),
         week_n > 24 & week_n < 40) %>% 
  group_by(strata, stock_group, sim_number) %>% 
  summarize(
    mean_sim_ppn = mean(sim_ppn)
  )
avg_obs_comp <- agg_dat %>% 
  filter(week_n > 24 & week_n < 40) %>% 
  group_by(strata, stock_group) %>% 
  summarize(
    mean_obs_ppn = mean(agg_ppn)
  )

post_sim <- ggplot() +
  geom_boxplot(data = avg_sim_comp ,
               aes(x = stock_group, y = mean_sim_ppn)) +
  geom_point(data = avg_obs_comp,
             aes(x = stock_group, y = mean_obs_ppn), col = "red") +
  ggsidekick::theme_sleek() +
  facet_wrap(~strata, ncol = 2) +
  labs(y = "Mean Stock Composition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank())


png(
  here::here("figs", "stock_comp_fishery", "posterior_sims.png"),
  height = 7.5, width = 4.5, units = "in", res = 250
)
post_sim
dev.off()



## PREDICT ---------------------------------------------------------------------

newdata <- expand.grid(
  strata = unique(dat$strata),
  week_n = unique(agg_dat$week_n),
  stock_group = levels(agg_dat$stock_group),
  year_n = unique(agg_dat$year_n)
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  mutate(
    year = as.factor(year_n),
    sg_year = paste(stock_group, year, sep = "_") %>% 
      as.factor(),
    strata = factor(strata, levels = levels(agg_dat$strata))
  ) %>% 
  filter(
    !strata == "Saanich"
  )


# fit 1 
# pred = predict(
#   fit,
#   se.fit = TRUE,
#   category_name = "stock_group",
#   origdata = agg_dat,
#   newdata = newdata
# )

# pred2 = predict(
#   fit2,
#   # se.fit = TRUE,
#   category_name = "stock_group",
#   origdata = agg_dat,
#   newdata = newdata
# )

# year-specific predictions (average strata)
pred3 = pred_dummy(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata
)

newdata_yr <- cbind( newdata, fit=pred3$fit, se.fit=pred3$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) 

year_preds <- ggplot(newdata_yr, aes(week_n, fit)) +
  geom_line(aes(colour = year)) +
  facet_grid(strata~stock_group, scales = "free_y") +
  scale_colour_discrete(name = "Stock Group") +
  scale_colour_discrete(name = "Stock Group") +
  coord_cartesian(xlim = c(24, 40)#, ylim = c(0, 1)
  ) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = c(25, 29.25, 33.5, 38),
    labels = c("Jun", "Jul", "Aug", "Sep")
  )



# predictions used for "average" effects integrating out year
newdata_b <- newdata %>% 
  filter(year == agg_dat$year[1])

# average seasonal proportion by stock
newdata_a <- newdata_b %>% 
  filter(strata == agg_dat$strata[1])
excl2 <- grepl("utm", gratia::smooths(fit2)) | 
  grepl("year", gratia::smooths(fit2))
utm_yr_coefs <- gratia::smooths(fit2)[excl2]
pred3a = pred_dummy(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata_a,
  exclude = utm_yr_coefs
)

newdata3a <- cbind( newdata_a, fit=pred3a$fit, se.fit=pred3a$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) 

season_preds <- ggplot(newdata3a, aes(week_n, fit)) +
  geom_line(aes(colour = stock_group)) +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = stock_group), alpha = 0.5) +
  facet_wrap(~stock_group, scales = "free_y") +
  scale_fill_manual(name = "Stock Group", values = smu_colour_pal) +
  scale_colour_manual(name = "Stock Group", values = smu_colour_pal) +
  coord_cartesian(xlim = c(24, 40)#, ylim = c(0, 1)
                  ) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = c(25, 29.25, 33.5, 38),
    labels = c("Jun", "Jul", "Aug", "Sep")
  )

season_preds_fullx <- season_preds +
  coord_cartesian(xlim = c(0, 52), expand = FALSE)  +
  scale_x_continuous()


# average across year predictions
excl <- grepl("year", gratia::smooths(fit2))
yr_coefs <- gratia::smooths(fit2)[excl]
pred3b = pred_dummy(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata_b,
  exclude = yr_coefs
)

newdata3b <- cbind( newdata_b, fit=pred3b$fit, se.fit=pred3b$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) 

# integrates out yearly variation
summer_preds <- ggplot(newdata3b, aes(week_n, fit)) +
  geom_point(
    data = agg_dat %>% 
      filter(week_n %in% newdata$week_n,
             !strata == "Saanich"),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "red") +
  facet_grid(stock_group~strata) +
  coord_cartesian(xlim = c(24, 40), ylim = c(0, 1)) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = c(25, 29.25, 33.5, 38),
    labels = c("Jun", "Jul", "Aug", "Sep")
  )

summer_preds_fullx <- summer_preds +
  coord_cartesian(xlim = c(0, 52)) +
  scale_x_continuous()


## stacked ribbon predictions
summer_pred_stacked <- ggplot(
  data = newdata3b %>% 
    filter(week_n > 21 & week_n < 39), 
  aes(x = week_n)
) +
  geom_area(aes(y = fit, colour = stock_group, fill = stock_group), 
            stat = "identity") +
  scale_fill_manual(name = "Stock Group", values = smu_colour_pal) +
  scale_colour_manual(name = "Stock Group", values = smu_colour_pal) +
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
  here::here("figs", "stock_comp_fishery", "smooth_preds_chinook_year.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
year_preds
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "season_preds_chinook.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
season_preds
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "season_preds_chinook_xaxis.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
season_preds_fullx
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "smooth_preds_chinook.png"),
  height = 8.5, width = 6.5, units = "in", res = 250
)
summer_preds
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "smooth_preds_chinook_xaxis.png"),
  height = 8.5, width = 6.5, units = "in", res = 250
)
summer_preds_fullx
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "smooth_preds_chinook_stacked.png"),
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
  stock_group = levels(agg_dat$stock_group)
) %>%
  mutate(
    fac = paste(week_n, as.factor(stock_group), sep = "_")
  ) %>%
  distinct() %>%
  split(., .$fac) %>%
  purrr::map(., function (x) {
    data.frame(
      sf::st_coordinates(pred_grid[ , 1]),
      week_n = x$week_n,
      stock_group = x$stock_group
    )
  }) %>%
  bind_rows %>%
  mutate(
    utm_y = Y / 1000,
    utm_x = X / 1000,
    year_n = unique(agg_dat$year_n)[1],
    sg_year = paste(stock_group, year_n, sep = "_") %>% 
      as.factor()
  ) %>% 
  filter(
    week_n == "29"
  )
  
# key for month labels
# month_key <- data.frame(
#   week_n = unique(new_dat_sp$week_n)
# ) %>% 
#   mutate(
#     month = c("May", "Jun", "Jul", "Aug", "Sep") %>% 
#       as.factor() %>% 
#       fct_reorder(., week_n)
#   )

excl3 <- grepl("week_n", gratia::smooths(fit2)) | 
  grepl("year", gratia::smooths(fit2))
week_yr_coefs <- gratia::smooths(fit2)[excl3]
pred_sp <- pred_dummy(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = new_dat_sp,
  exclude = week_yr_coefs
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
  group_by(stock_group) %>% 
  mutate(
    scaled_fit = fit / max(fit)
  ) %>% 
  ungroup() %>% 
  filter(week_n == "29")#%>% 
  # left_join(
  #   ., month_key, by = "week_n"
  # )

spatial_pred <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_wrap(
    ~ stock_group
  ) +
  # facet_grid(
  #   stock_group ~ month
  # ) +
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
  facet_wrap(
    ~ stock_group
  ) +
  # facet_grid(
  #   stock_group ~ month
  # ) +
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
  facet_wrap(
    ~ stock_group
  ) +
  # facet_grid(
  #   stock_group ~ month
  # ) +
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
  here::here("figs", "stock_comp_fishery", "spatial_preds.png"),
  height = 4, width = 6, units = "in", res = 250
)
spatial_pred
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "spatial_preds_scaled.png"),
  height = 4, width = 6, units = "in", res = 250
)
spatial_pred_scaled
dev.off()

png(
  here::here("figs", "stock_comp_fishery", "spatial_preds_se.png"),
  height = 4, width = 6, units = "in", res = 250
)
spatial_pred_se
dev.off()


## SENSITIVITY ANALYSES --------------------------------------------------------

# 1) Evaluate impact of slot limit by fitting model only to data west of Sooke
# and pre-2019
# 2) Evaluate impact of size preference of SRKW by fitting model only to samples
# greater than 750 mm fl


## Slot limit analysis
agg_dat_slot <- expand.grid(
  sample_id = unique(dat$sample_id),
  stock_group = unique(dat$stock_group)#,
  # slot_limit = unique(dat$slot_limit)
) %>% 
  left_join(., sample_key, by = "sample_id") %>% 
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
    sg_year = paste(stock_group, year) %>% as.factor(),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
  ) 

system.time(
  fit_slot <- gam(
    agg_prob ~ 0 + stock_group*slot_limit + 
      s(week_n, by = stock_group, k = 7, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25) +
      s(sg_year, bs = "re"),
      # s(year_n, by = stock_group, k = 4, bs = "tp"),
    data = agg_dat_slot, family = "tw", method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
class(fit_slot) = c( "mvtweedie", class(fit_slot) )
saveRDS(
  fit_slot,
  here::here(
    "data", "model_fits", "mvtweedie", "fit_slot.rds"
  )
)


## estimate slot limit period effects
slot_pars <- broom::tidy(fit_slot, parametric = TRUE, conf.int = TRUE) %>% 
  filter(grepl("slot", term)) %>% 
  mutate(
    stock_group = levels(agg_dat$stock_group) %>% 
      factor(., levels = levels(agg_dat$stock_group))
  ) 

slot_plot <- ggplot(slot_pars) +
  geom_pointrange(
    aes(x = stock_group, y = estimate,  ymin = conf.low, ymax = conf.high, 
        fill = stock_group), 
    shape = 21) +
  scale_fill_manual(values = smu_colour_pal) +
  geom_hline(aes(yintercept = 0), lty = 2) +
  ggsidekick::theme_sleek() +
  labs(y = "Effect of Slot Limit") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


png(
  here::here("figs", "stock_comp_fishery", "slot_limit_effect.png"),
  height = 3.5, width = 5, units = "in", res = 250
)
slot_plot
dev.off()



newdata_slot <- expand.grid(
  strata = levels(agg_dat_slot$strata),
  week_n = seq(25, 38, by = 0.25),
  stock_group = levels(agg_dat_slot$stock_group),
  slot_limit = unique(agg_dat_slot$slot_limit),
  year_n = unique(agg_dat_slot$year)[1]
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  mutate(
    strata = factor(strata, levels = levels(agg_dat_slot$strata)),
    sg_year = paste(stock_group, year_n, sep = "_") %>% as.factor()
  ) 

excl <- grepl("year", gratia::smooths(fit_slot))
yr_coefs <- gratia::smooths(fit_slot)[excl]
pred_slot = pred_dummy(
  fit_slot,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata_slot,
  exclude = yr_coefs
)
data_slot = cbind(newdata_slot, fit=pred_slot$fit, se.fit=pred_slot$se.fit )
data_slot$lower = data_slot$fit + (qnorm(0.025)*data_slot$se.fit)
data_slot$upper = data_slot$fit + (qnorm(0.975)*data_slot$se.fit)


# focus on strata that introduced slot limits during sampling period
data_slot2 <- data_slot %>% 
  filter(
    strata %in% c("Swiftsure", "Nitinat", "Renfrew")
  )

slot_pred_smooth <- ggplot(
  data_slot2,
  aes(week_n, fit, colour = slot_limit, fill = slot_limit)
) +
  geom_point(
    data = agg_dat_slot %>%
      filter(strata %in% data_slot2$strata),
    aes(x = week_n, y = agg_ppn, size = sample_id_n, colour = slot_limit),
    alpha = 0.3, position = position_dodge(width = 0.5) 
  ) +
  geom_line() +
  # geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  facet_grid(stock_group ~ strata) +
  coord_cartesian(ylim = c(0,1), xlim = c(25, 38)) +
  labs(
    y = "Predicted Proportion of Fishery Sample",
    fill = "Slot\nLimit",
    colour = "Slot\nLimit",
    size = "Sample\nSize"
  ) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) 

png(
  here::here("figs", "stock_comp_fishery", "smooth_preds_slot_limit.png"),
  height = 5, width = 5, units = "in", res = 250
)
slot_pred_smooth
dev.off()


## Slot limit analysis
large_dat <- dat %>%
  filter(fl >= 750) %>%
  group_by(sample_id) %>%
  mutate(sample_id_n = sum(prob)) %>%
  ungroup()
sample_key_large <- large_dat %>%
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x) %>%
  distinct()

agg_dat_large <- expand.grid(
  sample_id = unique(large_dat$sample_id),
  stock_group = unique(large_dat$stock_group)
) %>% 
  left_join(., sample_key_large, by = "sample_id") %>% 
  left_join(
    ., 
    large_dat %>% 
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
    sg_year = paste(stock_group, year) %>% as.factor()
  ) 

# Includes smooth for year by stock group
system.time(
  fit_large <- gam(
    agg_prob ~ 0 + stock_group + 
      s(week_n, by = stock_group, k = 7, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 25) +
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25) +
      s(sg_year, bs = "re"),
    data = agg_dat_large, family = "tw", method = "REML",
    knots = list(week_n = c(0, 52))
  )
)
class(fit_large) = c( "mvtweedie", class(fit_large) )
saveRDS(
  fit_large,
  here::here(
    "data", "model_fits", "mvtweedie", "fit_large.rds"
  )
)


excl <- grepl("year", gratia::smooths(fit_large))
yr_coefs <- gratia::smooths(fit_large)[excl]
pred_large = pred_dummy(
  fit_large,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat_large,
  newdata = newdata_b,
  exclude = yr_coefs
)


newdata_large <- cbind( newdata_b, fit=pred_large$fit, se.fit=pred_large$se.fit ) %>%
  mutate(
    lower = fit + (qnorm(0.025)*se.fit),
    upper = fit + (qnorm(0.975)*se.fit)
  ) %>%
  filter(!strata == "Saanich")


large_pred_smooth <- ggplot(newdata_large, aes(week_n, fit)) +
  geom_point(
    data = agg_dat_large %>% 
      filter(week_n %in% newdata_large$week_n,
             !strata == "Saanich"),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "red") +
  facet_grid(stock_group~strata) +
  coord_cartesian(xlim = c(24, 40), ylim = c(0, 1)) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top")

png(
  here::here("figs", "stock_comp_fishery", "smooth_preds_large.png"),
  height = 5, width = 5, units = "in", res = 250
)
large_pred_smooth
dev.off()


## compare predictions from all 3 models
fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_mvtw.rds"
  )
)
fit_slot <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_slot.rds"
  )
)
fit_large <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_large.rds"
  )
)
fit_list <- list(fit, fit_slot, fit_large)
model_names <- c("full", "slot", "large")

excl <- grepl("year", gratia::smooths(fit))
yr_coefs <- gratia::smooths(fit)[excl]

pred_list <- purrr::map(
  fit_list, 
  ~ pred_dummy(
    .x,
    se.fit = TRUE,
    category_name = "stock_group",
    origdata = .x$model,
    newdata = newdata_slot,
    exclude = yr_coefs
    )
)

new_dat <- purrr::map2(
  pred_list, model_names,
  ~ cbind(newdata_slot, fit = .x$fit, se.fit = .x$se.fit) %>% 
    mutate(
      model = .y,
      lower = fit + (qnorm(0.025)*se.fit),
      upper = fit + (qnorm(0.975)*se.fit)
    ) %>%
    filter(!strata == "Saanich",
           slot_limit == "yes")
) %>% 
  bind_rows()


model_comp_smooth <- ggplot(new_dat, aes(week_n, fit, colour = model)) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper, fill = model), alpha = 0.5) +
  facet_grid(stock_group~strata) +
  coord_cartesian(xlim = c(25, 38), ylim = c(0, 0.8)) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = "Sample\nSize") +
  theme(legend.position = "top",
        axis.title.x = element_blank()) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37),
    labels = c("Jun", "Jul", "Aug", "Sep"),
    expand = c(0, 0)
  ) 

png(
  here::here("figs", "stock_comp_fishery", "model_comp.png"),
  height = 6.5, width = 6.5, units = "in", res = 250
)
model_comp_smooth
dev.off()