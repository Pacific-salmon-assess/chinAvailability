## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Adaptation of rkw_strata_fit to fit a spatially explicit model using mvtweedie
# and mgcv or glmmtmb packages


library(tidyverse)
library(mvtweedie)
library(mgcv)
library(glmmTMB)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() 


dat <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    !legal == "sublegal",
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
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x) %>% 
  distinct()


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
    year = as.factor(year),
    stock_group = as.factor(stock_group),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
  ) 


system.time(
  fit <- gam(
    agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 5, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds") + 
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds"), 
    data = agg_dat, family = "tw",
    knots = list(week_n = c(0, 52))
  )
)
# ~850 seconds to converge
class(fit) = c( "mvtweedie", class(fit) )

saveRDS(
  fit,
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery.rds")
)

# system.time(
#   fit2 <- gam(
#     agg_prob ~ 0 + stock_group + 
#       s(week_n, by = stock_group, k = 5, bs = "cc") +
#       # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds") + 
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds"), 
#     data = agg_dat, family = "tw",
#     knots = list(week_n = c(0, 52))
#   )
# )
# class(fit2) = c( "mvtweedie", class(fit) )


# identify a single spatial location for each strata based on mean
loc_key <- dat %>% 
  group_by(strata) %>% 
  summarize(
    lat = mean(lat),
    lon = mean(lon)
  ) %>%
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "km",
    utm_names = c("utm_x", "utm_y")
  ) %>% 
  ungroup()

newdata <- expand.grid(
  strata = unique(dat$strata),
  week_n = unique(agg_dat$week_n),
  stock_group = levels(agg_dat$stock_group)
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  filter(
    !strata == "saanich"
  )

# fit 1 
pred = predict(
  fit,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata
)

# # fit 1 
# pred2 = predict(
#   fit2,
#   se.fit = TRUE,
#   category_name = "stock_group",
#   origdata = agg_dat,
#   newdata = newdata
# )

newdata = cbind( newdata, fit=pred$fit, se.fit=pred$se.fit )
newdata$lower = newdata$fit + (qnorm(0.025)*newdata$se.fit)
newdata$upper = newdata$fit + (qnorm(0.975)*newdata$se.fit)#+ newdata$se.fit
# newdata$fit2 <- pred2$fit
# newdata$se.fit2 <- pred2$se.fit
# newdata$lower2 = newdata$fit2 - newdata$se.fit2
# newdata$upper2 = newdata$fit2 + newdata$se.fit2


summer_preds <- ggplot(newdata, aes(week_n, fit)) +
  geom_point(
    data = agg_dat %>% 
      filter(week_n %in% newdata$week_n,
             !strata == "saanich"),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper, colour = "red"
  ), alpha = 0.5) +
  facet_grid(strata ~ stock_group) +
  ylim(0, 1) +
  xlim(24, 40) +
  labs(y="Predicted Proportion", x = "Sampling Week") +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top")

pdf(
  here::here("figs", "mvtweedie_preds", "summer_preds.pdf"),
  height = 5.5, width = 8.5
)
summer_preds
dev.off()


# ggplot(newdata, aes(week_n, fit2#, colour = era
# )) +
#   geom_point(
#     data = agg_dat %>% filter(week_n %in% newdata$week_n),
#     aes(x = week_n, y = agg_ppn, size = sample_id_n),
#     alpha = 0.3
#   ) +
#   geom_line(colour = "red") +
#   geom_ribbon(aes(ymin = lower2, ymax = upper2, colour = "red"
#   ), alpha = 0.5) +
#   # facet_wrap(vars(stock_group)) +
#   facet_grid(strata ~ stock_group) +
#   ylim(0,1) +
#   xlim(24, 40) +
#   labs(y="Predicted proportion") +
#   ggsidekick::theme_sleek() 



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
    stock_group = factor(
      stock_group, 
      levels = c("other", "PSD", "WCVI", "ECVI_SOMN", "Fraser_Spring_4.2",
                 "Fraser_Spring_5.2", "Fraser_Summer_5.2", "Fraser_Summer_4.1",
                 "Fraser_Fall"),
      labels = c("other", "PSD", "WCVI", "ECVI_SOMN", "FR_Spr_4.2",
                 "FR_Spr_5.2", "FR_Sum_5.2", "FR_Sum_4.1",
                 "FR_Fall")
    )
  )
  
pred_sp <- predict(
  fit,
  se.fit = FALSE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = new_dat_sp
)

new_dat_sp <- new_dat_sp %>% 
  mutate(
    fit = pred_sp$fit
  ) %>% 
  group_by(stock_group) %>% 
  mutate(
    scaled_fit = fit / max(fit)
  ) %>% 
  ungroup()

spatial_pred <- ggplot() +
  geom_raster(data = new_dat_sp, 
              aes(x = utm_x, y = utm_y, fill = fit)) +
  facet_grid(
    stock_group ~ week_n
  ) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    legend.position = "top"
  ) 


spatial_pred_scaled <- ggplot() +
  geom_raster(data = new_dat_sp, 
              aes(x = utm_x, y = utm_y, fill = scaled_fit)) +
  facet_grid(
    stock_group ~ week_n
  ) +
  scale_fill_viridis_c(option = "A") +
  ggsidekick::theme_sleek()  +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    legend.position = "top"
  ) 



pdf(
  here::here("figs", "mvtweedie_preds", "spatial_preds.pdf"),
  height = 7.5, width = 9.5
)
spatial_pred
dev.off()

pdf(
  here::here("figs", "mvtweedie_preds", "spatial_preds_scaled.pdf"),
  height = 7.5, width = 9.5
)
spatial_pred_scaled
dev.off()