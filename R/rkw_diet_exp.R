## Diet Explore
# Look at sample size and coverage of diet data


library(tidyverse)
library(mgcv)
library(mvtweedie)

raw_dat <- readRDS(
  here::here(
    "data", "rkw_diet", "RKW predation_chin samples_long_filtered.RDS"
    )
  ) 

dat <- raw_dat %>% 
  filter(
    #remove samples captured in central VI
    !latitude > 49.15
  ) %>%
  mutate(
    strata = ifelse(is.na(strata), "swiftsure", as.character(strata)),
    month = lubridate::month(date),
    stock_group = case_when(
      agg %in% c("CA_ORCST", "CR-lower_fa", "CR-lower_sp", "CR-upper_su/fa",
                     "WACST", "Russia", "CR-upper_sp", "NBC_SEAK") ~ "other",
      stock == "CAPILANO" | smu %in% c("ECVI", "SOMN") ~ "ECVI_SOMN",
      grepl("Fraser", smu) ~ smu,
      TRUE ~ agg
    ) %>%  factor(
      .,
      levels = c("other", "PSD", "WCVI", "ECVI_SOMN", "Fraser_Spring_4.2",
                 "Fraser_Spring_5.2", "Fraser_Summer_5.2", "Fraser_Summer_4.1",
                 "Fraser_Fall"),
      labels = c("other", "PSD", "WCVI", "ECVI_SOMN", "FR_Spr_4.2",
                 "FR_Spr_5.2", "FR_Sum_5.2", "FR_Sum_4.1",
                 "FR_Fall")
    ),
    era = ifelse(year < 2015, "old", "current")
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("longitude", "latitude"), ll_crs = 4326, units = "m"
  ) %>% 
  mutate(
    sample_id = paste(year, week, X, Y, sep = "_")
  ) %>%
  select(
    id, sample_id, era, year, month, week, strata, utm_y = Y, utm_x = X,
    stock, stock_prob, stock_group, lat = latitude, lon = longitude
  )

sample_key <- dat %>% 
  select(sample_id, year, month, week, utm_y, utm_x
         #,strata
         ) %>% 
  distinct()


colour_pal <- c("grey30", "#08306B", "#6A51A3", "#CBC9E2", "#67000D", "#A50F15",
                "#EF3B2C", "#FC9272", "#FCBBA1")
names(colour_pal) <- levels(dat$stock_group)


coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>%
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.75, ymax = 48.8) %>%
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m"))

  
ggplot() +
  geom_sf(data = coast) +
  geom_point(data = dat %>% 
               select(id, era, utm_y, utm_x) %>% 
               distinct(),
             aes(x = utm_x, y = utm_y, fill = era), 
             alpha = 0.7, shape = 21, width = 0.1) +
  ggsidekick::theme_sleek()


## RAW DATA FIGURES ##

# sample coverage through time and among strata
dat %>% 
  group_by(week, strata, year, era) %>% 
  summarize(n = length(unique(id))) %>% 
  ggplot(.) +
  geom_point(aes(x = week, y = year, size = n, fill = era), 
             alpha = 0.6, shape = 21) +
  facet_wrap(~strata) +
  ggsidekick::theme_sleek()


# stacked bar plots
ggplot(dat) +
  geom_bar(aes(x = month, y = stock_prob, fill = stock_group), 
           stat = "identity") +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = colour_pal)


# bubble plots of raw observations
ppn_dat <- dat %>% 
  mutate(  
    sample_id = paste(year, month, week, strata, sep = "_"),
  ) %>% 
  group_by(sample_id) %>% 
  mutate(n_samples = length(unique(id))) %>% 
  group_by(sample_id, era, year, month, strata, stock_group, n_samples,
           utm_x, utm_y) %>% 
  summarize(
    agg_count = sum(stock_prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 

ggplot(ppn_dat) +
  geom_jitter(aes(x = month, y = agg_prob, size = n_samples, fill = era), 
             alpha = 0.7, shape = 21, width = 0.1) +
  facet_grid(stock_group~strata) +
  ggsidekick::theme_sleek()

ggplot(ppn_dat) +
  geom_boxplot(aes(x = as.factor(month), y = agg_prob, fill = era), 
              alpha = 0.7) +
  facet_grid(stock_group~strata) +
  ggsidekick::theme_sleek()

ggplot(ppn_dat) +
  geom_point(aes(x = utm_x, y = utm_y, size = n_samples, fill = era), 
             alpha = 0.7, shape = 21, width = 0.1)+
  ggsidekick::theme_sleek()


## INITIAL FIT

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
        agg_prob = sum(stock_prob)
      ) %>% 
      ungroup(),
    by = c("sample_id", "stock_group")
  ) %>% 
  mutate(
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    era = ifelse(year < 2015, "old", "current"),
    year = as.factor(year),
    stock_group = as.factor(stock_group)
  ) %>% 
  droplevels()


fit <- gam(
  agg_prob ~ era*stock_group + s(month, by = stock_group, k = 4) + 
    s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 20), 
  data = agg_dat, family = "tw"
)
class(fit) = c( "mvtweedie", class(fit) )


## check predictive grid

pred_grid <- readRDS(
  here::here(
    "data", "spatial", "pred_bathy_grid_utm_no_bark.RDS"
  )
) %>% 
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>%
  sf::st_crop(
    ., 
    xmin = min(agg_dat$utm_x) - 500, 
    ymin = min(agg_dat$utm_y) - 500,
    xmax = max(agg_dat$utm_x) + 500, 
    ymax = max(agg_dat$utm_y) + 500
  ) 

newdata <- expand.grid(
  month = unique(agg_dat$month),
  stock_group = levels(agg_dat$stock_group)#,
  # era = unique(agg_dat$era)
) %>% 
  mutate(
    fac = paste(month, #era, 
                as.factor(stock_group), sep = "_")
  ) %>% 
  distinct() %>% 
  split(., .$fac) %>% 
  purrr::map(., function (x) {
    data.frame(
      sf::st_coordinates(pred_grid[ , 1]),
      month = x$month,
      stock_group = x$stock_group#,
      # era = x$era
    ) 
  }) %>% 
  bind_rows %>% 
  rename(utm_y = Y, utm_x = X)

loc_key <- data.frame(
  loc = c("swift", "renfrew", "vic"),
  lon = c(-124.703, -124.481, -123.189),
  lat = c(48.541, 48.520, 48.438)
) %>%
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "m",
    utm_names = c("utm_x", "utm_y")
  )

newdata <- expand.grid(
  loc = c("swift", "renfrew", "vic"),
  month = unique(agg_dat$month),
  stock_group = levels(agg_dat$stock_group)#,
  # era = unique(agg_dat$era)
  # strata = levels(agg_dat$strata)
) %>%
left_join(., loc_key, by = 'loc')

pred = predict(
  fit,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata
)
newdata = cbind( newdata, fit=pred$fit, se.fit=pred$se.fit )
newdata$lower = newdata$fit - newdata$se.fit
newdata$upper = newdata$fit + newdata$se.fit


# Plot
theme_set(theme_bw())
ggplot(newdata, aes(month, fit#, colour = era
                    )) +
  geom_pointrange(aes(ymin = lower, ymax = upper#, colour = era
                      )) +
  # facet_wrap(vars(stock_group)) +
  facet_grid(stock_group ~ loc) +
  ylim(0,1) +
  labs(y="Predicted proportion")


ggplot() +
  geom_raster(data = newdata, 
              aes(x = utm_x, y = utm_y, fill = fit)) +
  facet_grid(
    month ~ stock_group
  ) +
  scale_fill_viridis_c() +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank()
  ) 
