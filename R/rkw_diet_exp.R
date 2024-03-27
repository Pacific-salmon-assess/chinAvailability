## Diet Explore
# Clean and begin analysis of SRKW diet data
# 1) Raw data figs:
# - map of sample collection locations
# - temporal sampling coverage bar plots
# - stock comp bar plots
# 2) Model fitting and figs
# - smooth predictions
# - stacked ribbon plots



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
    strata = ifelse(is.na(strata), "swiftsure", as.character(strata)) %>%
      factor(
      .,
      levels = c("swiftsure", "swiftsure_nearshore", "renfrew", "cJDF",
                 "sooke", "vic"),
      labels = c("Swiftsure", "Nitinat", "Renfrew", "cJDF",
                 "Sooke/\nVictoria", "San Juan\nIslands")
    ),
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
    era = ifelse(year < 2015, "early", "current") %>% 
      fct_relevel(., "current", after = Inf),
    # sampling event = all samples collected in a given strata-year-week
    sample_id = paste(year, week, strata, sep = "_")
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("longitude", "latitude"), ll_crs = 4326, units = "m"
  ) %>% 
  select(
    id, sample_id, era, year, month, week, strata, utm_y = Y, utm_x = X,
    stock, stock_prob, stock_group, lat = latitude, lon = longitude
  )

# calculate mean location of each event
sample_key <- dat %>% 
  select(sample_id, id, utm_y, utm_x) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarize(
    n_samples = length(unique(id)),
    utm_y = mean(utm_y),
    utm_x = mean(utm_x)
  )

ppn_dat <- dat %>% 
  group_by(sample_id, era, year, week, strata, stock_group) %>% 
  summarize(
    agg_count = sum(stock_prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key, by = "sample_id") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 
saveRDS(
  ppn_dat, here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
)


# SMU colour palette
smu_colour_pal <- c("grey30", "#08306B", "#6A51A3", "#CBC9E2", "#67000D", 
                    "#A50F15", "#EF3B2C", "#FC9272", "#FCBBA1")
names(smu_colour_pal) <- levels(dat$stock_group)

era_pal <- c(15, 16)
names(era_pal) <- levels(ppn_dat)



## RAW DATA FIGURES ------------------------------------------------------------

# sample coverage through time and among strata
diet_samp_cov <- dat %>% 
  group_by(week, strata, year, era) %>% 
  summarize(n = length(unique(id))) %>% 
  ggplot(.) +
  geom_point(aes(x = week, y = year, size = n, shape = era), 
             alpha = 0.6) +
  facet_wrap(~strata) +
  scale_size_continuous(name = "Sample\nSize") +
  scale_shape_manual(values = era_pal, name = "Sample\nEra") +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )


# stacked bar plots
diet_samp_bar <- ggplot(dat) +
  geom_bar(aes(x = month, y = stock_prob, fill = stock_group), 
           stat = "identity") +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = "Stock\nGroup") +
  labs(
    y = "Prey Remains Composition"
  ) +
  scale_x_continuous(
    breaks = c(6, 7, 8, 9, 10),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )


# bubble plots of raw observations (replicated in smooths below)
ggplot(ppn_dat) +
  geom_jitter(aes(x = week, y = agg_prob, size = n_samples, fill = era), 
             alpha = 0.7, shape = 21, width = 0.1) +
  facet_grid(stock_group~strata) +
  ggsidekick::theme_sleek()


## export 
png(
  here::here("figs", "ms_figs", "temporal_sample_coverage.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
diet_samp_cov
dev.off()

png(
  here::here("figs", "ms_figs", "monthly_comp_bar.png"),
  height = 5, width = 8, units = "in", res = 250
)
diet_samp_bar
dev.off()


## FIT MODEL -------------------------------------------------------------------

# add zero observations
agg_dat <- expand.grid(
  sample_id = unique(dat$sample_id),
  stock_group = unique(dat$stock_group)
) %>% 
  # add n_samples separately to ensure missing stocks represented
  left_join(., 
            ppn_dat %>% 
              select(sample_id:strata, utm_y, utm_x, n_samples) %>% 
              distinct(), 
            by = "sample_id") %>%
  left_join(., 
            ppn_dat %>% 
              select(sample_id, stock_group, agg_prob), 
            by = c("sample_id", "stock_group")) %>%
  mutate(
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    era = ifelse(year < 2015, "early", "current"),
    year = as.factor(year),
    stock_group = as.factor(stock_group)
  ) %>% 
  droplevels()


fit <- gam(
  agg_prob ~ era*stock_group + 
    s(week, by = stock_group, k = 4) + 
    s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25), 
  data = agg_dat, family = "tw"
)
class(fit) = c( "mvtweedie", class(fit) )
saveRDS(
  fit,
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_diet_mvtw.rds")
)

# fit2 <- gam(
#   agg_prob ~ stock_group + 
#     s(week, by = stock_group, k = 4) + 
#     s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 25), 
#   data = agg_dat, family = "tw"
# )
# class(fit2) = c( "mvtweedie", class(fit2) )




# loc_key <- data.frame(
#   loc = c("swift", "renfrew", "vic"),
#   lon = c(-124.703, -124.481, -123.189),
#   lat = c(48.541, 48.520, 48.438)
# ) %>%
#   sdmTMB::add_utm_columns(
#     ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "m",
#     utm_names = c("utm_x", "utm_y")
#   )
loc_key <- agg_dat %>% 
  group_by(strata) %>% 
  summarize(utm_x = mean(utm_x),
            utm_y = mean(utm_y))

newdata <- expand.grid(
  strata = unique(agg_dat$strata),
  week = seq(25, 38, by = 0.25),
  stock_group = levels(agg_dat$stock_group),
  era = unique(agg_dat$era),
  year = levels(agg_dat$year)[1]
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  filter(
    !strata == "sooke"
  )

pred = predict(
  fit,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata
)
newdata = cbind( newdata, fit=pred$fit, se.fit=pred$se.fit )
newdata$lower = newdata$fit + (qnorm(0.025)*newdata$se.fit)
newdata$upper = newdata$fit + (qnorm(0.975)*newdata$se.fit)

# pred2 = predict(
#   fit2,
#   se.fit = TRUE,
#   category_name = "stock_group",
#   origdata = agg_dat,
#   newdata = newdata
# )
# newdata = cbind( newdata, fit2=pred2$fit, se.fit2=pred2$se.fit )
# newdata$lower2 = newdata$fit2 + (qnorm(0.025)*newdata$se.fit2)
# newdata$upper2 = newdata$fit2 + (qnorm(0.975)*newdata$se.fit2)
# newdata2 <- newdata %>% filter(era == "current")

# remove datalimited strata
newdata_trim <- newdata %>% 
  filter(
    !strata %in% c("Sooke/\nVictoria", "cJDF"),
    !(era == "current" & strata %in% c("San Juan\nIslands")) 
  )

diet_pred_smooth <- ggplot(newdata_trim,
       aes(week, fit, colour = era, fill = era)) +
  geom_point(
    data = agg_dat %>%
      filter(strata %in% newdata$strata,
             !strata %in% c("Sooke/\nVictoria", "cJDF")),
    aes(x = week, y = agg_prob, size = n_samples),
    alpha = 0.3
  ) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  facet_grid(stock_group ~ strata) +
  coord_cartesian(ylim = c(0,1), xlim = c(25, 38)) +
  labs(
    y = "Predicted Proportion of Diet Sample",
    fill = "Sampling\nEra",
    colour = "Sampling\nEra",
    size = "Sample\nSize"
  ) +
  ggsidekick::theme_sleek() +
  # scale_size_continuous(name = "Sample\nSize") +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) 


diet_pred_stacked <- ggplot(data = newdata_trim, 
       aes(x = week)) +
  geom_area(aes(y = fit, colour = stock_group, fill = stock_group), 
            stat = "identity") +
  scale_fill_manual(name = "Stock Group", values = smu_colour_pal) +
  scale_colour_manual(name = "Stock Group", values = smu_colour_pal) +
  labs(y = "Predicted Mean Composition of Diet Sample", x = "Week") +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "right",
    axis.text = element_text(size=9),
    plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points"),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA)) +
  facet_grid(strata~era) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  )


## export 
png(
  here::here("figs", "ms_figs", "smooth_preds.png"),
  height = 8, width = 6.5, units = "in", res = 250
)
diet_pred_smooth
dev.off()

png(
  here::here("figs", "rkw_diet", "stacked_pred.png"),
  height = 8, width = 6.5, units = "in", res = 250
)
diet_pred_stacked
dev.off()



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
