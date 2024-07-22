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
library(sf)

raw_dat <- readRDS(
  here::here(
    "data", "rkw_diet", "RKW predation_chin samples_long_filtered.RDS"
    )
  )


stock_key <- readRDS(
  here::here("data", "rec", "finalStockList_Jul2024.rds")
  ) %>%
  janitor::clean_names() %>% 
  mutate(
    # adjust stock groups to match priorities
    stock_group = case_when(
      pst_agg == "CR-upper_sp" | region1name == "Willamette_R" ~ "Col_Spring",
      pst_agg %in% c("CR-lower_fa", "CR-lower_sp", "CR-upper_su/fa") ~ 
        "Col_Summer_Fall",
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN", "NEVI") ~
        "ECVI_SOMN",
      pst_agg %in% c("CA_ORCST", "WACST", "Russia", "NBC_SEAK", 
                     "Yukon") ~ "other",
      grepl("Fraser", region1name) ~ region1name,
      TRUE ~ pst_agg
    ) %>% 
      factor(
        .,
        levels = c("other", "Col_Spring", "Col_Summer_Fall", "PSD",  
                   "WCVI", "ECVI_SOMN", "Fraser_Spring_4.2",
                   "Fraser_Spring_5.2", "Fraser_Summer_5.2", 
                   "Fraser_Summer_4.1", "Fraser_Fall"),
        labels = c("other", "Col_Spr", "Col_Sum/Fall", "PSD", "WCVI", 
                   "ECVI_SOMN", "FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2", 
                   "FR_Sum_4.1", "FR_Fall")
      )
  ) %>% 
  select(
    stock, stock_group
  )


dat <- raw_dat %>% 
  # correct weird stock 
  mutate(
    stock = ifelse(stock == "BIGQUL@LANG", "BIG_QUALICUM", stock)
  ) %>% 
  left_join(., stock_key, by = "stock") %>%
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
    # in-fill missing ages based on dominant life history strategy
    total_year = case_when(
      grepl("M", gr_age) & grepl(".2", smu) ~ sw_year + 2,
      grepl("M", gr_age) & !grepl(".2", smu) ~ sw_year + 1,
      TRUE ~ total_year
    ) %>% 
      as.factor(),
    era = ifelse(year < 2015, "early", "current") %>% 
      fct_relevel(., "current", after = Inf),
    # sampling event = all samples collected in a given strata-year-month
    # week not feasible given sample sizes
    sample_id = paste(year, week, strata, sep = "_"),
    sw_age = as.factor(sw_year),
    age_stock_group = case_when(
      grepl("Fraser", smu) ~ smu,
      stock == "CAPILANO" ~ "Fraser_Fall",
      agg == "SOG" ~ "ECVI_SOMN",
      TRUE ~ agg
    ) %>% 
      as.factor()
  ) %>% 
  # since weeks may span multiple months, calculate median month for each week
  group_by(sample_id) %>% 
  mutate(month = median(month)) %>% 
  ungroup() %>% 
  mutate(sample_id_pooled = paste(era, month, strata, sep = "_")) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("longitude", "latitude"), ll_crs = 4326, units = "km"
  ) %>% 
  select(
    id, sample_id, sample_id_pooled,
    fw_year, sw_age, total_year, age_f = total_year,
    era, year, month, week_n = week, strata, utm_y = Y, utm_x = X,
    stock, stock_prob, stock_group, age_stock_group,
    lat = latitude, lon = longitude
  )


## aggregate data (calculate mean location and sample size of each event) 
# at two scales: 
# 1) week-year-strata for modeling
# 2) month-strata for simulation sampling

# week-year scale
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
  group_by(sample_id, era, year, week_n, month, strata, stock_group) %>% 
  summarize(
    agg_count = sum(stock_prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key, by = "sample_id") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 
# saveRDS(
#   ppn_dat, here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
# )


# month scale
sample_key_pooled <- dat %>% 
  select(sample_id_pooled, id, utm_y, utm_x, week_n) %>% 
  distinct() %>% 
  group_by(sample_id_pooled) %>% 
  summarize(
    n_samples = length(unique(id)),
    utm_y = mean(utm_y),
    utm_x = mean(utm_x),
    week_n = mean(week_n)
  )

ppn_dat_pooled <- dat %>% 
  group_by(sample_id_pooled, era, month, strata, stock_group) %>% 
  summarize(
    agg_count = sum(stock_prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key_pooled, by = "sample_id_pooled") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 
# saveRDS(
#   ppn_dat_pooled, here::here("data", "rkw_diet", "cleaned_ppn_dat_pooled.rds")
# )


# calculate hatchery contribution by stock (uses mean values from 
# mvtweedie_fit.R)
# hatchery_dat <- readRDS(here::here("data", "rec", "hatchery_stock_df.rds")) %>% 
#   select(stock_group, origin2, ppn) %>% 
#   ungroup() %>% 
#   full_join(., 
#             dat %>% select(month, stock_prob, stock_group, era, strata),
#             by = "stock_group",
#             relationship = "many-to-many") %>% 
#   arrange(month, strata) %>% 
#   mutate(scaled_prob = stock_prob * ppn)



## CALCULATE MEAN SIZE ---------------------------------------------------------

## use model fit in size_by_stock.R 
size_fit <- readRDS(here::here("data", "rec", "size_at_age_fit.rds"))

size_pred_dat <- dat %>% 
  filter(!is.na(sw_age)) %>% 
  mutate(
    slot_limit = "no"
  ) 

pp <- predict(size_fit, size_pred_dat)

# since each sample includes multiple stock IDs calculate mean size 
size_pred_dat2 <- size_pred_dat %>% 
  mutate(pred_fl = pp) %>% 
  group_by(id) %>% 
  summarize(
    pred_fl = mean(pred_fl)
  ) %>% 
  ungroup

dat2 <- left_join(
  dat, 
  size_pred_dat2 %>% 
    select(id, pred_fl),
  by = "id"
) %>% 
  mutate(
    size_bin = cut(
      pred_fl, 
      breaks = c(-Inf, 651, 751, 851, Inf), 
      labels = c("55-65", "65-75", "75-85", ">85")
    )
  )

dd <- dat2 %>% 
  select(id, strata, era, month, fw_year, sw_age, age_f, pred_fl, size_bin) %>% 
  distinct() 

ppn_dat_size <- dat2 %>% 
  filter(!is.na(size_bin)) %>% 
  group_by(sample_id, era, year, week_n, month, strata, size_bin) %>% 
  summarize(
    agg_count = sum(length(unique(id))),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key, by = "sample_id") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 
saveRDS(
  ppn_dat_size, here::here("data", "rkw_diet", "cleaned_ppn_dat_size.rds")
)


## PALETTES --------------------------------------------------------------------

age_pal <- c(
  "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c"
)
names(age_pal) <- levels(size_fit$model$sw_age)

# SMU colour palette
smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(dat$stock_group)

# era shape palette
era_pal <- c(15, 16)
names(era_pal) <- levels(ppn_dat)

# size colour palette
size_colour_pal <- c("grey30", "#8c510a", "#f6e8c3", "#c7eae5", "#01665e")
names(size_colour_pal) <- c(NA, levels(dd$size_bin))

# hatchery origin colour palette
# hatchery_colour_pal <- c("#006d2c", "#bae4b3", "grey30", "grey60", "#7a0177",
#                          "#fbb4b9")
# names(hatchery_colour_pal) <- levels(hatchery_dat$origin2)


## RAW DATA FIGURES ------------------------------------------------------------

# sample coverage through time and among strata
diet_samp_cov <- dat %>% 
  group_by(week_n, strata, year, era) %>% 
  summarize(n = length(unique(id))) %>% 
  ggplot(.) +
  geom_point(aes(x = week_n, y = year, size = n, shape = era), 
             alpha = 0.6) +
  facet_wrap(~strata) +
  geom_hline(aes(yintercept = 2013), col = "red", lty = 2) +
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

samp_size <- dat %>% 
  group_by(era, strata) %>% 
  summarize(
    n = length(unique(id))
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
  ) +
  geom_text(
    data = samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

age_samp_bar <- ggplot(dd) +
  geom_bar(aes(x = month, fill = sw_age)) +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(
    name = "Marine\nAge", values = age_pal, na.value = "grey60" 
    ) +
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
  )+
  geom_text(
    data = samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

# hatchery_samp_bar <- ggplot(hatchery_dat) +
#   geom_bar(aes(x = month, y = scaled_prob, fill = origin2), stat = "identity") +
#   facet_grid(era~strata) +
#   ggsidekick::theme_sleek() +
#   scale_fill_manual(
#     name = "Hatchery\nContribution", values = hatchery_colour_pal) +
#   labs(
#     y = "Prey Remains Composition"
#   ) +
#   scale_x_continuous(
#     breaks = c(6, 7, 8, 9, 10),
#     labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
#   ) +
#   theme(
#     legend.position = "top",
#     axis.title.x = element_blank()
#   )

# box plots of size
size_samp_bar <- ggplot(dd) +
  geom_bar(aes(x = month, fill = size_bin)) +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(name = "Size\nClass", values = size_colour_pal, 
                    na.value = "grey60" ) +
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
  ) +
  geom_text(
    data = samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )



## export 
png(
  here::here("figs", "rkw_diet", "temporal_sample_coverage.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
diet_samp_cov
dev.off()

png(
  here::here("figs", "rkw_diet", "monthly_comp_bar.png"),
  height = 5, width = 8, units = "in", res = 250
)
diet_samp_bar
dev.off()

png(
  here::here("figs", "rkw_diet", "comp_bar_prey_age.png"),
  height = 5, width = 8, units = "in", res = 250
)
age_samp_bar
dev.off()

png(
  here::here("figs", "rkw_diet", "comp_bar_prey_size.png"),
  height = 5, width = 8, units = "in", res = 250
)
size_samp_bar
dev.off()
# 
# png(
#   here::here("figs", "rkw_diet", "comp_bar_prey_hatchery.png"),
#   height = 5, width = 8, units = "in", res = 250
# )
# hatchery_samp_bar
# dev.off()


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
              select(sample_id, stock_group, agg_count, agg_prob), 
            by = c("sample_id", "stock_group")) %>%
  mutate(
    month_f = as.factor(month),
    agg_count = ifelse(is.na(agg_prob), 0, agg_count),
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    era = ifelse(year < 2015, "early", "current") %>% 
      as.factor(),
    year_n = year,
    year = as.factor(year),
    stock_group = as.factor(stock_group)
  ) %>% 
  droplevels()

# fails to converge
# fit <- gam(
#   agg_count ~ 0 + stock_group * era * strata +
#     s(week_n, by = stock_group, k = 7),
#     # s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 15),
#   data = agg_dat, family = "tw"
# )
# fit2 <- gam(
#   agg_count ~ 0 + stock_group * era  +#* strata +
#     s(week_n, by = stock_group, k = 7) +
#     s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 15),
#   data = agg_dat, family = "tw"
# )
# class(fit) = c( "mvtweedie", class(fit) )
# saveRDS(
#   fit,
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_diet_mvtw.rds")
# )
# fit <- readRDS(
#   here::here("data", "model_fits", "mvtweedie", "fit_spatial_diet_mvtw.rds")
# )


## parameter estimates plot

# era_pars <- broom::tidy(fit, parametric = TRUE, conf.int = TRUE) %>% 
#   filter(grepl("era", term)) %>% 
#   mutate(
#     stock_group = levels(agg_dat$stock_group) %>% 
#       factor(., levels = levels(agg_dat$stock_group))
#   ) 
# 
# era_plot <- ggplot(era_pars) +
#   geom_pointrange(
#     aes(x = stock_group, y = estimate,  ymin = conf.low, ymax = conf.high, 
#         fill = stock_group), 
#     shape = 21) +
#   scale_fill_manual(values = smu_colour_pal) +
#   geom_hline(aes(yintercept = 0), lty = 2) +
#   ggsidekick::theme_sleek() +
#   labs(y = "Sampling Period Effect Size") +
#   theme(
#     legend.position = "none",
#     axis.title.x = element_blank(),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
#   
#   
# png(
#   here::here("figs", "rkw_diet", "sampling_period_effect.png"),
#   height = 3.5, width = 5, units = "in", res = 250
# )
# era_plot
# dev.off()


## CHECK -----------------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_count == 0) / nrow(agg_dat)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
library(sdmTMB)
fit_sdmTMB <- sdmTMB(
  agg_count ~ 0 + stock_group +#*era +
    s(week_n, by = stock_group, k = 6) +
    s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 15),
  data = agg_dat,
  spatial = "off",
  spatiotemporal = "off",
  family = tweedie(link = "log")
)
saveRDS(
  fit_sdmTMB,
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_sdmTMB.rds")
)
fit_sdmTMB <- readRDS(
  here::here("data", "model_fits", "mvtweedie", "fit_spatial_fishery_ri_sdmTMB.rds")
)



# fit spatial model
mesh <- make_mesh(agg_dat, c("utm_x", "utm_y"), n_knots = 50)

fit_sdmTMB2 <-  sdmTMB(
  agg_count ~ 0 + stock_group*era #+
    # s(week_n, by = stock_group, k = 6)
    ,
  data = agg_dat,
  mesh = mesh,
  spatial = "on",
  groups = "stock_group",
  spatiotemporal = "off",
  family = tweedie(link = "log"),
  control = sdmTMBcontrol(
    map = list(
      ln_tau_Z = factor(
        rep(1, times = length(unique(agg_dat$stock_group)))
      )
    )
  ),
  silent = FALSE,
  do_fit = FALSE
)



## PREDICT ---------------------------------------------------------------------

# import location key made in sampling_maps.R with representative locations
# for each strata
loc_key <- readRDS(
  here::here("data", "spatial", "strata_key.rds")
) %>% 
  # match scale of fitted model
  mutate(
    utm_x = utm_x_m / 1000,
    utm_y = utm_y_m / 1000
  )


newdata <- expand.grid(
  strata = levels(agg_dat$strata),
  week_n = seq(25, 38, by = 0.25),
  stock_group = levels(agg_dat$stock_group),
  # era = unique(agg_dat$era),
  year = levels(agg_dat$year)[1]
) %>%
  left_join(., loc_key, by = 'strata') %>% 
  mutate(
    strata = factor(strata, levels = levels(agg_dat$strata))
  ) %>% 
  filter(
    !strata == "sooke",
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
    # remove strata that are sparsely sampled
    !strata %in% c("Sooke/\nVictoria", "cJDF", "San Juan\nIslands")
  )

diet_pred_smooth <- ggplot(newdata_trim,
       aes(week_n, fit
           #, colour = era, fill = era
           )) +
  geom_point(
    data = agg_dat %>%
      filter(strata %in% newdata_trim$strata),
    aes(x = week_n, y = agg_prob, size = n_samples),
    alpha = 0.3
  ) +
  geom_line() +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5) +
  facet_grid(stock_group ~ strata) +
  coord_cartesian(ylim = c(0,1), xlim = c(25, 38)) +
  labs(
    y = "Predicted Proportion of Diet Sample",
    # fill = "Sampling\nEra",
    # colour = "Sampling\nEra",
    size = "Sample\nSize"
  ) +
  ggsidekick::theme_sleek() +
  # scale_size_continuous(name = "Sample\nSize") +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  ) 

## NOTE RENFREW ECVI/SOMN ESTIAMTES UNRELIABLE
diet_pred_stacked <- ggplot(data = newdata_trim, 
       aes(x = week_n)) +
  geom_area(aes(y = fit, colour = stock_group, fill = stock_group), 
            stat = "identity") +
  scale_fill_manual(name = "Stock Group", values = smu_colour_pal) +
  scale_colour_manual(name = "Stock Group", values = smu_colour_pal) +
  labs(y = "Predicted Mean Composition of Diet Sample", x = "Week") +
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
    breaks = c(25, 29, 33, 37, 41),
    labels = c("Jun", "Jul", "Aug", "Sep", "Oct")
  )


## export 
png(
  here::here("figs", "rkw_diet", "smooth_preds.png"),
  height = 8, width = 6.5, units = "in", res = 250
)
diet_pred_smooth
dev.off()

png(
  here::here("figs", "rkw_diet", "stacked_pred.png"),
  height = 4, width = 6.5, units = "in", res = 250
)
diet_pred_stacked
dev.off()


## SIDE BY SIDE ---------------------------------------------------------------


# generate monthly predictions using both models to compare expected comp
# Given high uncertainty for era-specific model, use simplified version then
# check again when more samples added
fit2 <- gam(
  agg_count ~ 0 + strata + s(week_n, by = stock_group, k = 6) +
  s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 15),
  data = agg_dat, family = "tw"
)
class(fit2) = c( "mvtweedie", class(fit2) )


newdata_both1 <- expand.grid(
  strata = levels(agg_dat$strata),
  # exclude Sep/Oct samples due to small sample size
  week_n = c(25, 29, 33),
  stock_group = levels(agg_dat$stock_group),
  era = unique(agg_dat$era),
  year_n = max(agg_dat$year_n)
) %>%
  left_join(., loc_key, by = 'strata') %>%
  mutate(
    strata = factor(strata, levels = levels(agg_dat$strata)),
    month = factor(week_n, labels = c("Jun", "Jul", "Aug"))
  ) %>% 
  filter(
    strata %in% c("Swiftsure", "Nitinat", "Renfrew"),
    era == "current"
  )

pred_rkw = predict(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
  origdata = agg_dat,
  newdata = newdata_both1
)
newdata1 = cbind( newdata_both1, fit=pred_rkw$fit, se.fit=pred_rkw$se.fit )
newdata1$lower = newdata1$fit + (qnorm(0.025)*newdata1$se.fit)
newdata1$upper = newdata1$fit + (qnorm(0.975)*newdata1$se.fit)
newdata1$model <- "srkw"


# import fishery model fit in mvtweedie_fit.R
rec_fit <- readRDS(
  here::here(
    "data", "model_fits", "mvtweedie", "fit_spatial_fishery_yr_s_mvtw.rds"
  )
)

# modified prediction function that allows for exclude 
source(here::here("R", "functions", "pred_mvtweedie2.R"))

# integrate out year smooths
excl <- grepl("year_n", gratia::smooths(rec_fit))
yr_coefs <- gratia::smooths(rec_fit)[excl]
pred_rec = pred_dummy(
  rec_fit,
  se.fit = TRUE,
  category_name = "stock_group",
  # original dataset used to fit fishery model
  origdata = rec_fit$model,
  newdata = newdata_both1,
  exclude = yr_coefs
)

newdata2 = cbind( newdata_both1, fit=pred_rec$fit, se.fit=pred_rec$se.fit )
newdata2$lower = newdata2$fit + (qnorm(0.025)*newdata2$se.fit)
newdata2$upper = newdata2$fit + (qnorm(0.975)*newdata2$se.fit)
newdata2$model <- "rec"

newdata_both <- rbind(newdata1, newdata2)

png(
  here::here("figs", "ms_figs", "paired_predictions.png"),
  height = 8, width = 6.5, units = "in", res = 250
)
ggplot(newdata_both) +
  geom_pointrange(
    aes(x = month, y = fit, ymin = lower, ymax = upper, fill = model),
    shape = 21,
    position = position_dodge(width = 0.5)
  ) +
  facet_grid(stock_group ~ strata, scales =  "free_y") +
  coord_cartesian(ylim = c(0,1)) +
  labs(
    y = "Predicted Proportion of Sample",
    fill = "Data Source"
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )
dev.off()