## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Adaptation of rkw_strata_fit to fit a spatially explicit model using mvtweedie
# and mgcv or glmmtmb packages
# Includes sensitivity analyses to evaluate effects of a) slot limits and b)
# excluding fish < 75 cm

# Set French language option
FRENCH <- TRUE

# Create appropriate figure directories
if (FRENCH) {
  dir.create("figs-french", showWarnings = FALSE)
  dir.create("figs-french/stock_comp_fishery", showWarnings = FALSE)
  fig_dir <- "figs-french"
} else {
  dir.create("figs/stock_comp_fishery", showWarnings = FALSE)
  fig_dir <- "figs"
}

# Translation helper function
tr <- function(english, french) {
  if (FRENCH) french else english
}

# Helper function for figure paths
fig_path <- function(filename) {
  file.path(fig_dir, filename)
}

# Translation function for strata labels
translate_strata <- function(strata_values) {
  if (FRENCH) {
    case_when(
      strata_values == "Swiftsure\nBank" ~ "Banc\nSwiftsure",
      strata_values == "Port\nRenfrew" ~ "Port\nRenfrew",
      strata_values == "Sooke/\nVictoria" ~ "Sooke/\nVictoria",
      strata_values == "S. Gulf\nIslands" ~ "Îles du Golfe\ndu Sud",
      strata_values == "Saanich" ~ "Saanich",
      strata_values == "Juan\nde Fuca" ~ "Juan\nde Fuca",
      strata_values == "San Juan\nIslands" ~ "Îles\nSan Juan",
      strata_values == "Nitinat" ~ "Nitinat",
      TRUE ~ strata_values
    )
  } else {
    strata_values
  }
}

# Translation function for management measures
translate_management <- function(management_values) {
  if (FRENCH) {
    case_when(
      management_values == "non-retention" ~ "non-rétention",
      management_values == "size limit" ~ "limite de taille",
      management_values == "standard" ~ "standard",
      TRUE ~ management_values
    )
  } else {
    management_values
  }
}

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
    !strata == "other",
    # remove Bessette and Duteau stocks due to assignment concerns
    !(grepl("BESS", stock) | grepl("DUTEA", stock))
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
    ),
    stock_group = fct_relevel(stock_group, "PSD", after = 3)#
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "km",
    utm_names = c("utm_x", "utm_y")
  ) %>% 
  droplevels() %>% 
  group_by(sample_id) %>% 
  mutate(sample_id_n = sum(prob)) %>% 
  ungroup()


## export stock key
# stks <- dat %>%
#   select(stock, region1name, stock_group) %>%
#   distinct() %>%
#   arrange(stock_group, region1name)
# write.csv(stks, here::here("data", "stock_key.csv"), row.names = FALSE)


sample_key <- dat %>% 
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x) %>% 
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
      as.factor(),
    week_z = scale(week_n) %>% as.numeric()
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
hatchery_colour_pal <- c("#006d2c", "#bae4b3", "grey30", "grey60", "#7a0177",
                         "#fbb4b9")
names(hatchery_colour_pal) <- levels(dat$origin2)


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

sample_gap_poly <- data.frame(
  x = c(1, 1, 52, 52),
  y = c(2014, 2016, 2016, 2014)
)

# sampling coverage 
rec_samp_cov <- sample_key %>%
  mutate(strata = translate_strata(strata)) %>%
  ggplot(.) +
  geom_polygon(data = sample_gap_poly, aes(x = x, y = y), 
               fill = "red", alpha = 0.3) +
  geom_jitter(aes(x = week_n, y = year, size = sample_id_n),
              alpha = 0.4
  ) +
  facet_wrap(~strata) +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  labs(
    size = tr("Sample Size", "Taille d'échantillon")
  ) +
  scale_x_continuous(
    breaks = c(2, 20, 36, 50),
    labels = c(tr("Jan", "jan"), tr("May", "mai"), tr("Sep", "sep"), tr("Dec", "déc"))
  ) + 
  lims(y = c(2003, 2023.5)) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )


## supplementary figure showing locations/times of closures
## week_n == 25 is June 15, == 29 is July 15, == 31 is Aug 1, == 36 is Aug 31
sample_key2 <- sample_key %>%
  mutate(
    restricted = case_when(
      ## Vic area strata
      # size restrictions before July 15 and prior to 2019
      strata == "Sooke/\nVictoria" & year == 2018 & week_n < 31 ~ "size limit",
      strata == "Sooke/\nVictoria" & year < 2018 & week_n < 29 ~ "size limit",
      # closed before Aug 1 and prior to 2019
      strata %in% c("Sooke/\nVictoria", "Saanich", "S. Gulf\nIslands") & 
        year >= 2021 & week_n < 31 ~ "non-retention",
      strata %in% c("Sooke/\nVictoria", "Saanich", "S. Gulf\nIslands") & 
        year >= 2019 & week_n < 29 ~ "non-retention",
      strata %in% c("Sooke/\nVictoria", "Saanich", "S. Gulf\nIslands") & 
        year > 2019 & week_n < 36 ~ "size limit",
      ## Swiftsure/Nitinat/Renfrew strata
      # unrestricted prior to 2019
      strata %in% c("Swiftsure\nBank", "Port\nRenfrew") & year < 2019 ~ 
        "standard",
      # closed prior to July 15 with 80 cm slot limit in place until Aug 1
      strata %in% c("Swiftsure\nBank", "Port\nRenfrew") & year >= 2019 & 
        week_n < 28 ~ "non-retention",
      strata %in% c("Swiftsure\nBank", "Port\nRenfrew") & year >= 2019 & 
        week_n < 32 ~ "size limit",
      ## S Gulf Islands area strata
      # size restrictions before July 15 and prior to 2019
      # strata == "S. Gulf\nIslands" & year < 2019 & (week_n > 23 & week_n < 29) ~ 
      #   "size limit",
      # closed before Aug 1 and prior to 2019
      # strata == "S. Gulf\nIslands" & year >= 2019 & week_n < 29 ~
      #   "non-retention",
      strata == "S. Gulf\nIslands" & year == 2018 &  
        (week_n > 26 & week_n < 32) ~ "size limit",
      strata == "S. Gulf\nIslands" & year == 2017 & 
        (week_n > 28 & week_n < 32) ~ "size limit",
      # strata == "S. Gulf\nIslands" & year >= 2019 & week_n < 32 ~ "size limit",
      ## S Gulf Islands area strata
      # size restrictions before July 15 and prior to 2019
      # strata == "Saanich" & year < 2019 & (week_n > 23 & week_n < 29) ~ 
      #   "size limit",
      # closed before Aug 1 and prior to 2019
      strata == "Saanich" & year >= 2019 & week_n < 29 ~
        "non-retention",
      strata == "Saanich" & year >= 2019 & week_n < 32 ~ "size limit",
      TRUE ~ "standard"
    ) 
  ) %>% 
  filter(week_n > 17 & week_n < 41)

# Management measures color palette - consistent colors for English/French
management_colour_pal <- c("#e31a1c", "#33a02c", "#1f78b4")
# Set palette names based on language setting
if (FRENCH) {
  names(management_colour_pal) <- translate_management(
    c("non-retention", "size limit", "standard")
    )
} else {
  names(management_colour_pal) <- c("non-retention", "size limit", "standard")
}

rec_samp_cov2 <- sample_key2 %>%
  mutate(restricted = translate_management(restricted)) %>%
  ggplot(.) +
  geom_jitter(aes(x = week_n, y = year, size = sample_id_n, colour = restricted),
              alpha = 0.4
  ) +
  facet_wrap(~strata, labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  scale_colour_manual(values = management_colour_pal, name = tr("Management\nMeasures", "Mesures de\ngestion")) +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  labs(
    colour = tr("Management Measures", "Mesures de gestion"),
    size = tr("Sample Size", "Taille d'échantillon")
  ) +
  scale_x_continuous(
    breaks = c(19.5, 24, 28.5, 32.5, 37),
    labels = c(tr("May 15", "15 mai"), tr("Jun 15", "15 juin"), tr("Jul 15", "15 juil"), tr("Aug 15", "15 août"), tr("Sep 15", "15 sep"))
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank(),
    axis.text = element_text(size = rel(0.6))
  )


# stacked bar plot
rec_samp_bar <- dat %>%
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = stock_group), 
           stat = "identity") +
  facet_wrap(~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = tr("Stock", "Stock")) +
  labs(
    y = tr("Recreational Fishery Composition\n(Individual Samples)", "Composition de la pêche récréative\n(Échantillons individuels)"),
    fill = tr("Stock", "Stock")
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 9, 12),
    labels = c(tr("Jan", "jan"), tr("May", "mai"), tr("Sep", "sep"), tr("Dec", "déc"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  geom_text(
    data = full_samp_size, 
    aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

# subset of monthly samples that matches RKW diet
rec_samp_bar_summer <- dat %>% 
  filter(month_n %in% c("5", "6", "7", "8", "9", "10")) %>% 
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = stock_group), 
           stat = "identity") +
  facet_wrap(~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = tr("Stock", "Stock")) +
  labs(
    y = tr("Recreational Fishery Composition\n(Individual Samples)", "Composition de la pêche récréative\n(Échantillons individuels)"),
    fill = tr("Stock", "Stock")
  ) +
  scale_x_continuous(
    breaks = c(5, 6, 7, 8, 9, 10),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  geom_text(
    data = summer_samp_size, 
    aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )


# as above but for hatchery origin
rec_samp_bar_h <- dat %>%
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = origin2),
           stat = "identity") +
  facet_wrap(~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = hatchery_colour_pal, name = tr("Hatchery\nOrigin", "Origine\nÉcloserie")) +
  labs(
    y = tr("Recreational Fishery Composition\n(Individual Samples)", "Composition de la pêche récréative\n(Échantillons individuels)"),
    fill = tr("Hatchery Origin", "Origine d'écloserie")
  ) +
  scale_x_continuous(
    breaks = c(1, 5, 9, 12),
    labels = c(tr("Jan", "jan"), tr("May", "mai"), tr("Sep", "sep"), tr("Dec", "déc"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  )+
  geom_text(
    data = full_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  ) 

# subset of monthly samples that matches RKW diet
rec_samp_bar_summer_h <- dat %>%
  filter(month_n %in% c("5", "6", "7", "8", "9", "10")) %>%
  ggplot(.) +
  geom_bar(aes(x = month_n, y = prob, fill = origin2),
           stat = "identity") +
  facet_wrap(~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = hatchery_colour_pal, 
                    name = tr("Hatchery\nOrigin", "Origine\nÉcloserie")) +
  labs(
    y = tr("Recreational Fishery Composition\n(Individual Samples)", 
           "Composition de la pêche récréative\n(Échantillons individuels)"),
    fill = tr("Hatchery Origin", "Origine d'écloserie")
  ) +
  scale_x_continuous(
    breaks = c(5, 6, 7, 8, 9, 10),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  geom_text(
    data = summer_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )


# proportion hatchery within each stock group (also constrained to post2019 with
# similar patterns)
hatchery_stock_bar <- dat %>%
  group_by(stock_group, origin2) %>%
  summarize(
    sum_prob = sum(prob), #/ total_n,
    .groups = "drop"
  ) %>%
  group_by(stock_group) %>%
  mutate(
    ppn = sum_prob / sum(sum_prob)
  ) %>% 
  ggplot(.) +
  geom_bar(aes(x = stock_group, y = ppn, fill = origin2),
           stat = "identity") +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = hatchery_colour_pal, name = tr("Hatchery\nOrigin", "Origine\nÉcloserie")) +
  labs(
    y = tr("Recreational Fishery Composition\n(Proportion of Individual Samples)", "Composition de la pêche récréative\n(Proportion des échantillons individuels)"),
    fill = tr("Hatchery Origin", "Origine d'écloserie")
  ) +
 theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text.x = element_text(size = 9, angle = 45, hjust = 1)
  ) 


# stock-specific size
median_size <- dat %>% 
  mutate(fl = fl / 10) %>% 
  group_by(stock_group) %>% 
  summarise(
    median_value = median(fl),
    lower_quantile = quantile(fl, 0.1),
    upper_quantile = quantile(fl, 0.9),
    cv = sd(fl) / mean(fl)
  )

size_density <- ggplot() +
  geom_density(data = dat %>%
                 filter(month_n %in% c("5", "6", "7", "8", "9", "10")), 
               aes(x=fl / 10, group=stock_group, fill=stock_group),
               adjust=1.5) +
  geom_vline(data = median_size, aes(xintercept = median_value), lty = 2) +
  geom_vline(data = median_size, aes(xintercept = lower_quantile), lty = 3) +
  geom_vline(data = median_size, aes(xintercept = upper_quantile), lty = 3) +
  facet_wrap(~stock_group) +
  scale_fill_manual(values = smu_colour_pal) +
  ggsidekick::theme_sleek() +
  labs(x = tr("Fork Length (cm)", "Longueur à la fourche (cm)")) +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    axis.ticks=element_blank(),
    axis.text.y = element_blank(),
    axis.title.y = element_blank()
  )



## export
png(
  fig_path(file.path("stock_comp_fishery", "rec_temporal_sample_coverage.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_cov
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_temporal_sample_coverage_management.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_cov2
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_monthly_comp_bar.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_monthly_comp_bar_summer.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_summer
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_monthly_hatchery_comp_bar.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_h
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_monthly_hatchery_comp_bar_summer.png")),
  height = 5, width = 7.5, units = "in", res = 250
)
rec_samp_bar_summer_h
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "rec_stock_hatchery_bar.png")),
  height = 5, width = 8.25, units = "in", res = 250
)
hatchery_stock_bar
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "size_density.png")),
  height = 5, width = 8.25, units = "in", res = 250
)
size_density
dev.off()


## FIT MODEL -------------------------------------------------------------------

# Includes year/stock as RE; remove global smooth after sdmTMB v fails to 
# converge
# system.time(
#   fit2 <- gam(
#     agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 20, bs = "cc") +
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 35) +
#       s(sg_year, bs = "re"),
#     data = agg_dat, family = "tw", method = "REML",
#     knots = list(week_n = c(0, 52))
#   )
# )
# class(fit2) = c( "mvtweedie", class(fit2) )
# saveRDS(
#   fit2,
#   here::here("data", "model_fits", "fit_spatial_fishery_ri_mvtw.rds")
# )

fit2 <- readRDS(
  here::here(
    "data", "model_fits", "fit_spatial_fishery_ri_mvtw.rds")
)


## CHECK -----------------------------------------------------------------------

ppn_zero_obs <- sum(agg_dat$agg_prob == 0) / nrow(agg_dat)


# predictions from equivalent models are identical in both families
# fit_m <- gam(
#   agg_prob ~ 0 + stock_group + s(week_n, k = 4, bs = "cc") +
#     s(sg_year, bs = "re"),
#   data = agg_dat, family = "tw", method = "REML")
# fit_s <- sdmTMB(
#   agg_prob ~ 0 + stock_group +  s(week_n, k = 4, bs = "cc") + 
#     (1 | sg_year),
#   data = agg_dat,
#   spatial = "off",
#   spatiotemporal = "off",
#   family = tweedie(link = "log"))
# p_m <- predict(fit_m)
# p_s <- predict(fit_s)
# plot(p_m ~ p_s$est)
# abline(0, 1)


# simulate by fitting sdmTMB equivalent of univariate Tweedie
library(sdmTMB)
fit_sdmTMB <- sdmTMB(
  agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 20, bs = "cc") +
    # s(utm_y, utm_x, m = c(0.5, 1), bs = "ds", k = 15) +
    s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 35) +
    (1 | sg_year)
  ,
  data = agg_dat,
  spatial = "off",
  spatiotemporal = "off",
  family = tweedie(link = "log"),
  knots = list(week_n = c(0, 52))
)

saveRDS(
  fit_sdmTMB,
  here::here("data", "model_fits", "fit_spatial_fishery_ri_sdmTMB.rds")
)
fit_sdmTMB <- readRDS(
  here::here("data", "model_fits", "fit_spatial_fishery_ri_sdmTMB.rds")
)


# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
s_sdmTMB <- simulate(fit_sdmTMB, nsim = 500)
sum(s_sdmTMB == 0) / length(s_sdmTMB)

samp <- sdmTMBextra::predict_mle_mcmc(
  fit_sdmTMB, mcmc_iter = 101, mcmc_warmup = 100
  )
mcmc_res <- residuals(fit_sdmTMB, type = "mle-mcmc", mcmc_samples = samp)

png(
  fig_path(file.path("stock_comp_fishery", "qq_plot_stock.png")),
  height = 4.5, width = 4.5, units = "in", res = 250
)
qqnorm(mcmc_res,
       xlab = tr("Theoretical Quantiles", "Quantiles théoriques"),
       ylab = tr("Sample Quantiles", "Quantiles échantillonnaux")
       )
qqline(mcmc_res)
dev.off()


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
      week_n, breaks = 5, 
      labels = c(tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), 
                 tr("Sep", "sep"), tr("Oct", "oct")))
  )

colnames(s_sdmTMB) <- paste("sim", seq(1:500), sep = "_")
sim_comp <- cbind(agg_dat, s_sdmTMB) %>% 
  pivot_longer(cols = starts_with("sim_"), names_to = "sim_number", 
               values_to = "conc") 

avg_sim_comp1 <- sim_comp %>% 
  group_by(sim_number, sample_id) %>% 
  mutate(
    # calculate number of fish simulated per simulation and sampling event
    sim_sample_conc = sum(conc),
    sim_ppn = ifelse(sim_sample_conc == "0", 0, conc / sim_sample_conc)
  ) %>% 
  filter(week_n %in% sum_samps$week_n) %>% 
  left_join(., week_key, by = "week_n") 
avg_sim_comp <- avg_sim_comp1 %>% 
  group_by(strata, month, stock_group, sim_number) %>% 
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
  group_by(strata, month, strata_n, stock_group) %>% 
  summarize(
    mean_obs_ppn = mean(agg_ppn),
    se_obs_ppn = sqrt(mean_obs_ppn * (1 - mean_obs_ppn) / strata_n),
    lo = pmax(0, mean_obs_ppn - (1.96 * se_obs_ppn)),
    up = pmin(1, mean_obs_ppn + (1.96 * se_obs_ppn))
  ) %>% 
  distinct()

# manage size by plotting 4 at a time
g1 <- c("other", "Col_Spr", "Col_Sum/Fall", "PSD")
g2 <- c("WCVI", "ECVI_SOMN", "FR_Spr_4sub2", "FR_Spr_5sub2")
g3 <- c("FR_Sum_5sub2", "FR_Sum_4sub1", "FR_Fall")
post_sim_list <- purrr::map(
  list(g1, g2, g3), function (x) {
    ggplot() +
      geom_boxplot(
        data = avg_sim_comp %>% 
          filter(stock_group %in% x) %>% 
          droplevels(),
        aes(x = month, y = mean_sim_ppn)
      ) +
      geom_pointrange(
        data = avg_obs_comp %>% 
          filter(stock_group %in% x) %>%
          droplevels(),
        aes(x = as.numeric(month) - 0.2, y = mean_obs_ppn, ymin = lo, ymax = up), 
        col = "red", alpha = 0.6) +
      ggsidekick::theme_sleek() +
      facet_grid(strata~stock_group) +
      labs(y = tr("Mean Composition", "Composition moyenne")) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.title.x = element_blank()) +
      geom_text(
        data = summer_samp_size, aes(x = Inf, y = Inf, label = paste(n)),
        hjust = 1.1, vjust = 1.1
      )
  }
)


png(
  fig_path(file.path("stock_comp_fishery", "posterior_sims_stock1.png")),
  height = 7.5, width = 7.5, units = "in", res = 250
)
post_sim_list[[1]]
dev.off()
png(
  fig_path(file.path("stock_comp_fishery", "posterior_sims_stock2.png")),
  height = 7.5, width = 7.5, units = "in", res = 250
)
post_sim_list[[2]]
dev.off()
png(
  fig_path(file.path("stock_comp_fishery", "posterior_sims_stock3.png")),
  height = 7.5, width = 7.5, units = "in", res = 250
)
post_sim_list[[3]]
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
    upper = fit + (qnorm(0.975)*se.fit),
    stock_group2 = gsub("_", "\n", stock_group)
  ) %>% 
  filter(
    !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & 
        week_n < 22),
    !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & 
        week_n > 40)
  )

year_preds <- newdata_yr %>%
  ggplot(., aes(week_n, fit)) +
  geom_line(aes(colour = year)) +
  facet_grid(stock_group2~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  scale_colour_discrete() +
  coord_cartesian(xlim = c(20, 43)
  ) +
  labs(y=tr("Predicted Proportion", "Proportion prédite"),
       colour = tr("Year", "Année")) +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  labs(
    size = tr("Sample Size", "Taille d'échantillon")
  ) +
  theme(legend.position = "top",
        axis.title.x = element_blank()) +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), 
               tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
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
  scale_fill_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
  scale_colour_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
  labs(
    fill = tr("Stock Group", "Groupe de stocks"),
    colour = tr("Stock Group", "Groupe de stocks")
  ) +
  coord_cartesian(xlim = c(20, 43)#, ylim = c(0, 1)
                  ) +
  labs(y=tr("Predicted Proportion", "Proportion prédite"), 
       x = tr("Sampling Week", "Semaine d'échantillonnage")) +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  labs(
    size = tr("Sample Size", "Taille d'échantillon")
  ) +
  theme(legend.position = "top") +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
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
summer_preds <- newdata3b %>%
  ggplot(., aes(week_n, fit)) +
  geom_point(
    data = agg_dat %>% 
      filter(week_n %in% newdata$week_n,
             !strata == "Saanich"
             ),
    aes(x = week_n, y = agg_ppn, size = sample_id_n),
    alpha = 0.3
  ) +
  geom_line(colour = "red") +
  geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.5, fill = "red") +
  facet_grid(stock_group~strata, 
             labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  coord_cartesian(xlim = c(20, 43), ylim = c(0, 1)) +
  labs(y=tr("Predicted Proportion", "Proportion prédite"), 
       x = tr("Sampling Week", "Semaine d'échantillonnage"),
       size = tr("Sample Size", "Taille d'échantillon")) +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  theme(legend.position = "top",
        strip.text = element_text(size = 7)) +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
  ) 

summer_preds_fullx <- summer_preds +
  coord_cartesian(xlim = c(0, 52)) +
  scale_x_continuous()


## stacked ribbon predictions
summer_pred_stacked <- newdata3b %>% 
    filter(week_n > 19 & week_n < 44,
           !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & week_n < 22),
           !(strata %in% c("Swiftsure\nBank", "Nitinat", "Port\nRenfrew") & week_n > 40)) %>%
  ggplot(., aes(x = week_n)) +
  geom_area(aes(y = fit, colour = stock_group, fill = stock_group), 
            stat = "identity") +
  scale_fill_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
  scale_colour_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
  labs(y = tr("Predicted Mean Composition of Fishery Sample", 
              "Composition moyenne prédite de l'échantillon de pêche"), 
       x = tr("Week", "Semaine"),
       fill = tr("Stock Group", "Groupe de stocks"),
       colour = tr("Stock Group", "Groupe de stocks")) +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "none",
    axis.text = element_text(size=9),
    plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points"),
    axis.title.x = element_blank()
  ) +
  coord_cartesian(expand = FALSE, ylim = c(0, NA)) +
  facet_wrap(~strata, labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  scale_x_continuous(
    breaks = c(20.75, 25, 29.25, 33.5, 38, 42.5),
    labels = c(tr("May", "mai"), tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep"), tr("Oct", "oct"))
  ) 
summer_pred_legend <- cowplot::get_legend(
  ggplot(
    data = newdata3b %>% filter(strata== "Swiftsure\nBank"), aes(x = week_n)
  ) +
    geom_area(aes(y = fit, colour = stock_group, fill = stock_group), 
              stat = "identity") +
    scale_fill_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
    scale_colour_manual(name = tr("Stock Group", "Groupe de stocks"), values = smu_colour_pal) +
    guides(fill = guide_legend(ncol = 2),
           colour = guide_legend(ncol = 2))
)


png(
  fig_path(file.path("stock_comp_fishery", "smooth_preds_chinook_year.png")),
  height = 7.5, width = 6.5, units = "in", res = 250
)
year_preds
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "season_preds_chinook.png")),
  height = 6.5, width = 7.5, units = "in", res = 250
)
season_preds
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "smooth_preds_chinook.png")),
  height = 8.5, width = 7.5, units = "in", res = 250
)
summer_preds
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "smooth_preds_chinook_xaxis.png")),
  height = 8.5, width = 7.5, units = "in", res = 250
)
summer_preds_fullx
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "smooth_preds_chinook_stacked.png")),
  height = 6.5, width = 8.5, units = "in", res = 250
)
cowplot::ggdraw() +
  cowplot::draw_plot(summer_pred_stacked) +
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
  ) 

excl3 <- #grepl("week_n", gratia::smooths(fit2)) | 
  grepl("year", gratia::smooths(fit2))
yr_coefs <- gratia::smooths(fit2)[excl3]
pred_sp <- pred_dummy(
  fit2,
  se.fit = TRUE,
  category_name = "stock_group",
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
# 3) calculate scaled_fit for each stock_group
new_dat_sp_plot <- cbind(
  new_dat_sp, fit=pred_sp$fit#, se.fit=pred_sp$se.fit 
) %>% 
  group_by(stock_group, X, Y, utm_y, utm_x) %>% 
  summarize(
    mean_fit = mean(fit)
  ) %>% 
  group_by(X, Y, utm_y, utm_x) %>% 
  mutate(fit = mean_fit / sum(mean_fit)) %>% 
  group_by(stock_group) %>% 
  mutate(
    scaled_fit = fit / max(fit)
  ) %>% 
  ungroup() 

spatial_pred <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_wrap(
    ~ stock_group, ncol = 3
  ) +
  scale_fill_viridis_c(
    name = tr("Predicted Proportion\nof Rec Catch", "Proportion prédite\ndes prises récréatives")
  ) +
  labs(
    fill = tr("Predicted Proportion of Rec Catch", "Proportion prédite des prises récréatives")
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 9)
  ) 
spatial_legend <- cowplot::get_legend(
  ggplot() +
    geom_raster(data = new_dat_sp_plot, 
                aes(x = X, y = Y, fill = fit)) +
    geom_sf(data = coast, color = "black", fill = "grey") +
    scale_fill_viridis_c(
      name = tr("Predicted Proportion\nof Rec Catch", "Proportion prédite\ndes prises récréatives")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(
      legend.key.size = unit(0.5, "lines"),  # Adjust the size of the legend keys
      legend.text = element_text(size = 8)   # Adjust the text size of the legend
    )
)


spatial_pred_scaled <- ggplot() +
  geom_raster(data = new_dat_sp_plot, 
              aes(x = X, y = Y, fill = scaled_fit)) +
  geom_sf(data = coast, color = "black", fill = "grey") +
  facet_wrap(
    ~ stock_group, ncol = 3
  ) +
  scale_fill_viridis_c(
    option = "A",
    name = tr("Predicted Scaled\nProportion\nof Rec Catch", 
              "Proportion mise à l'échelle prédite\ndes prises récréatives")
  ) +
  labs(
    fill = tr("Predicted Scaled Proportion of Rec Catch", 
              "Proportion mise à l'échelle prédite des prises récréatives")
  ) +
  ggsidekick::theme_sleek()  +
  theme(
    axis.title = element_blank(),
    axis.text = element_blank(), 
    axis.ticks = element_blank(),
    legend.position = "none",
    strip.text = element_text(size = 9)
  ) 
spatial_scaled_legend <- cowplot::get_legend(
  ggplot() +
    geom_raster(data = new_dat_sp_plot, 
                aes(x = X, y = Y, fill = scaled_fit)) +
    geom_sf(data = coast, color = "black", fill = "grey") +
    scale_fill_viridis_c(
      option = "A",
      name = tr("Predicted Scaled Proportion\nof Rec Catch", 
                "Proportion mise à l'échelle prédite\ndes prises récréatives")
    ) +
    guides(fill = guide_legend(nrow = 1)) +
    theme(
      legend.key.size = unit(0.5, "lines"),  # Adjust the size of the legend keys
      legend.text = element_text(size = 8)   # Adjust the text size of the legend
    )
)

# spatial_pred_se <- ggplot() +
#   geom_raster(data = new_dat_sp_plot, 
#               aes(x = X, y = Y, fill = se.fit)) +
#   geom_sf(data = coast, color = "black", fill = "grey") +
#   facet_wrap(
#     ~ stock_group
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
  fig_path(file.path("stock_comp_fishery", "spatial_preds.png")),
  height = 4.5, width = 7.5, units = "in", res = 250
)
cowplot::ggdraw() +
  cowplot::draw_plot(spatial_pred) +
  cowplot::draw_plot(spatial_legend,
                     height = 0.17, x = 0.33, y = 0.05)
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "spatial_preds_scaled.png")),
  height = 4.5, width = 7.5, units = "in", res = 250
)
cowplot::ggdraw(spatial_pred_scaled) +
  cowplot::draw_plot(spatial_scaled_legend,
                     height = 0.17, x = 0.33, y = 0.05
                     )
dev.off()

# png(
#   fig_path(file.path("stock_comp_fishery", "spatial_preds_se.png")),
#   height = 4, width = 6, units = "in", res = 250
# )
# spatial_pred_se
# dev.off()


## SENSITIVITY ANALYSES --------------------------------------------------------

# 1) Evaluate impact of slot limit by fitting model only to data west of Sooke
# and pre-2019 (REMOVE UNLESS MODIFIED)
# 2) Evaluate impact of size preference of SRKW by fitting model only to samples
# greater than 750 mm fl

# Kept model converges but so few samples that predictions are functionally 
# identical

newdata <- expand.grid(
  strata = unique(agg_dat$strata),
  week_n = unique(agg_dat$week_n),
  stock_group = levels(agg_dat$stock_group),
  year_n = unique(agg_dat$year_n)#,
  # slot_limit = c("yes", "no")
) %>%
  left_join(., loc_key, by = 'strata') %>%
  mutate(
    year = as.factor(year_n),
    sg_year = paste(stock_group, year, sep = "_") %>%
      as.factor(),
    strata = factor(strata, levels = levels(agg_dat$strata))
  ) #%>%
  # filter(
  #   strata %in% c("Renfrew", "Swiftsure", "Nitinat")
  # )


## Exclude released individuals 
# kept_dat <- dat %>%
#   filter(disposition == "Kept") %>%
#   group_by(sample_id) %>%
#   mutate(sample_id_n = sum(prob)) %>%
#   ungroup()
# sample_key_kept <- kept_dat %>%
#   select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x) %>%
#   distinct()
# 
# agg_dat_kept <- expand.grid(
#   sample_id = unique(kept_dat$sample_id),
#   stock_group = unique(kept_dat$stock_group)
# ) %>% 
#   left_join(., sample_key_kept, by = "sample_id") %>% 
#   left_join(
#     ., 
#     kept_dat %>% 
#       group_by(sample_id, stock_group) %>% 
#       summarize(
#         agg_prob = sum(prob)
#       ) %>% 
#       ungroup(),
#     by = c("sample_id", "stock_group")
#   ) %>% 
#   mutate(
#     agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
#     agg_ppn = agg_prob / sample_id_n,
#     year_n = as.numeric(year),
#     year = as.factor(year),
#     stock_group = as.factor(stock_group),
#     utm_x_m = utm_x * 1000,
#     utm_y_m = utm_y * 1000,
#     sg_year = paste(stock_group, year, sep = "_") %>% as.factor()
#   ) 
# 
# # Includes smooth for year by stock group
# system.time(
#   fit_kept <- gam(
#     agg_prob ~ 0 + stock_group + 
#       s(week_n, by = stock_group, k = 20, bs = "cc") +
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 35) +
#       s(sg_year, bs = "re"),
#     data = agg_dat_kept, family = "tw", method = "REML",
#     knots = list(week_n = c(0, 52))
#   )
# )
# class(fit_kept) = c( "mvtweedie", class(fit_kept) )
# saveRDS(
#   fit_kept,
#   here::here(
#     "data", "model_fits", "fit_kept.rds"
#   )
# )



## Size analysis
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
    sg_year = paste(stock_group, year, sep = "_") %>% as.factor()
  ) 

# Includes smooth for year by stock group
# system.time(
#   fit_large <- gam(
#     agg_prob ~ 0 + stock_group + 
#       s(week_n, by = stock_group, k = 20, bs = "cc") +
#       s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds", k = 35) +
#       s(sg_year, bs = "re"),
#     data = agg_dat_large, family = "tw", method = "REML",
#     knots = list(week_n = c(0, 52))
#   )
# )
# class(fit_large) = c( "mvtweedie", class(fit_large) )
# saveRDS(
#   fit_large,
#   here::here(
#     "data", "model_fits", "fit_large.rds"
#   )
# )


## compare predictions from all 3 models
fit <- readRDS(
  here::here(
    "data", "model_fits", "fit_spatial_fishery_ri_mvtw.rds"
  )
)
# fit_kept <- readRDS(
#   here::here(
#     "data", "model_fits", "fit_kept.rds"
#   )
# )
fit_large <- readRDS(
  here::here(
    "data", "model_fits", "fit_large.rds"
  )
)
fit_list <- list(fit, fit_large)
model_names <- c("full", "large")

excl <- grepl("year", gratia::smooths(fit))
yr_coefs <- gratia::smooths(fit)[excl]

newdata3 <- newdata %>% 
  filter(year == "2014")

pred_list <- purrr::map(
  fit_list, 
  ~ pred_dummy(
    .x,
    se.fit = TRUE,
    category_name = "stock_group",
    origdata = .x$model,
    newdata = newdata3,
    exclude = yr_coefs
    )
)

new_dat <- purrr::map2(
  pred_list, model_names,
  ~ cbind(newdata3, fit = .x$fit, se.fit = .x$se.fit) %>% 
    mutate(
      model = .y,
      lower = fit + (qnorm(0.025)*se.fit),
      upper = fit + (qnorm(0.975)*se.fit)
    ) 
) %>% 
  bind_rows() %>% 
  mutate(
    # model = case_when(
    #   model == "slot" ~ paste(model, slot_limit, sep = " "),
    #   TRUE ~ model
    # ),
    model = factor(
      model,
      levels = c("full", #"slot yes", "slot no", 
                 "large"),
      labels = c("standard",# "post-2019\nmanagement", "pre-2019\nmanagement", 
                 "large")
    )
  ) %>% 
  filter(
    week_n > 23 & week_n < 39,
    !strata == "Saanich"
    # !(model %in% c("full", "large") & slot_limit == "yes")
    )

model_pal <- c("#e41a1c", "#377eb8")#, "#4daf4a", "#984ea3")
names(model_pal) <- levels(new_dat$model)

model_comp_smooth1 <- new_dat %>% 
    filter(!(grepl("FR", stock_group) | stock_group == "ECVI_SOMN")) %>%
  ggplot(., aes(week_n, fit, colour = model)) +
  geom_line() +
  facet_grid(stock_group~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  labs(y=tr("Predicted Proportion", "Proportion prédite"), 
       x = tr("Sampling Week", "Semaine d'échantillonnage"),
       colour = tr("Model", "Modèle")) +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  scale_colour_manual(
    name = tr("Model", "Modèle"),
    values = model_pal
  ) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_text(size = 8)) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37),
    labels = c(tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep")),
    expand = c(0, 0)
  ) 
model_comp_smooth2 <- new_dat %>% 
    filter((grepl("FR", stock_group) | stock_group == "ECVI_SOMN")) %>%
  ggplot(., aes(week_n, fit, colour = model)) +
  geom_line() +
  facet_grid(stock_group~strata, scales = "free_y", labeller = labeller(strata = function(x) if(FRENCH) translate_strata(x) else x)) +
  labs(y=tr("Predicted Proportion", "Proportion prédite"), 
       x = tr("Sampling Week", "Semaine d'échantillonnage"),
       colour = tr("Model", "Modèle")) +
  ggsidekick::theme_sleek() +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  scale_colour_manual(
    name = tr("Model", "Modèle"),
    values = model_pal
  ) +
  theme(legend.position = "none",
        axis.title.x = element_blank(),
        strip.text = element_text(size = 8)) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37),
    labels = c(tr("Jun", "juin"), tr("Jul", "juil"), tr("Aug", "août"), tr("Sep", "sep")),
    expand = c(0, 0)
  ) 


png(
  fig_path(file.path("stock_comp_fishery", "model_comp_stock1.png")),
  height = 6.5, width = 6.5, units = "in", res = 250
)
model_comp_smooth1
dev.off()

png(
  fig_path(file.path("stock_comp_fishery", "model_comp_stock2.png")),
  height = 6.5, width = 6.5, units = "in", res = 250
)
model_comp_smooth2
dev.off()
