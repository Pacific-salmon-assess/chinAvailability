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
  here::here("data", "rec", "finalStockList_Sep2024.rds")
  ) %>%
  janitor::clean_names() %>% 
  mutate(
    # adjust stock groups to match priorities
    stock_group = case_when(
      pst_agg %in% c("CR-lower_sp", "CR-upper_sp")  ~ "Col_Spring",
      pst_agg %in% c("CR-lower_fa", "CR-upper_su/fa") | 
        region1name == "Willamette_R" ~  "Col_Summer_Fall",
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN", "NEVI") ~
        "ECVI_SOMN",
      pst_agg %in% c("CA_ORCST", "WACST", "Russia", "NBC_SEAK", "Yukon") ~
        "other",
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
    stock, stock_prob, stock_group, age_stock_group, stock_id_method,
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
# saveRDS(
#   ppn_dat_size, here::here("data", "rkw_diet", "cleaned_ppn_dat_size.rds")
# )


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


## CHECK PBT COVERAGE ----------------------------------------------------------

# import PBT estimates to estimate coverage 
pbt_rate <- readRDS(
  here::here("data", "sep", "mean_pbt_rate.rds")
) %>% 
  mutate(
    pbt_stock = ifelse(
      collection_extract == "SHUSWAP_RIVER-LOWER (includes Kingfisher)",
      "SHUSWAP_RIVER_LOWER",
      gsub("-", "_", collection_extract) %>% 
        toupper()
    )
  ) %>% 
  select(
    pbt_stock, pbt_brood_year_n = year, tag_rate
  )


dat_pbt <- dat %>% 
  mutate(
    pbt_stock = gsub("-", "_", stock) %>% 
      gsub(" H", "", .) %>% 
      gsub(" ", "_", .) %>% 
      toupper(),
    # replace stocks with mismatched names
    pbt_stock = case_when(
      pbt_stock %in% c("CHILLIWACK_RIVER", "CHILLIWACK_VEDDER_RIVER",
                       "H_CHILLIWACK_RIVER", "S_CHILLIWACK_R") ~ 
        "CHILLIWACK_RIVER_FALL",
      pbt_stock == "COLDWATER_RIVER_UPPER" ~ "COLWDWATER_RIVER",
      pbt_stock == "SHUSWAP_RIVER,_MIDDLE," ~ "SHUSWAP_RIVER_MIDDLE",
      pbt_stock == "CHEHALIS_RIVER" ~ "CHEHALIS_RIVER_SUMMER",
      pbt_stock == "BIG_QUALICUM_RIVER" ~ "QUALICUM_RIVER",
      pbt_stock == "ATNARKO_RIVER_UPPER" ~ "ATNARKO_RIVER",
      pbt_stock == "H_NITINAT_RIVER" ~ "NITINAT_RIVER",
      pbt_stock == "S_BURMAN_R" ~ "BURMAN_RIVER",        
      pbt_stock %in% c("S_CONUMA_R", "H_CONUMA_R") ~ "CONUMA_RIVER",        
      pbt_stock == "S_COWICHAN_R" ~ "COWICHAN_RIVER",
      pbt_stock %in%  c("S_FIRST_LK/GSVI", "S_NANAIMO_R", 
                        "NANAIMO_RIVER_FALL")   ~ "NANAIMO_RIVER_FALL", 
      pbt_stock == "S_GOLD_R" ~ "GOLD_RIVER",
      pbt_stock == "NANAIMO_RIVER_S" ~ "NANAIMO_RIVER_SUMMER",       
      pbt_stock == "S_NITINAT_R" ~ "NITINAT_RIVER",       
      pbt_stock %in% c("S_ROBERTSON_CR", "H_ROBERTSON_CR", 
                       "ROVERTSON_CREEK") ~ "ROBERTSON_CREEK",    
      grepl("TAHSIS", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("TAHSISH", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("PUNTLEDGE", pbt_stock) ~ "PUNTLEDGE_RIVER_FALL",
      pbt_stock == "H_SOOKE_RIVER" ~ "SOOKE_RIVER",
      pbt_stock == "S_SALMON_R/JNST" ~ "SALMON_RIVER_JNST",
      pbt_stock == "S_SARITA_R" ~ "SARITA_RIVER",        
      pbt_stock == "S_SUCWOA/TLUPANA_R" ~ "TLUPANA_RIVER",
      pbt_stock == "CARIBOO_RIVER" ~ "CARIBOO",
      pbt_stock == "S_LEINER_R" ~ "LEINER_RIVER",
      pbt_stock == "S_MARBLE_R" ~ "MARBLE_RIVER",
      pbt_stock == "S_NAHMINT_R" ~ "NAHMINT_RIVER",
      pbt_stock == "S_QUINSAM_R" ~ "QUINSAM_RIVER",
      pbt_stock == "S_SAN_JUAN_R" ~ "SAN_JUAN_RIVER",
      pbt_stock == "WOSS_LAKE"~ "WOSS_RIVER",
      TRUE ~ pbt_stock
    ),
    pbt_brood_year_n = year - as.numeric(as.character(age_f))
  ) %>% 
  left_join(
    ., pbt_rate, by = c("pbt_brood_year_n", "pbt_stock")
  )

dat_pbt %>% 
  mutate(
    high_phos = ifelse(
      grepl("ROBERT", stock) | grepl("NITINAT", stock) | grepl("SARITA", stock) |
        grepl("SOOKE", stock) | grepl("QUALIC", stock) | 
        grepl("NANAIMO", stock) | grepl("PUNTLE", stock) | 
        grepl("QUINS", stock),
      TRUE, FALSE
    )
  ) %>% 
  filter(high_phos == TRUE)
