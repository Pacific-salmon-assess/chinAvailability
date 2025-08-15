## Diet Explore
# Clean and begin analysis of SRKW diet data
# - temporal sampling coverage bar plots
# - stock comp bar plots
# - size comp bar plots

# Set French language option
FRENCH <- TRUE

# Create appropriate figure directories
if (FRENCH) {
  dir.create("figs-french", showWarnings = FALSE)
  dir.create("figs-french/rkw_diet", showWarnings = FALSE)
  fig_dir <- "figs-french"
} else {
  dir.create("figs/rkw_diet", showWarnings = FALSE)
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

# Translation functions for facet labels and legends
translate_era <- function(era_values) {
  if (FRENCH) {
    case_when(
      era_values == "current" ~ "actuel",
      era_values == "early" ~ "précoce", 
      TRUE ~ era_values
    )
  } else {
    era_values
  }
}

translate_strata <- function(strata_values) {
  if (FRENCH) {
    case_when(
      strata_values == "Juan\nde Fuca" ~ "Juan\nde Fuca",
      strata_values == "S. Gulf\nIslands" ~ "Îles du Golfe\ndu Sud",
      strata_values == "Swiftsure\nBank" ~ "Banc\nSwiftsure",
      strata_values == "Port\nRenfrew" ~ "Port\nRenfrew",
      strata_values == "Sooke/\nVictoria" ~ "Sooke/\nVictoria",
      strata_values == "Saanich" ~ "Saanich",
      strata_values == "San Juan\nIslands" ~ "Îles\nSan Juan",
      TRUE ~ strata_values
    )
  } else {
    strata_values
  }
}

library(tidyverse)
library(mgcv)
library(mvtweedie)
library(sf)

# import cleaned data 
stock_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_diet_samples_stock.rds")
  ) %>% 
  # remove uncertain Spring 4.2 assignments
  filter(
    !(grepl("BESS", stock) | grepl("DUTEA", stock))
  ) 

size_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_diet_samples_size.rds")
  )


## aggregate data (calculate mean location and sample size of each event) 
# by week-year-strata for simulation sampling; size and stock data separately
# due to slightly different sample sizes (not all samples could be aged)
sample_key_stock <- stock_dat %>% 
  select(sample_id, id, utm_y, utm_x) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarize(
    n_samples = length(unique(id)),
    utm_y = mean(utm_y),
    utm_x = mean(utm_x)
  )
ppn_dat_stock <- stock_dat %>% 
  group_by(sample_id, era, year, week_n, month, strata, stock_group) %>% 
  summarize(
    agg_count = sum(stock_prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key_stock, by = "sample_id") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 


sample_key_size <- size_dat %>% 
  select(sample_id, id, utm_y, utm_x) %>% 
  distinct() %>% 
  group_by(sample_id) %>% 
  summarize(
    n_samples = length(unique(id)),
    utm_y = mean(utm_y),
    utm_x = mean(utm_x)
  )
ppn_dat_size <- size_dat %>% 
  group_by(sample_id, era, year, week_n, month, strata, size_bin) %>% 
  summarize(
    agg_count = sum(prob),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  left_join(., sample_key_size, by = "sample_id") %>% 
  mutate(
    agg_prob = agg_count / n_samples
  ) 

# ppn_dat_list <- list(ppn_dat_stock, ppn_dat_size)
# names(ppn_dat_list) <- c("stock", "size")
# 
# saveRDS(
#   ppn_dat_list, here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
# )



## PALETTES --------------------------------------------------------------------
age_pal <- c("grey30", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c")
names(age_pal) <- c(NA, "1", "2", "3", "4")

# SMU colour palette
smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(stock_dat$stock_group)

# era shape palette
era_pal <- c(15, 16)
# Set palette names based on language setting
if (FRENCH) {
  names(era_pal) <- translate_era(levels(stock_dat$era))
} else {
  names(era_pal) <- levels(stock_dat$era)
}

# size colour palette
size_colour_pal <- c("grey30", "#8c510a", "#f6e8c3", "#c7eae5", "#01665e")
names(size_colour_pal) <- c(NA, levels(ppn_dat_size$size_bin))


## RAW DATA FIGURES ------------------------------------------------------------

sample_gap_poly <- data.frame(
  x = c(25, 25, 41, 41),
  y = c(2014, 2016, 2016, 2014)
)

# sample coverage through time and among strata
diet_samp_cov <- stock_dat %>% 
  group_by(week_n, strata, year, era) %>% 
  summarize(n = length(unique(id))) %>% 
  mutate(
    era = translate_era(era),
    strata = translate_strata(strata)
  ) %>%
  ggplot(.) +
  geom_point(aes(x = week_n, y = year, size = n, shape = era), 
             alpha = 0.6) +
  facet_wrap(~strata) +
  geom_polygon(data = sample_gap_poly, aes(x = x, y = y), 
               fill = "red", alpha = 0.3) +
  scale_size_continuous(name = tr("Sample\nSize", "Taille\nd'échantillon")) +
  scale_shape_manual(values = era_pal, name = tr("Sample\nEra", "Période\nd'échantillonnage")) +
  scale_x_continuous(
    breaks = c(25, 29, 33, 37, 41),
    labels = tr(c("Jun", "Jul", "Aug", "Sep", "Oct"), c("juin", "juil", "août", "sep", "oct"))
  ) +
  labs(
    size = tr("Sample Size", "Taille d'échantillon"),
    shape = tr("Sample Era", "Période d'échantillonnage")
  ) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title = element_blank()
  )


# sample size labels
samp_size_stock <- stock_dat %>% 
  group_by(era, strata) %>% 
  summarize(
    n = length(unique(id))
  ) %>%
  mutate(
    era = translate_era(era),
    strata = translate_strata(strata)
  )

samp_size_age <- size_dat %>% 
  group_by(era, strata) %>% 
  summarize(
    n = length(unique(id))
  ) %>%
  mutate(
    era = translate_era(era),
    strata = translate_strata(strata)
  )


# stacked bar plots
diet_samp_bar <- ggplot(stock_dat %>% 
                        mutate(era = translate_era(era),
                               strata = translate_strata(strata))) +
  geom_bar(aes(x = month, y = stock_prob, fill = stock_group), 
           stat = "identity") +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(values = smu_colour_pal, name = tr("Stock", "Stock")) +
  labs(
    y = tr("Prey Remains Composition\n(Individual Samples)", "Composition des restes de proie\n(Échantillons individuels)"),
    fill = tr("Stock", "Stock")
  ) +
  scale_x_continuous(
    breaks = c(6, 7, 8, 9, 10),
    labels = tr(c("Jun", "Jul", "Aug", "Sep", "Oct"), c("juin", "juil", "août", "sep", "oct"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7.5)
  ) +
  geom_text(
    data = samp_size_stock, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

age_samp_bar <- size_dat %>% 
  select(id, strata, era, month, fw_year, sw_year, total_year) %>%
  distinct() %>% 
  mutate(era = translate_era(era),
         strata = translate_strata(strata)) %>%
  ggplot(.) +
  geom_bar(aes(x = month, fill = as.factor(sw_year))) +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(
    name = tr("Marine\nAge", "Âge\nmarin"), values = age_pal, na.value = "grey60" 
    ) +
  labs(
    y = tr("Prey Remains Composition\n(Individual Samples)", "Composition des restes de proie\n(Échantillons individuels)"),
    fill = tr("Marine Age", "Âge marin")
  ) +
  scale_x_continuous(
    breaks = c(6, 7, 8, 9, 10),
    labels = tr(c("Jun", "Jul", "Aug", "Sep", "Oct"), c("juin", "juil", "août", "sep", "oct"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7.5)
  )+
  geom_text(
    data = samp_size_age, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )

# box plots of size
size_samp_bar <- ggplot(size_dat %>% 
                        mutate(era = translate_era(era),
                               strata = translate_strata(strata))) +
  geom_bar(aes(x = month, y = prob, fill = size_bin), 
           stat = "identity") +
  facet_grid(era~strata) +
  ggsidekick::theme_sleek() +
  scale_fill_manual(name = tr("Size\nClass", "Classe\nde taille"), values = size_colour_pal, 
                    na.value = "grey60" ) +
  labs(
    y = tr("Prey Remains Composition\n(Individual Samples)", "Composition des restes de proie\n(Échantillons individuels)"),
    fill = tr("Size Class", "Classe de taille")
  ) +
  scale_x_continuous(
    breaks = c(6, 7, 8, 9, 10),
    labels = tr(c("Jun", "Jul", "Aug", "Sep", "Oct"), c("juin", "juil", "août", "sep", "oct"))
  ) +
  theme(
    legend.position = "top",
    axis.title.x = element_blank(),
    axis.text = element_text(size = 7.5)
  ) +
  geom_text(
    data = samp_size_age, aes(x = Inf, y = Inf, label = paste(n)),
    hjust = 1.1, vjust = 1.1
  )



## export 
png(
  fig_path("rkw_diet/temporal_sample_coverage.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
diet_samp_cov
dev.off()

png(
  fig_path("rkw_diet/monthly_comp_bar.png"),
  height = 5, width = 8, units = "in", res = 250
)
diet_samp_bar
dev.off()

png(
  fig_path("rkw_diet/comp_bar_prey_age.png"),
  height = 5, width = 8, units = "in", res = 250
)
age_samp_bar
dev.off()

png(
  fig_path("rkw_diet/comp_bar_prey_size.png"),
  height = 5, width = 8, units = "in", res = 250
)
size_samp_bar
dev.off()


## CHECK PBT COVERAGE ----------------------------------------------------------

# import PBT estimates to estimate coverage 
pbt_rate <- readRDS(here::here("data", "sep", "cleaned_pbt.rds")) %>%
  rename(
    pbt_brood_year_n = brood_year
  )


dat_pbt <- dat %>% 
  mutate(
    pbt_stock = gsub("-", "_", stock) %>% 
      gsub(" H", "", .) %>% 
      gsub(" ", "_", .) %>% 
      toupper(),
    # replace stocks with mismatched names
    pbt_stock = case_when(
      grepl("CAPILANO", pbt_stock) ~ "CAPILANO_RIVER",
      stock %in% c("CHILLIWACK_RIVER_SUMMER") ~ "CHILLIWACK_RIVER_SUMMER",
      grepl("CHILLIWAC", pbt_stock) ~ "CHILLIWACK_RIVER_FALL",
      grepl("HARRISON", pbt_stock) ~ "HARRISON_RIVER",
      grepl("HORSEFLY", pbt_stock) ~ "HORSEFLY_RIVER",
      grepl("COLDWATER", pbt_stock) ~ "COLWDWATER_RIVER",
      grepl("INDIANPOINT", pbt_stock) ~ "INDIANPOINT_CREEK",
      pbt_stock %in% c("SHUSWAP_RIVER,_MIDDLE,") ~ "SHUSWAP_RIVER_MIDDLE",
      grepl("CHEAKAMUS", pbt_stock) ~ "CHEAKAMUS_RIVER",
      grepl("CHEHALIS_RIVER", pbt_stock) ~ "CHEHALIS_RIVER_SUMMER",
      pbt_stock %in% c("L_QUALICUM") ~ "LITTLE_QUALICUM_RIVER",
      grepl("QUALICUM", pbt_stock) ~ "QUALICUM_RIVER",
      grepl("ATNARKO", pbt_stock) ~ "ATNARKO_RIVER",
      grepl("BURMAN", pbt_stock) ~ "BURMAN_RIVER",
      grepl("MCGREGOR", pbt_stock) ~ "MCGREGOR_RIVER",
      grepl("MORKILL", pbt_stock) ~ "MORKILL_RIVER",
      grepl("PHILLIPS", pbt_stock) ~ "PHILLIPS_RIVER",
      grepl("NAZKO", pbt_stock) ~ "NAZKO_RIVER",
      grepl("CONUMA", pbt_stock) ~ "CONUMA_RIVER",
      grepl("NECHAKO", pbt_stock) ~ "NECHAKO_RIVER",
      grepl("NICOLA", pbt_stock) ~ "NICOLA_RIVER",
      grepl("PORTAGE", pbt_stock) ~ "PORTAGE_CREEK",
      grepl("COWICHAN", pbt_stock) ~ "COWICHAN_RIVER",
      grepl("QUESNEL", pbt_stock) ~ "QUESNEL_RIVER",
      grepl("SERPENTINE", pbt_stock) ~ "SERPENTINE_RIVER",
      grepl("TORPY", pbt_stock) ~ "TORPY_RIVER",
      grepl("SLIM_C", pbt_stock) ~ "SLIM_CREEK",
      grepl("SPIUS", pbt_stock) ~ "SPIUS_CREEK",
      grepl("SWIFT", pbt_stock) ~ "SWIFT_CREEK",
      pbt_stock == "SALMON_RIVER(THOM)" ~ "SALMON_RIVER_SOTH",
      pbt_stock %in% 
        c("S_FIRST_LK/GSVI", "S_NANAIMO_R", "NANAIMO_RIVER_FALL", "NANAIMO",
          "NANAIMO_RIVER_UPPER", "NANAIMO_RIVER_F", "NANAIMO_F", 
          "NANAIMOUPPER") ~ 
        "NANAIMO_RIVER_FALL", 
      pbt_stock %in% c("NANAIMO_RIVER_S", "NANAIMO_SU", "NANAIMOSU_BACKX") ~ 
        "NANAIMO_RIVER_SUMMER",       
      grepl("ASHLU", pbt_stock) ~ "ASHLU_CREEK",
      grepl("NAHMINT", pbt_stock) ~ "NAHMINT_RIVER",
      grepl("THORNTON", pbt_stock) ~ "THORNTON_CREEK",
      grepl("GOLD", pbt_stock) ~ "GOLD_RIVER",
      grepl("SOOKE", pbt_stock) ~ "SOOKE_RIVER",
      grepl("NITINAT", pbt_stock) ~ "NITINAT_RIVER",
      grepl("MAMQUAM", pbt_stock) ~ "MAMQUAM_RIVER",
      grepl("SHOVEL", pbt_stock) ~ "SHOVELNOSE_CREEK",
      grepl("ROBERTSON", pbt_stock) ~ "ROBERTSON_CREEK",
      grepl("TAHSIS", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("TAHSISH", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("CHEMAINUS", pbt_stock) ~ "CHEMAINUS_RIVER",
      grepl("PUNTLEDGE", pbt_stock) ~ "PUNTLEDGE_RIVER_FALL",
      pbt_stock == "S_SALMON_R/JNST" ~ "SALMON_RIVER_JNST",
      grepl("SARITA", pbt_stock) ~ "SARITA_RIVER",
      grepl("TLUPANA", pbt_stock) ~ "TLUPANA_RIVER",
      grepl("COLDWATER", stock) ~ "COLDWATER_RIVER",
      grepl("DEADMAN", stock) ~ "DEADMAN_RIVER",
      pbt_stock == "COTTONWOOD_RIVER_LOWER" ~ "COTTONWOOD_RIVER_(LOWER)",
      grepl("CARIBOO", pbt_stock) ~ "CARIBOO",
      grepl("ENDAKO", pbt_stock) ~ "ENDAKO_RIVER",
      grepl("LEINER", pbt_stock) ~ "LEINER_RIVER",
      grepl("WILLOW", pbt_stock) ~ "WILLOW_RIVER",
      grepl("MARBLE", pbt_stock) ~ "MARBLE_RIVER",
      grepl("QUINSAM", pbt_stock) ~ "QUINSAM_RIVER",
      grepl("WOSS", pbt_stock) ~ "WOSS_RIVER",
      grepl("KENNEDY", pbt_stock) ~ "KENNEDY_RIVER_LOWER",
      grepl("NIMPKISH", pbt_stock) ~ "NIMPKISH_RIVER",
      grepl("TRANQUIL", pbt_stock) ~ "TRANQUIL_CREEK",
      grepl("TOQUART", pbt_stock) ~ "TOQUART_RIVER",
      grepl("BEDWELL", pbt_stock) ~ "BEDWELL_RIVER",
      grepl("SAN_JUAN", pbt_stock) ~ "SAN_JUAN_RIVER",
      grepl("BOWRON", pbt_stock) ~ "BOWRON_RIVER",
      grepl("BRIDGE", pbt_stock) ~ "BRIDGE_RIVER",
      stock %in% c("CARIBOO_RIVER-UPPER", "U_CARIBOO") ~ 
        "CARIBOO_RIVER_UPPER",
      grepl("CHILAKO", pbt_stock) ~ "CHILAKO_RIVER",
      grepl("CHILCOTIN", pbt_stock) ~ "CHILCOTIN_RIVER_UPPER",
      grepl("CHILKO", pbt_stock) ~ "CHILKO_RIVER",
      TRUE ~ pbt_stock
    ),
    pbt_brood_year_n = year - as.numeric(as.character(age_f))
  ) %>% 
  left_join(
    ., pbt_rate, by = c("pbt_brood_year_n", "pbt_stock")
  )

dat_pbt %>% 
  filter(high == TRUE) %>% 
  group_by(id) %>% 
  summarize(sum(stock_prob))

dat_pbt %>% 
  mutate(
    high_phos = ifelse(
      grepl("ROBERT", stock) | grepl("NITINAT", stock) | grepl("SARITA", stock) |
        grepl("SOOKE", stock) | grepl("QUALIC", stock) | 
        grepl("NANAIMO", stock) | grepl("PUNTLE", stock) | 
        grepl("QUINS", stock)  | grepl("JUAN", stock),
      TRUE, FALSE
    )
  ) %>% 
  filter(high_phos == TRUE) %>% 
  group_by(id, stock, stock_id_method, pbt_brood_year_n) %>% 
  summarize(sum_prob = sum(stock_prob)) %>% 
  print(n = Inf)

