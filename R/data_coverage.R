## Sample Coverage
# Use GSI samples through 2020 to explore composition coverage and determine
# spatial scale of models; equivalent process with creel data

library(tidyverse)
library(maptools)
library(rmapshaper)
library(mapdata)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  # redefine region based on analysis
  mutate(
    cap_region = case_when(
      lat < 48.8 & lon > -125.25 & lon < -124.25 ~ "swiftsure",
      lat < 48.45 & lon < -123.4 & lon > -124.25 ~ "sooke",
      TRUE ~ "outside"
    )
  )


# map data
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48, xmax = -122, ymax = 48.8)

## convert back to WGS84 and export
hab_sf <- readRDS(
   here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
)


# spatial distribution ---------------------------------------------------------

# all observed locations
obs_stations <- rec_raw %>%
  filter(lat < 48.8) %>% 
  select(id, lat, lon, rkw_habitat, cap_region) %>% 
  distinct() %>% 
  group_by(lat, lon, rkw_habitat, cap_region) %>% 
  tally()


# map including observed locations
base_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = NA) +
  geom_sf(data = hab_sf, color = "red") +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

pdf(here::here("figs", "data_coverage", "harvest_locations.pdf"))
base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, colour = rkw_habitat, shape = cap_region),
    alpha = 0.4
  ) +
  theme(
    legend.position = "top"
  )
dev.off()


# sample coverage for GSI ------------------------------------------------------

# number of samples stratified by location and whether spatio-temporal overlap
# w/ RKW samples
pdf(here::here("figs", "data_coverage", "samples_by_location.pdf"))
rec_raw %>%
  filter(!legal == "sublegal", 
         !cap_region == "outside") %>% 
  group_by(rkw_habitat, cap_region, year, month_n) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>% 
  mutate(
    whale_samples_time = ifelse(
      (year < 2011 | year > 2017) & month_n %in% c("6", "7", "8") & 
        rkw_habitat == "yes",
      "yes",
      "no"
      )
  ) %>% 
  ungroup() %>% 
  ggplot(., aes(x = as.factor(month_n), y = n, fill = whale_samples_time)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_grid(rkw_habitat ~ cap_region, scales = "free_y") +
  labs(y = "Number of Samples")
dev.off()


# stock comp stratified by location and whether spatio-temporal overlap
# w/ RKW samples
pdf(here::here("figs", "data_coverage", "comp_by_location.pdf"))
rec_raw %>%
  filter(!legal == "sublegal", 
         !cap_region == "outside") %>% 
  group_by(cap_region, month_n, rkw_habitat) %>% 
  mutate(total_samples = sum(prob),
         month = as.factor(month_n)) %>% 
  group_by(cap_region, rkw_habitat, month, total_samples, 
           stock_group) %>% 
  summarize(
    agg_prob = sum(prob),
    agg_ppn = agg_prob / total_samples,
    .groups = "drop"
  ) %>%  
  distinct() %>% 
  ggplot(., aes(fill = stock_group, y = agg_ppn, x = month)) + 
  geom_bar(position="stack", stat="identity") +
  # scale_fill_manual(name = "Stock", values = stock_pal) +
  facet_grid(rkw_habitat ~ cap_region) +
  labs(x = "Month", y = "Agg Probability") +
  ggsidekick::theme_sleek()
dev.off()


# visualize seasonal changes in composition by area (NOTE excludes sublegal
# samples)
stock_comp_dat <- rec_raw %>%
  filter(legal == "legal") %>% 
  group_by(area, month) %>% 
  mutate(total_samples = sum(prob)) %>% 
  group_by(area, region, month, month_n, total_samples, pst_agg) %>% 
  summarize(
    agg_prob = sum(prob),
    agg_ppn = agg_prob / total_samples,
    .groups = "drop"
  ) %>%  
  distinct() 
stock_comp_list <- split(stock_comp_dat, stock_comp_dat$region)

  
stock_pal <- disco::disco("rainbow", n = length(unique(stock_comp_dat$pst_agg)))
  
seasonal_stock_comp <- purrr::map2(
  stock_comp_list, names(stock_comp_list),
  function (x, y) {
    p <- ggplot(x, aes(fill = pst_agg, y = agg_ppn, x = month_n)) + 
      geom_bar(position="stack", stat="identity") +
      scale_fill_manual(name = "Stock", values = stock_pal) +
      facet_wrap(~area) +
      labs(x = "Month", y = "Agg Probability", 
           title = y) +
      ggsidekick::theme_sleek()
    print(p)
  }
)

# as above but for subareas, restricted to regions of focus and aggregated 
# post-hoc

subarea_comp_dat <- rec_raw %>%
  filter(legal == "legal",
         !reg_f == "out") %>%
  group_by(subarea, month) %>% 
  mutate(total_samples = sum(prob)) %>% 
  group_by(subarea, reg_f, month, month_n, total_samples, pst_agg, core_area) %>% 
  summarize(
    agg_prob = sum(prob),
    agg_ppn = agg_prob / total_samples,
    .groups = "drop"
  ) %>%  
  distinct() %>% 
  mutate(strata = paste(reg_f, subarea))

stock_pal <- disco::disco("rainbow", n = length(unique(subarea_comp_dat$pst_agg)))


labs <- subarea_comp_dat %>% 
  dplyr::select(month_n, strata, total_samples) %>% 
  distinct() 


subset_stock_comp <- ggplot() + 
  geom_bar(data = subarea_comp_dat, 
           aes(fill = pst_agg, y = agg_ppn, x = month_n, alpha = core_area),
           position="stack", stat="identity") +
  scale_fill_discrete(name = "Stock") +
  scale_alpha_manual(name = "Core Area", values = alpha_scale) +
  geom_text(data = labs, aes(x = month_n, y = -Inf, label = total_samples),
            position = position_dodge(width = 0.9), size = 2.5,
            vjust = -0.5) +
  facet_wrap(~strata) +
  labs(x = "Month", y = "Agg Probability") +
  ggsidekick::theme_sleek()


pdf(here::here("figs", "data_coverage", "mean_monthly_stock_comp.pdf"))
seasonal_stock_comp
dev.off()

png(here::here("figs", "data_coverage", "mean_monthly_stock_comp_subarea.png"),
    height = 7, width = 10, units = "in", res = 250)
subset_stock_comp
dev.off()


## bubble plot
# alpha_scale <- c(0.3, 0.95)
# names(alpha_scale) <- c("no", "yes")
png(here::here("figs", "data_coverage", "weekly_bubble_subarea.png"),
    height = 7, width = 7, units = "in", res = 250)
rec_raw %>%
  filter(!reg == "out") %>% 
  mutate(
    week = lubridate::week(date),
    month_n = lubridate::month(date),
    sample_id = paste(month_n, subarea, week, year, sep = "_")
    ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id)) %>% as.numeric) %>% 
  ungroup() %>% 
  select(sample_id, year, reg, subarea, month_n, nn) %>%
  distinct() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(subarea, as.numeric(reg)),
             ncol = 3) +
  ggsidekick::theme_sleek() +
  theme(
    legend.position = "null"
  ) +
  labs(x = "Month", y = "Year")
dev.off()


# sample coverage for size -----------------------------------------------------

area_list_fl <- rec_raw %>%
  filter(!is.na(fl)) %>% 
  group_by(area, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>%
  ungroup() %>% 
  split(., .$region)

area_coverage_fl <- purrr::map2(
  area_list_fl, names(area_list_fl),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = legal)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_grid(area~year) +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)

pdf(here::here("figs", "data_coverage", "fl_samples_new_area.pdf"))
area_coverage_fl
dev.off()


# coverage for catch/effort ----------------------------------------------------

catch <- readRDS(here::here("data", "rec", "rec_creel_subarea.rds"))

catch_dist_list <- catch %>%
  group_by(region, area_n, subarea, month_n, year, effort) %>%
  summarize(
    sum_catch = sum(catch),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  mutate(
    samples = case_when(
      is.na(sum_catch) & is.na(effort) ~ "none",
      is.na(sum_catch) ~ "effort",
      is.na(effort) ~ "catch",
      TRUE ~ "catch & effort"
    )
  ) %>% 
  split(., .$region) %>% 
  purrr::map2(., names(.), function (x, y) {
    ggplot(x) +
      geom_tile(aes(x = month_n, y = year, fill = samples)) +
      facet_wrap(~subarea) +
      ggsidekick::theme_sleek() +
      labs(x = "Year", y = "Samples", title = y)
  })


pdf(here::here("figs", "data_coverage", "creel_subarea.pdf"))
catch_dist_list
dev.off()


catch2 <- readRDS(here::here("data", "rec", "rec_creel_area.rds"))

catch_dist_list <- catch2 %>%
  group_by(region, area_n, month_n, year, effort) %>%
  summarize(
    sum_catch = sum(catch),
    .groups = "drop"
  ) %>% 
  ungroup() %>% 
  mutate(
    samples = case_when(
      is.na(sum_catch) & is.na(effort) ~ "none",
      is.na(sum_catch) ~ "effort",
      is.na(effort) ~ "catch",
      TRUE ~ "catch & effort"
    )
  ) %>% 
  split(., .$region) %>% 
  purrr::map2(., names(.), function (x, y) {
    ggplot(x) +
      geom_tile(aes(x = month_n, y = year, fill = samples)) +
      facet_wrap(~area_n) +
      ggsidekick::theme_sleek() +
      labs(x = "Year", y = "Samples", title = y)
  })


pdf(here::here("figs", "data_coverage", "creel_area.pdf"))
catch_dist_list
dev.off()