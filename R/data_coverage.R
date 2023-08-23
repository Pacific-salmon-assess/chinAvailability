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
  filter(!is.na(lat), !is.na(lon)) 


# map data
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48, xmax = -122, ymax = 48.8)

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

shape_pal <- c(21, 22, 23)
names(shape_pal) <- unique(obs_stations$cap_region)

pdf(here::here("figs", "data_coverage", "harvest_locations.pdf"), width = 9,
    height = 4)
base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, fill = rkw_habitat, shape = cap_region),
    alpha = 0.4
  ) +
  scale_shape_manual(values = shape_pal) +
  theme(
    legend.position = "top"
  ) +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21)
    )
  )
dev.off()


# sample coverage for GSI ------------------------------------------------------

# number of samples stratified by location and whether spatio-temporal overlap
# w/ RKW samples
pdf(here::here("figs", "data_coverage", "samples_by_location.pdf"))
# stacked bar
rec_raw %>%
  filter(!legal == "sublegal", 
         !cap_region == "outside") %>% 
  group_by(whale_samples_time, rkw_habitat, cap_region, year, month_n) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>% 
  ungroup() %>% 
  ggplot(., aes(x = as.factor(month_n), y = n, fill = whale_samples_time)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_grid(rkw_habitat ~ cap_region, scales = "free_y") +
  labs(y = "Number of Samples")

# bubbles
rec_raw %>%
  filter(!legal == "sublegal", 
         !cap_region == "outside") %>% 
  group_by(cap_region, month_n, rkw_habitat) %>% 
  mutate(total_samples = sum(prob),
         month = as.factor(month_n),
         sample_id = paste(month, cap_region, rkw_habitat, year, sep = "_")) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id)) %>% as.numeric) %>% 
  ungroup() %>% 
  select(sample_id, month_n, cap_region, rkw_habitat, year, nn) %>% 
  distinct() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = rkw_habitat),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~cap_region) +
  ggsidekick::theme_sleek() +
  labs(x = "Month", y = "Year")
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
  facet_grid(rkw_habitat ~ cap_region) +
  labs(x = "Month", y = "Agg Probability") +
  ggsidekick::theme_sleek()
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