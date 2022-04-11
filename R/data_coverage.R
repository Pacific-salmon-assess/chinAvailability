## Sample Coverage
# Use GSI samples through 2020 to explore composition coverage and determine
# spatial scale of models; equivalent process with creel data

library(tidyverse)

# samples up to 2019 generated in Chinook distribution project
# rec_raw <- readRDS(here::here("data", "rec", "recIndProbsLong.rds"))

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date)) 


# sample coverage for GSI ------------------------------------------------------

subarea_list <- rec_raw %>%
  group_by(area, subarea, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>% 
  mutate(year = as.factor(year),
         core_area = case_when(
           subarea %in% c("121A", "121B", "21A") ~ "yes",
           area %in% "19GST" ~ "yes",
           TRUE ~ "no"
         )) %>% 
  ungroup() %>% 
  split(., .$region)
area_list <- rec_raw %>%
  group_by(area, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>%
  ungroup() %>% 
  split(., .$region)


# visualize coverage by subarea or area
subarea_coverage <- purrr::map2(
  subarea_list, names(subarea_list),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = year)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_wrap(~subarea, scales = "free_y") +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)

area_coverage <- purrr::map2(
  area_list, names(area_list),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = legal)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_grid(area~year) +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)


# as above but focused only on core regions
subarea_trim <- subarea_list %>% 
  bind_rows() %>% 
  filter(area %in% c("19GST", "19JDF", "20E", "20W", "21", "121")) 

alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")

png(here::here("figs", "data_coverage", "gsi_samples_subarea_trim.png"))
ggplot(subarea_trim, 
       aes(x = as.factor(month_n), y = n, fill = year, alpha = core_area)) +
  geom_bar(position="stack", stat="identity") +
  scale_alpha_manual(values = alpha_scale) +
  ggsidekick::theme_sleek() +
  facet_wrap(~subarea, scales = "free_y") +
  labs(x = "Month", y = "Samples")
dev.off()




# export
pdf(here::here("figs", "data_coverage", "gsi_samples_new_subarea.pdf"))
subarea_coverage
dev.off()

pdf(here::here("figs", "data_coverage", "gsi_samples_new_area.pdf"))
area_coverage
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

# subset to focus on areas of interest
subset_stock_comp <- stock_comp_dat %>% 
  filter(
    area %in% c("121", "21", "20E", "20W", "19JDF", "19GST", "18", "123")
  ) %>% 
  mutate(
    reg = abbreviate(region, minlength = 4),
    area = paste(reg, area, sep = "_")
  ) %>% 
  ggplot(., aes(fill = pst_agg, y = agg_ppn, x = month_n)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_manual(name = "Stock", values = stock_pal) +
  facet_wrap(~area) +
  labs(x = "Month", y = "Agg Probability", title = "Core Study Area") +
  ggsidekick::theme_sleek()


pdf(here::here("figs", "data_coverage", "mean_monthly_stock_comp.pdf"))
seasonal_stock_comp
subset_stock_comp
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