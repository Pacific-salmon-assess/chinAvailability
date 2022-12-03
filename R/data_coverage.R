## Sample Coverage
# Use GSI samples through 2020 to explore composition coverage and determine
# spatial scale of models; equivalent process with creel data

library(tidyverse)

# samples up to 2019 generated in Chinook distribution project
# rec_raw <- readRDS(here::here("data", "rec", "recIndProbsLong.rds"))

### TEMPORARY AMALGAMATION ONLY
rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date),
         subarea = case_when(
           # subarea %in% c("18B", "18D", "18E") ~ "18BDE",
           subarea %in% c("21B") ~ "21A",
           subarea %in% c("US7") ~ "19B",
           subarea %in% c("29-11", "29I") ~ "29B",
           subarea %in% c("29DW", "29DE") ~ "29D",
           # subarea %in% c("121C") ~ "121BC",
           # consolidate subareas due to small sample sizes and doc issues w/
           # convergence
           # subarea %in% c("19C", "19D", "19E") ~ "19CDE",
           # subarea %in% c("20A", "20E") ~ "20AE",
           # subarea %in% c("20C", "20D") ~ "20CD",
           # subarea %in% c("29G", "29F") ~ "29FG",
           # subarea == "29E" ~ "29DE",
           # grepl("29D", subarea) ~ "29DE",
           # subarea %in% c("29B", "29C") ~ "29BC",
           TRUE ~ subarea
         ),
         reg = case_when(
           subarea == "21A" ~ "SWVI",
           area %in% c("121", "21") ~ "SWVI",
           subarea == "19C" ~ "SSoG",
           area %in% c("20W", "20E", "19JDF") ~ "JdFS",
           area %in% c("18", "19GST", "29") ~ "SSoG",
           TRUE ~ "out"
         ),
         reg = as.factor(reg),
         core_area = case_when(
           subarea == "121B" ~ "no",
           subarea == "29DE" ~ "yes",
           reg %in% c("SWVI") ~ "yes", 
           subarea %in% c("18BDE", "19B", "19C", "20AE") ~ "yes",
           TRUE ~ "no"
         )
  ) 

alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")

# sample coverage for GSI ------------------------------------------------------

subarea_list <- rec_raw %>%
  group_by(area, subarea, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>% 
  mutate(year = as.factor(year)) %>% 
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