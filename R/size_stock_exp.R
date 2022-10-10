## Size by stock 
# Oct. 9, 2022
# Explore patterns in size by areas for different stocks assuming stock
# correct stock ID

library(tidyverse)


comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!is.na(fl)) %>% 
  mutate(
    month = lubridate::month(date),
    mid_reg = case_when(
      pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                     "CR-lower_sp", "CR-upper_sp", "PSD", 
                     "NBC_SEAK", "WACST") ~ pst_agg,
      TRUE ~ Region1Name
    ),
    fl = as.numeric(fl),
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 450, 600, 750, 900, Inf), 
      labels=c("<45", "45-60", "60-75", "75-90", ">90")
    )
  ) %>% 
  #assign fixed probabilities
  group_by(id, year, month, cap_region, area, subarea, fl, size_bin, mid_reg) %>% 
  summarize(
    sum_prob = sum(prob),
    .groups = "drop"
  ) %>% 
  filter(!sum_prob < 0.75) %>% 
  droplevels()

month_ppn <- comp1 %>%
  group_by(size_bin, month, mid_reg, cap_region) %>%
  summarize(n = length(unique(id))) %>%
  group_by(month, mid_reg, cap_region) %>% 
  mutate(total_n = sum(n),
         ppn_obs = n / total_n) %>% 
  ggplot(., aes(x = as.factor(month), y = ppn_obs, fill = size_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_grid(mid_reg~cap_region)

pdf(here::here("figs", "size_exp", "size_ppn_region.pdf"), height = 9, width = 9)
month_ppn
dev.off()


# subset to May-Sep in SRKW regions
trim_comp <- comp1 %>% 
  filter(!subarea %in% c("19A", "18A"),
         area %in% c("18", "19", "20", "21", "121")) %>% 
  mutate(subarea_original = subarea,
       subarea = case_when(
         subarea %in% c("18B", "18D", "18E") ~ "18BDE",
         # consolidate subareas due to small sample sizes and doc issues w/
         # convergence
         subarea %in% c("19C", "19D", "19E") ~ "19CDE",
         subarea %in% c("20A", "20E") ~ "20AE",
         subarea %in% c("20C", "20D") ~ "20CD",
         TRUE ~ subarea
       ),
       reg = case_when(
         # subarea == "123I" ~ "SWVI", # qu
         #correction for subarea 21A IDd as SSoG
         subarea == "21A" ~ "SWVI",
         area %in% c("121", "21") ~ "SWVI",
         subarea == "19C" ~ "SSoG",
         area == "20" ~ "JdFS",
         subarea %in% c("19D", "19E") ~ "JdFS",
         area %in% c("18", "19") ~ "SSoG",
         # area  ~ "SSoG",
         TRUE ~ "out"
       ),
       reg = as.factor(reg),
       core_area = case_when(
         subarea == "121B" ~ "no",
         reg %in% c("SWVI") ~ "yes", 
         subarea %in% c("18BDE", "19B", "19C", "20AE") ~ "yes",
         TRUE ~ "no"
       )) %>% 
  filter(!reg == "out",
         month > 4 & month < 10) %>% 
  droplevels()


month_ppn_core_area <- trim_comp %>%
  group_by(size_bin, month, mid_reg, reg) %>%
  summarize(n = length(unique(id))) %>%
  group_by(month, mid_reg, reg) %>% 
  mutate(total_n = sum(n),
         ppn_obs = n / total_n) %>% 
  ggplot(., aes(x = as.factor(month), y = ppn_obs, fill = size_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_grid(mid_reg~reg)


pdf(here::here("figs", "size_exp", "size_ppn_core_area.pdf"), 
    height = 9, width = 9)
month_ppn_core_area
dev.off()
