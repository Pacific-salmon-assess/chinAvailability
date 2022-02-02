## Data Clean
# For now utilize datasets from chinook distribution manuscript
# Dec. 1, 2021

library(tidyverse)

# COMPOSITION HELPER FUNCTIONS -------------------------------------------------

# helper function to calculate aggregate probs
calc_agg_prob <- function(grouped_data, full_data) {
  grouped_data %>% 
    summarize(agg_prob = sum(adj_prob), .groups = 'drop') %>% 
    arrange(sample_id, desc(agg_prob)) %>%
    distinct() %>% 
    left_join(.,
              full_data %>% 
                #important to subset appropriately to remove individual traits
                #that lead to duplicates
                select(sample_id, region, year, month, gear, month_n) %>% 
                distinct(),
              by = "sample_id")
}

# helper function to pool non-focal stocks
pool_aggs <- function(full_data) {
  full_data %>% 
    mutate(
      reg1 = case_when(
        Region1Name %in% c("Fraser_Fall", "ECVI", "Fraser_Summer_4.1", 
                           "Fraser_Spring_4.2", "Fraser_Spring_5.2", 
                           "Fraser_Summer_5.2", "WCVI", "SOMN") ~ Region1Name,
        TRUE ~ "Other"
      )
    )
}

clean_comp <- function(grouping_col, raw_data, ...) {
  temp <- raw_data %>% 
    pool_aggs() %>% 
    rename(agg = grouping_col) %>%
    group_by(sample_id, agg) %>%
    calc_agg_prob(., raw_data) %>% 
    group_by(region, month, year) %>%
    mutate(nn = sum(agg_prob)) %>%
    #remove strata with less than 10 individuals total
    filter(!nn < 10) %>%
    ungroup() %>%
    droplevels() %>%
    mutate(region_c = as.character(region),
           region = factor(abbreviate(region, minlength = 4))) %>% 
    select(sample_id, gear, region, region_c, year, month, month_n, agg,
           agg_prob, nn) %>%
    distinct()
  if (raw_data$gear[1] == "sport") {
    temp %>%
      mutate(region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG"))
  } else {
    temp
  }
}

# DATA CLEAN -------------------------------------------------------------------

# recreational composition data through 2019
rec_raw <- readRDS(here::here("data", "rec", 
                          "recIndProbsLong.rds")) 
rec <- rec_raw %>% 
  filter(legal == "legal",
         !region %in% c("Queen Charlotte Sound")) %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    min_m = case_when(
      region %in% c("N. Strait of Georgia", "S. Strait of Georgia") ~ 1,
      region == "Queen Charlotte and\nJohnstone Straits" ~ 6,
      region == c("Juan de Fuca Strait") ~ 3
    ),
    max_m = case_when(
      region %in% c("Juan de Fuca Strait", 
                    "N. Strait of Georgia") ~ 9,
      region == "S. Strait of Georgia" ~ 12,
      region == "Queen Charlotte and\nJohnstone Straits" ~ 8
    ),
    coarse_agg = case_when(
      pst_agg %in% c("WCVI", "NBC_SEAK") ~ "BC-coastal",
      pst_agg %in% c("CA_ORCST", "WACST") ~ "US-coastal",
      pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa") ~ "CR-fall",
      pst_agg %in% c("CR-lower_sp", "CR-upper_sp") ~ "CR-spring",
      TRUE ~ pst_agg
    )
  ) %>% 
  group_by(region) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>%
  ungroup() %>% 
  droplevels()

clean_rec <- rec %>% 
  clean_comp(., grouping_col = "coarse_agg")



# recreational composition data since through 2021 (clean to match rec_raw)
rec_raw_new <- read.csv(here::here("data", "rec", "sc_biodata_jul8_21.csv"),
                        stringsAsFactors = FALSE, na.strings=c("","NA")) %>% 
  janitor::clean_names(.) 

tt <- rec_raw_new %>% 
  # change US area 7 (near San Juan island) to SSoG
  mutate(
    area = ifelse(area == "US7", "19GST", area),
    # region = case_when(
    area_n = as.numeric(str_replace_all(area, "[:letter:]", ""))) %>% 
  glimpse()

,
    #   # separate northern areas of 13 (normally in JS) and add to NSoG
    #   subarea %in% c("13M", "13N") ~ "N. Strait of Georgia",
    #   area > 124 ~ "NWVI",
    #   area < 28 & area > 24 ~ "NWVI",
    #   area %in% c("20", "121", "21") ~ "Juan de Fuca Strait",
    #   is.na(area) ~ "Juan de Fuca Strait",
    #   area < 125 & area > 120 ~ "SWVI",
    #   area < 25 & area > 20 ~ "SWVI",
    #   area %in% c("14", "15", "16") ~ "N. Strait of Georgia",
    #   area %in% c("17", "18", "19", "28", "29") ~ "S. Strait of Georgia",
    #   area %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
    #   area %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
    # ),
    legal_lim = case_when(
      area < 20 & area > 11 ~ 620,
      area %in% c("28, 29") ~ 620,
      TRUE ~ 450
    ),
    legal = case_when(
      length_mm >= legal_lim ~ "legal",
      length_mm < legal_lim ~ "sublegal",
      disposition == "Kept" ~ "legal",
      KEPTREL == "Rel" ~ "sublegal"
    ),
  ) %>% 
  select(
    id = biokey, region 
  )
  


# stock key to 
stockKey <- readRDS(here::here("data", "rec", "finalStockList_Nov2020.rds"))






# recreational catch data
rec_catch <- readRDS(here::here("data", "rec", "month_area_recCatch.rds")) %>% 
  # identify months to exclude based on minimal catch or comp estimates 
  #(make region specific)
  mutate(
    min_m = case_when(
      region %in% c("NSoG", "SSoG") ~ 5,
      region %in% c("JdFS", "QCaJS") ~ 6
    ),
    max_m = case_when(
      region == "QCaJS" ~ 8, 
      region %in% c("JdFS", "NSoG", "SSoG")  ~ 9
    ),
    region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG")
  ) %>% 
  # drop areas with fewer than 10 datapoints (month-year-area observation = 1)
  group_by(area_n) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>%
  droplevels()


saveRDS(clean_rec, here::here("data", "rec", "coarse_rec_comp.rds"))
saveRDS(rec_catch, here::here("data", "rec", "month_area_recCatch_clean.rds"))


## EXPLORE STAT AREA COVERAGE --------------------------------------------------

stat_area_samples <- rec_raw %>% 
  filter(region == c("Juan de Fuca Strait", "S. Strait of Georgia")) %>% 
  group_by(region, area, month, year) %>% 
  summarize(n_samples = length(unique(id)))

ggplot(stat_area_samples) +
  geom_boxplot(aes(x = month, y = n_samples, fill = region)) +
  facet_wrap(~area)

dum <- stat_area_samples %>% 
  filter(area %in% c("21", "121", "19"),
         month %in% c("5", "6", "7", "8", "9"))

sum(dum$n_samples) / (6*5*3)

rec_catch %>% 
  filter(area %in% c("21", "121", "19")) %>% 
  group_by(area, month) %>% 
  tally()


## EXPLORE SIZE DATA -----------------------------------------------------------

size_dat <- rec_raw %>% 
  select(id, fl, jDay, month_n, year, area, region, temp_strata) %>% 
  distinct() %>% 
  filter(!is.na(fl)) %>% 
  mutate(
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 450, 600, 750, 900, Inf), 
      labels=c("<45", "45-60", "60-75", "75-90", ">90")
    )
  ) 


# visualize sample coverage through space and time
size_n <- size_dat %>% 
  group_by(month_n, year, region, temp_strata) %>% 
  tally() 

ggplot(size_n) +
  geom_raster(aes(x = month_n, y = year, fill = n)) +
  facet_wrap(~region)


# visualize changes in size composition
size_dat %>%
  group_by(size_bin, year, region) %>%
  summarize(n = length(unique(id))) %>%
  group_by(year, region) %>% 
  mutate(total_n = sum(n),
            ppn_obs = n / total_n) %>% 
  ggplot(., aes(x = as.factor(year), y = ppn_obs, fill = size_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_wrap(~region)


# identify months/years with sufficient sample sizes
size_n %>% 
  group_by(month_n, region, temp_strata) %>% 
  summarize(max(n)) %>% 
  arrange(region) %>% 
  print(n = Inf)


# export size data
rec_size_out <- size_dat %>%
  mutate(
    min_m = ifelse(region == "Queen Charlotte and\nJohnstone Straits" , 5, 1),
    max_m = ifelse(region == "Queen Charlotte and\nJohnstone Straits" , 9, 12),
  ) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>%
  select(-min_m, -max_m) 


# export size data as proportions
rec_size_ppn_out <- size_dat %>%
  mutate(sample_id = paste(month_n, region, jDay, year, sep = "_"),
         region_c = region,
         region = factor(abbreviate(region, minlength = 4)),
         region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG")) %>%  
  group_by(sample_id, region, region_c, year, month_n, size_bin) %>% 
  summarize(sum_count = length(unique(id)), .groups = "drop") 


saveRDS(rec_size_out, here::here("data", "rec", "rec_size.rds"))
saveRDS(rec_size_ppn_out, here::here("data", "rec", "rec_size_ppn.rds"))
