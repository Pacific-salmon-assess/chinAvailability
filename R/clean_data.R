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

# recreational composition data
rec <- readRDS(here::here("data", "rec", 
                          "recIndProbsLong.rds")) %>% 
  filter(legal == "legal",
         !region %in% c("Queen Charlotte Sound")) %>%
  mutate(
    temp_strata = paste(month, region, sep = "_"),
    sample_id = paste(temp_strata, jDay, year, sep = "_"),
    # min_m = case_when(
    #   region %in% c("N. Strait of Georgia", "S. Strait of Georgia") ~ 5,
    #   region %in% c("Juan de Fuca Strait",
    #                 "Queen Charlotte and\nJohnstone Straits") ~ 6
    # ),
    # max_m = case_when(
    #   region %in% c("Juan de Fuca Strait", "N. Strait of Georgia",
    #                 "S. Strait of Georgia") ~ 9,
    #   region == "Queen Charlotte and\nJohnstone Straits" ~ 8
    # ),
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

