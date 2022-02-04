## Data Clean
# For now utilize datasets from chinook distribution manuscript
# Dec. 1, 2021

library(tidyverse)

# stock key to 
stock_key <- readRDS(here::here("data", "rec", "finalStockList_Jan2022.rds"))


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


# INDIVIDUAL DATA CLEAN --------------------------------------------------------

# recreational composition data since through 2021 (clean to match rec_raw)
rec_raw_new <- read.csv(here::here("data", "rec", "sc_biodata_jul8_21.csv"),
                        stringsAsFactors = FALSE, na.strings=c("","NA")) %>% 
  janitor::clean_names(.) 

wide_rec <- rec_raw_new %>% 
  # change US area 7 (near San Juan island) to SSoG
  mutate(
    area = ifelse(area == "US7", "19GST", area),
    # region = case_when(
    area_n = as.numeric(str_replace_all(area, "[:letter:]", "")),
    # separate northern areas of 13 (normally in JS) and add to NSoG
    cap_region = case_when(
      subarea %in% c("13M", "13N") ~ "N. Strait of Georgia",
      area_n > 124 ~ "NWVI",
      area_n < 28 & area > 24 ~ "NWVI",
      area %in% c("20W", "20E", "20", "121", "21", "19JDF") ~ 
        "Juan de Fuca Strait",
      area_n < 125 & area_n > 120 ~ "SWVI",
      area_n < 25 & area_n > 20 ~ "SWVI",
      area %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      area %in% c("17", "18", "19", "19GST", "28", "29") ~ 
        "S. Strait of Georgia",
      area %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      area %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
    ),
    legal_lim = case_when(
      area_n < 20 & area_n > 11 ~ 620,
      area_n %in% c("28, 29") ~ 620,
      TRUE ~ 450
    ),
    legal = case_when(
      length_mm >= legal_lim ~ "legal",
      length_mm < legal_lim ~ "sublegal",
      disposition == "Kept" ~ "legal",
      disposition == "Released" ~ "sublegal"
    ),
    date = as.POSIXct(collection_date, format="%d-%b-%Y"),
    month_n = lubridate::month(date)
    )


# correct some size entries
# weird_sizes <- wide_rec %>% 
#   filter(length_mm < 150 | length_mm > 1500) %>%
#   select(biokey, length_mm, new_disposition, contains("size")) 
# write.csv(weird_sizes, here::here("data", "rec", "southcoast_size_errors.csv"),
#           row.names = FALSE)

corrected_sizes <- read.csv(
  here::here("data", "rec", "southcoast_size_errors_corrected.csv")
) %>% 
  mutate(
    new_length_mm = ifelse(is.na(new_length_mm), "remove", new_length_mm)
  ) %>% 
  select(biokey, new_length_mm)

wide_rec2 <- full_join(wide_rec, corrected_sizes, by = "biokey") %>%
  mutate(
    fl = ifelse(is.na(new_length_mm), length_mm, new_length_mm) 
  ) 


# GSI CLEAN --------------------------------------------------------------------

# trim for GSI purposes
wide_rec2_trim <- wide_rec2 %>% 
  filter(
    !is.na(resolved_stock_source),
    # for now remove all samples that don't also DNA data (unrepresentative 
    # sampling)
    !is.na(dna_results_stock_1)
  ) %>% 
  mutate(
    # make oto/dna/cwt equivalent
    stock_1 = case_when(
      resolved_stock_source == "DNA" ~ dna_results_stock_1,
      resolved_stock_source == "CWT" ~ cwt_result,
      resolved_stock_source == "Otolith Stock" ~ oto_stock,
    ),
    prob_1 = case_when(
      resolved_stock_source == "DNA" ~ prob_1,
      resolved_stock_source == "CWT" ~ 1.0,
      resolved_stock_source == "Otolith Stock" ~ 1.0
    )
  ) %>% 
  select(
    id = biokey, date, month, year, cap_region, area, area_n, subarea, legal, 
    fl = length_mm, sex, ad = adipose_fin_clipped, resolved_stock_source, 
    stock_1, stock_2 = dna_stock_2, stock_3 = dna_stock_3,
    stock_4 = dna_stock_4, stock_5 = dna_stock_5,
    starts_with("prob"),
    # retain region data to help build stock key (will eventually be removed)
    ends_with("rollup")
  ) 

# replace 0 probabilities with v. small values (just to be safe); then recalc
# ppns
prbs <- wide_rec2_trim %>% 
  select(starts_with("prob")) %>% 
  as.matrix()
prbs[prbs == 0] <- .00001
row_sums <- apply(prbs, 1, sum, na.rm = T)
new_prbs <- prbs / row_sums


#pivot to long (probs and stock IDs separately) and join 
probs <- wide_rec2_trim %>% 
  # replace with updated probabilities from above
  select(-starts_with("prob")) %>% 
  cbind(., new_prbs) %>% 
  pivot_longer(., cols = starts_with("prob"), names_to = "rank", 
               names_pattern = "prob_(.+)",
               values_to = "prob") %>%
  select(id, rank, prob)

regions <- wide_rec2_trim %>% 
  select(id, starts_with("region")) %>% 
  pivot_longer(., cols = starts_with("region"), names_to = "rank", 
               names_pattern = "region_(.+)_rollup",
               values_to = "region") %>%
  select(id, rank, region)

long_rec <- wide_rec2_trim %>% 
  select(-starts_with("prob"), -starts_with("region")) %>% 
  pivot_longer(., cols = starts_with("stock"), names_to = "rank", 
               names_pattern = "stock_(.+)",
               values_to = "stock") %>% 
  left_join(., probs, by = c("id", "rank")) %>% 
  left_join(., regions, by = c("id", "rank")) %>% 
  arrange(desc(date), id, desc(prob)) %>% 
  select(-rank) %>% 
  mutate(
    stock = toupper(stock)
  ) %>% 
  filter(!is.na(stock)) %>% 
  left_join(., stock_key, by = "stock") 

# stocks_to_add <- long_dat %>% 
#   filter(
#     is.na(Region1Name)
#   ) %>% 
#   select(
#     temp_stock, region
#   ) %>% 
#   arrange(temp_stock) %>% 
#   distinct()
# saveRDS(stocks_to_add, 
#         here::here("data", "rec", "southcoast_missing_stocks.rds"))


## export
saveRDS(long_dat, here::here("data", "rec", "rec_gsi.rds"))


## CLEAN SIZE ------------------------------------------------------------------

wide_size <- wide_rec2 %>% 
  select(id = biokey, fl, month_n, year, area, area_n, 
         region = cap_region) %>% 
  #remove missing and non-sensical size_classes
  filter(
    !is.na(fl),
    !fl == "remove"
  ) %>% 
  mutate(
    fl = as.numeric(fl),
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 450, 600, 750, 900, Inf), 
      labels=c("<45", "45-60", "60-75", "75-90", ">90")
    )
  ) 

# visualize sample coverage through space and time
size_n <- wide_size %>% 
  group_by(month_n, year, region) %>% 
  tally() 

ggplot(size_n) +
  geom_raster(aes(x = month_n, y = year, fill = n)) +
  facet_wrap(~region)

wide_size %>%
  group_by(size_bin, year, region) %>%
  summarize(n = length(unique(id))) %>%
  group_by(year, region) %>% 
  mutate(total_n = sum(n),
         ppn_obs = n / total_n) %>% 
  ggplot(., aes(x = as.factor(year), y = ppn_obs, fill = size_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_wrap(~region)


# # export size data as proportions
# rec_size_ppn_out <- size_dat %>%
#   mutate(sample_id = paste(month_n, region, jDay, year, sep = "_"),
#          region_c = region,
#          region = factor(abbreviate(region, minlength = 4)),
#          region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG")) %>%  
#   group_by(sample_id, region, region_c, year, month_n, size_bin) %>% 
#   summarize(sum_count = length(unique(id)), .groups = "drop") 

saveRDS(wide_size, here::here("data", "rec", "rec_size.rds"))


## CLEAN CATCH -----------------------------------------------------------------

rec_catch_raw <- read.csv(
  here::here("data", "rec", "rec_creel_nov15_21.csv"),
  stringsAsFactors = FALSE,
  na.strings=c("","NA"),
  header = FALSE
)
rec_catch1 <- rec_catch_raw[5:nrow(rec_catch_raw), ] 
names(rec_catch1) <- rec_catch_raw[4, ]
rec_catch1 <- janitor::clean_names(rec_catch1)

rec_catch <- rec_catch1 %>% 
  mutate(
    month = as.factor(month),
    month_n = match(month, month.name),
    year = as.numeric(year),
    region = case_when(
      subarea %in% c("13M", "13N") ~ "N. Strait of Georgia",
      area_n > 124 ~ "NWVI",
      area_n < 28 & area > 24 ~ "NWVI",
      area %in% c("20W", "20E", "20", "121", "21", "19JDF") ~ 
        "Juan de Fuca Strait",
      area_n < 125 & area_n > 120 ~ "SWVI",
      area_n < 25 & area_n > 20 ~ "SWVI",
      area %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      area %in% c("17", "18", "19", "19GST", "28", "29") ~ 
        "S. Strait of Georgia",
      area %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      area %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
    ),
    legal = ifelse(
      disposition == "Kept" | disposition == "Released Legal",
      "legal",
      "sublegal"
    ),
    kept_legal = case_when(
      disposition == "Kept" ~ "y",
      disposition == "Released Legal" ~ "n",
      legal == "sublegal" ~ NA
    )
  ) %>% 
  select(month, month_n, year, area, subarea, region, legal, kept_legal,
         adipose_clip)

tst <- match(unique(rec_catch1$month), month.name)

temp <- readRDS(here::here("data", "rec", "month_subarea_recCatch.rds"))

## OLD VERSIONS ----------------------------------------------------------------

# recreational composition data through 2019
# rec_raw <- readRDS(here::here("data", "rec", 
#                           "recIndProbsLong.rds")) 
# rec <- rec_raw %>% 
#   filter(legal == "legal",
#          !region %in% c("Queen Charlotte Sound")) %>%
#   mutate(
#     temp_strata = paste(month, region, sep = "_"),
#     sample_id = paste(temp_strata, jDay, year, sep = "_"),
#     min_m = case_when(
#       region %in% c("N. Strait of Georgia", "S. Strait of Georgia") ~ 1,
#       region == "Queen Charlotte and\nJohnstone Straits" ~ 6,
#       region == c("Juan de Fuca Strait") ~ 3
#     ),
#     max_m = case_when(
#       region %in% c("Juan de Fuca Strait", 
#                     "N. Strait of Georgia") ~ 9,
#       region == "S. Strait of Georgia" ~ 12,
#       region == "Queen Charlotte and\nJohnstone Straits" ~ 8
#     ),
#     coarse_agg = case_when(
#       pst_agg %in% c("WCVI", "NBC_SEAK") ~ "BC-coastal",
#       pst_agg %in% c("CA_ORCST", "WACST") ~ "US-coastal",
#       pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa") ~ "CR-fall",
#       pst_agg %in% c("CR-lower_sp", "CR-upper_sp") ~ "CR-spring",
#       TRUE ~ pst_agg
#     )
#   ) %>% 
#   group_by(region) %>% 
#   filter(!month_n < min_m,
#          !month_n > max_m) %>%
#   ungroup() %>% 
#   droplevels()
# 
# clean_rec <- rec %>% 
#   clean_comp(., grouping_col = "coarse_agg")


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


