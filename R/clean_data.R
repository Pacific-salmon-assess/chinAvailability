## Data Clean
# Dec. 1, 2021

library(tidyverse)

# stock key 
stock_key <- readRDS(here::here("data", "rec", "finalStockList_Jan2022.rds"))

# spatial data for creel subareas in southern bc to use as a covariate
creel_spatial <- readRDS(
  here::here("data", "rec", "creel_subarea_spatial.rds")
) %>% 
  sf::st_drop_geometry() %>% 
  #coerce 20-D to equal specific subarea
  mutate(
    creelsub = ifelse(creelsub == "20-DO", "20-D", creelsub)
  )


# INDIVIDUAL DATA CLEAN --------------------------------------------------------

# recreational composition data since through 2021 (clean to match rec_raw)
rec_raw_new <- read.csv(here::here("data", "rec", "sc_biodata_jul8_21.csv"),
                        stringsAsFactors = FALSE, na.strings=c("","NA")) %>% 
  janitor::clean_names(.) 

wide_rec <- rec_raw_new %>% 
  # change US area 7 (near San Juan island) to SSoG
  mutate(
    area = ifelse(area == "US7", "19GST", area),
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
    stock = toupper(stock),
    subarea = ifelse(subarea %in% c("US7"), "19B", subarea),
    creelsub = gsub("(\\d+)([[:alpha:]])", "\\1-\\2", subarea)
  ) %>% 
  filter(!is.na(stock),
         !is.na(prob)) %>% 
  left_join(., stock_key, by = "stock") %>% 
  # add distance from juan de fuca mouth from creel subarea spatial file
  left_join(., creel_spatial %>% select(creelsub, dist_123i),
            by = "creelsub")

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
saveRDS(long_rec, here::here("data", "rec", "rec_gsi.rds"))


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

rec_creel_raw <- read.csv(
  here::here("data", "rec", "rec_creel_aug10_21.csv"),
  stringsAsFactors = FALSE,
  na.strings=c("","NA")#,
  # header = FALSE
)
# rec_catch1 <- rec_creel_raw[5:nrow(rec_creel_raw), ] 
# names(rec_catch1) <- rec_creel_raw[4, ]
rec_creel_raw <- janitor::clean_names(rec_creel_raw) %>% 
  #focus only on chinook and effort estimates
  filter(
    species %in% c("BOAT TRIPS", "CHINOOK SALMON")
  )

# note that subarea-specific estimates are not always available 
# (e.g. in 2014 estimates for 125A and 125G available for Jul/Aug,
# 125 total (?) estimates for June/Sep)
rec_creel <- rec_creel_raw %>% 
  mutate(
    month = as.factor(month),
    month_n = match(month, month.name),
    area_n = as.numeric(str_replace_all(pfma, "[:letter:]", "")),
    year = as.numeric(year),
    region = case_when(
      creel_sub_area %in% c("13M", "13N") ~ "N. Strait of Georgia",
      area_n > 124 ~ "NWVI",
      area_n < 28 & area_n > 24 ~ "NWVI",
      area_n %in% c("20", "21", "22", "121") | 
        creel_sub_area %in% c("Area 19 (JDF)", "19D", "19E", "19C") ~ 
        "Juan de Fuca Strait",
      area_n < 125 & area_n > 120 ~ "SWVI",
      area_n < 25 & area_n > 20 ~ "SWVI",
      area_n %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      area_n %in% c("17", "18", "28", "29") | 
        creel_sub_area %in% c("19A", "19B", "Area 19 (GS)") ~ 
        "S. Strait of Georgia",
      area_n %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      area_n %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
    ),
    reg_f = as.factor(abbreviate(region, 4)),
    area = ifelse(
      area_n %in% c("13", "19"),
      paste(area_n, reg_f, sep = "_"),
      as.character(area_n)
    ),
    #is a subarea specific estimate available
    subarea_est = ifelse(grepl("Area", creel_sub_area), "n", "y")
  ) 
  

# separate catch/effort and rejoin
effort_subarea <- rec_creel %>% 
  filter(species == "BOAT TRIPS",
         !is.na(estimate)) %>% 
  mutate(effort = as.numeric(estimate)) %>%
  select(month, year, subarea = creel_sub_area, area, effort)
  
catch_subarea <- rec_creel %>% 
  filter(species == "CHINOOK SALMON") %>% 
  mutate(
    legal = ifelse(
      disposition == "Kept" | disposition == "Released Legal",
      "legal",
      "sublegal"
    ),
    kept_legal = case_when(
      disposition == "Kept" ~ "y",
      disposition == "Released Legal" ~ "n",
      legal == "sublegal" ~ NA_character_
    ),
    catch = as.numeric(estimate)
  ) %>% 
  select(month, month_n, year, area_n, area, subarea = creel_sub_area, region,
         subarea_est, legal, kept_legal, adipose_mark, catch) %>% 
  left_join(., effort_subarea, by = c("month", "year", "subarea", "area"))


# check no duplicates
catch_subarea %>% 
  group_by(kept_legal, legal, adipose_mark, subarea, month, year) %>% 
  tally() %>% 
  filter(n > 1)


# calculate area totals; effort replicated among catch categories so sum 
# separately then rejoin
# NOTE pools adipose marked and kept/released LEGAL fish
effort_area <- effort_subarea %>% 
  group_by(month, year, area) %>% 
  summarize(
    effort = sum(effort, na.rm = T),
    .groups = "drop"
  )
sum(effort_area$effort) == sum(effort_subarea$effort)
  
catch_area <- catch_subarea %>% 
  mutate(reg = abbreviate(region, 4)) %>% 
  group_by(month, month_n, year, area, reg, region, legal) %>% 
  summarize(
    catch = sum(catch, na.rm = T),
    .groups = "drop"
  ) %>% 
  left_join(., effort_area, by = c("month", "year", "area"))
sum(catch_subarea$catch, na.rm = T) == sum(catch_area$catch, na.rm = T)


# check to ensure one estimate per month
catch_area %>% 
  group_by(legal, area, month, year) %>% 
  tally() %>% 
  filter(n > 1)


saveRDS(catch_subarea, here::here("data", "rec", "rec_creel_subarea.rds"))
saveRDS(catch_area, here::here("data", "rec", "rec_creel_area.rds"))


