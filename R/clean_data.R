## Data Clean
# Dec. 1, 2021
# Updated: May 8, 2023

library(tidyverse)
library(sf)
library(ggplot2)


# CRITICAL HABITAT CLEAN -------------------------------------------------------

if (Sys.info()['sysname'] == "Windows") {
  shp_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/"
} else {
  shp_path <- "/Users/cam/Google Drive/spatial"
}


# shapefile habitat boundaries (most conservative 70% boundary based on 
# posterior draws)
swift_shp <- st_read(
  here::here(shp_path, "srkw_foraging_areas", 
             "swiftsure.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp"))
haro_shp <- st_read(
  here::here(shp_path, "srkw_foraging_areas", 
             "haro.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp"))
hab_sf <- rbind(swift_shp, haro_shp) %>%
  # sf::st_transform(., crs = sp::CRS("+init=epsg:3005")) %>% 
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))


# coastline cropped 
coast_albers <- rbind(rnaturalearth::ne_states( "United States of America", 
                                                returnclass = "sf"), 
                      rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -125.5, ymin = 48.15, xmax = -122.25, ymax = 49.25) 


# pfma boundaries
pfma_areas <- st_read(
  here::here(shp_path, "pfma_subareas", "shape_files",
             "PFMA_Subareas_50k.shp")) %>% 
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  janitor::clean_names() %>% 
  #subset to southern
  filter(mgnt_area %in% c(123, 121, 21, 20, 19, 18)) 

# identify subareas w/ habitat
pfma_areas_sub <- pfma_areas %>% 
  # select pts in habitat
  st_intersection(., hab_sf)
pfma_areas$rkw_overlap <- ifelse(
  !(pfma_areas$label %in% pfma_areas_sub$label | pfma_areas$label == "19-3") |
    # exclude areass w/ very minimal overlap
    pfma_areas$label %in% c("19-3", "123-3", "18-6"),
  "no",
  "yes"
)

# check by plotting
ggplot() +
  geom_sf(data = coast_albers, colour = "black", fill = "white") + 
  geom_sf(data = hab_sf, colour = "red") +
  # geom_sf(data = pfma_areas %>% filter(mgnt_area %in% c("19", "18")), 
  #         aes(fill = as.factor(label))) +
  geom_sf(data = pfma_areas, aes(colour = rkw_overlap), fill = NA) +
  ggsidekick::theme_sleek()

saveRDS(
  hab_sf ,
  here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
)
saveRDS(
  pfma_areas ,
  here::here("data", "spatial", "pfma_subareas_sBC.rds")
)


# INDIVIDUAL DATA CLEAN --------------------------------------------------------

# recreational composition data since through 2021 (clean to match rec_raw)
rec_raw_new <- read_csv(
  here::here("data", "rec", "sc_biodata_aug_23.csv"),
  # here::here("data", "rec", "sc_biodata_jan_23.csv"),
  na = c("","NA"),
  # remove header
  skip = 6
) %>% 
  janitor::clean_names(.) %>% 
  # remove blanks
  filter(!is.na(biokey))


# remove duplicates, id numbers w/ multiple in columns, check unique vals are
# correct, merge with id2
wide_rec <- rec_raw_new %>% 
  mutate(
    new_pfma = str_replace_all(new_pfma, ".00", ""),
    area = case_when(
      is.na(area) ~ as.character(new_pfma),
      new_area == "US7" ~ "US7",
      TRUE ~ as.character(area)
    ),
    area_n = as.numeric(area),
    # separate northern areas of 13 (normally in JS) and add to NSoG
    cap_region = case_when(
      new_creel_subarea %in% c("13M", "13N") ~ "N. Strait of Georgia",
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
    date = as.POSIXct(collection_date, format="%d-%b-%Y"),
    month_n = lubridate::month(date),
    week_n = lubridate::week(date),
    # temporary key for identifying weird sizes below
    temp_key = paste(biokey, date, sep = "_"),
    # convert lat/lon to numeric
    lat = as.numeric(lat),
    lon = as.numeric(long) %>%
      ifelse(. > 0, -1 * ., .)
  ) 


# check to see if true duplicates by grouping by biokey then checking to see if 
# fork lengths and resolved stock id match
# dups <- wide_rec %>% 
#   group_by(biokey) %>%
#   filter(n() > 1) %>% 
#   mutate(row_id = row_number()) 
# 
# biokey_seq <- unique(dups$biokey)
# code_vec <- NULL
# for (i in seq_along(biokey_seq)) {
#   dum <- dups %>% filter(biokey == biokey_seq[i])
#   if (!is.na(dum$length_mm[1]) & dum$length_mm[1] != dum$length_mm[2]) {
#     code_vec <- c(code_vec, biokey_seq[i])
#   } else if (dum$resolved_stock_origin[1] != dum$resolved_stock_origin[2]) {
#     code_vec <- c(code_vec, biokey_seq[i])
#   }
# }
# 
# # all duplicates appear to be valid so remove second entry
# wide_rec <- wide_rec %>% 
#   group_by(biokey) %>% 
#   mutate(row_id = row_number()) %>% 
#   filter(!row_id == "2") %>% 
#   ungroup()


## correct some size entries
weird_sizes <- wide_rec %>%
  filter(length_mm < 150 | length_mm > 1500) %>%
  select(temp_key, length_mm, new_disposition, contains("size"))

corrected_sizes <- read.csv(
  here::here("data", "rec", "southcoast_size_errors_corrected.csv"),
  colClasses = "character"
) %>% 
  mutate(
    new_length_mm = ifelse(is.na(new_length_mm), "remove", new_length_mm)
  ) %>% 
  select(temp_key, new_length_mm)

# export then paste/add corrections into southcoast_size_errors_corrected.csv
write.csv(weird_sizes %>%
            filter(!temp_key %in% corrected_sizes$temp_key)
          ,
          here::here("data", "rec", "southcoast_size_errors.csv"),
          row.names = FALSE)

wide_rec3 <- full_join(wide_rec, corrected_sizes, by = "temp_key") %>%
  mutate(
    fl = ifelse(is.na(new_length_mm), length_mm, new_length_mm) %>% 
      as.numeric(),
    legal_lim = case_when(
      area_n < 20 & area_n > 11 ~ 620,
      area_n %in% c("28, 29") ~ 620,
      TRUE ~ 450
    ),
    legal = case_when(
      disposition == "Kept" ~ "legal",
      fl >= legal_lim ~ "legal",
      fl < legal_lim ~ "sublegal"
    ),
    fl = ifelse((fl < 150 | fl > 1500), NaN, fl),
    # shift locations to overlap with subareas
    lat = ifelse(new_location == "Cullite", 48.505, lat),
    lon = ifelse(new_location == "Moresby I.", -123.287, lon)
  ) 


rec_sf <- wide_rec3 %>%
  filter(!is.na(lat), !is.na(lon)) %>% 
  st_as_sf(
    ., 
    coords = c("lon", "lat"), 
    crs = sp::CRS("+proj=longlat +datum=WGS84")
  ) 

## identify whether samples caught inside/outside critical habitat
hab_sf <- readRDS(
  here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
)
rec_sf_sub <- st_intersection(rec_sf, hab_sf)

# assign areas with creel names
pfma_areas <- readRDS(
  here::here("data", "spatial", "pfma_subareas_sBC.rds")
)
rec_sf_areas <- st_intersection(rec_sf, pfma_areas)
rec_sf_areas_trim <- rec_sf_areas %>% 
  sf::st_drop_geometry() %>% 
  select(biokey, #new_location,
         subarea_new = label, 
         subarea_inc_rkw = rkw_overlap)


wide_rec4 <- wide_rec3 %>% 
  left_join(., rec_sf_areas_trim, by = "biokey") %>% 
  mutate(
    rkw_habitat1 = ifelse(
      biokey %in% rec_sf_sub$biokey, "yes", "no"
    ),
    # identify regions based on subarea boundaries and coherent habitat
    cap_region = case_when(
      subarea_inc_rkw == "yes" & lon < -124.25 ~ "swiftsure",
      subarea_inc_rkw == "yes" & lon < -123.5 & lon > -124.25 ~ "sooke",
      subarea_inc_rkw == "yes" & lon > -123.5 ~ "haro",
      TRUE ~ "outside"
    ),
    strata = case_when(
      subarea_new == "21-0" ~ "nitinat_nearshore",
      rkw_habitat1 == "yes" & subarea_new %in% c("121-1", "121-2") ~
        "nitinat_midshore",
      rkw_habitat1 == "no" & subarea_new %in% c("121-1", "121-2", "121-3") ~ "swiftsure",
      subarea_new %in% c("123-1", "123-2", "123-3") ~ "barkley_corner",
      rkw_habitat1 == "yes" & subarea_new %in% c("20-1", "20-3") ~ 
        "renfrew_habitat",
      rkw_habitat1 == "no" & subarea_new %in% c("20-1", "20-2", "20-3") ~ 
        "renfrew_nonhabitat",
      rkw_habitat1 == "yes" & subarea_new %in% c("20-5") ~ "sooke_habitat",
      rkw_habitat1 == "no" & subarea_new %in% c("20-4", "20-5", "20-6", "20-7") ~
        "sooke_nonhabitat",
      subarea_new %in% c("19-1", "19-2", "19-3") ~ "victoria",
      # combine s_haro habitat and nonhabitat given small number of samps
      # rkw_habitat1 == "yes" & subarea_new %in% c("19-4") ~ 
      #   "s_haro_habitat",
      subarea_new %in% c("19-4") ~ "s_haro",
      rkw_habitat1 == "yes" & subarea_new %in% c("19-5") ~ 
        "n_haro_habitat",
      rkw_habitat1 == "no" & subarea_new %in% 
        c("19-5", "19-6", "18-6", "18-10", "18-4", "18-5", "18-11", "18-2") ~ 
        "n_haro_nonhabitat",
      subarea_new %in% c("19-7", "18-7") ~ "saanich",
      area %in% c("18", "19") ~ "sog_outside",
      area == "123" ~ "wcvi_outside"
    ),
    rkw_habitat = ifelse(grepl("outside", strata), "outside", rkw_habitat1),
    whale_samples_time = ifelse(
      (year < 2011 | year > 2017) & month_n %in% c("6", "7", "8") & 
        rkw_habitat == "yes",
      "yes",
      "no"
    ),
    year = as.factor(year),
    yday = lubridate::yday(date),
    day_samp = paste(yday, year)
  )

# sample sizes
nrow(wide_rec4) #103k samples
nrow(wide_rec4 %>% 
       filter(!is.na(resolved_stock_origin))) #79k stock
nrow(wide_rec4 %>% 
       filter(!is.na(resolved_stock_origin) & !is.na(lat))) #69k stock/loc 

wide_rec4 %>% 
  filter(!is.na(resolved_stock_origin) & !is.na(lat)) %>% 
  group_by(cap_region) %>% 
  tally()


saveRDS(wide_rec4, here::here("data", "rec", "wide_rec.rds"))


# GSI CLEAN --------------------------------------------------------------------

wide_rec4 <- readRDS(here::here("data", "rec", "wide_rec.rds"))

# stock key
stock_key <- readRDS(here::here("data", "rec", "finalStockList_Jul2023.rds")) %>%
  janitor::clean_names() %>% 
  mutate(
    stock_group = case_when(
      pst_agg %in% c("CA_ORCST", "CR-lower_fa", "CR-lower_sp", "CR-upper_su/fa",
                     "WACST", "Russia", "CR-upper_sp", "NBC_SEAK") ~ "other",
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN") ~ "ECVI_SOMN",
      # grepl(".2", region1name) ~ "Fraser Yearling",
      grepl("Fraser", region1name) ~ region1name,
      # region1name == "Fraser_Summer_4.1" ~ "Fraser Summer 4.1",
      # region1name == "Fraser_Fall" ~ "Fraser Fall",
      TRUE ~ pst_agg
    )
  )

# trim for GSI purposes
wide_rec4_trim <- wide_rec4 %>% 
  filter(
    !is.na(resolved_stock_source)
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
    id = biokey, date, week_n, month_n, year, cap_region, area, area_n, 
    fishing_site = new_location, subarea = subarea_new, strata,
    lat, lon, rkw_habitat, cap_region, subarea_inc_rkw, whale_samples_time,
    legal, fl, sex, ad = adipose_fin_clipped,
    age = resolved_age, resolved_stock_source, 
    stock_1, stock_2 = dna_stock_2, stock_3 = dna_stock_3,
    stock_4 = dna_stock_4, stock_5 = dna_stock_5,
    starts_with("prob"),
    # retain region data to help build stock key (will eventually be removed)
    ends_with("rollup")
  ) 

# replace 0 probabilities with v. small values (just to be safe); then recalc
# ppns
prbs <- wide_rec4_trim %>% 
  select(starts_with("prob")) %>% 
  as.matrix()
prbs[prbs == 0] <- .00001
row_sums <- apply(prbs, 1, sum, na.rm = T)
new_prbs <- prbs / row_sums


#pivot to long (probs and stock IDs separately) and join 
probs <- wide_rec4_trim %>% 
  # replace with updated probabilities from above
  select(-starts_with("prob")) %>% 
  cbind(., new_prbs) %>% 
  pivot_longer(., cols = starts_with("prob"), names_to = "rank", 
               names_pattern = "prob_(.+)",
               values_to = "prob") %>%
  select(id, rank, prob)

regions <- wide_rec4_trim %>% 
  select(id, starts_with("region")) %>% 
  pivot_longer(., cols = starts_with("region"), names_to = "rank", 
               names_pattern = "region_(.+)_rollup",
               values_to = "region") %>%
  select(id, rank, region)

long_rec <- wide_rec4_trim %>% 
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
    creelsub = gsub("(\\d+)([[:alpha:]])", "\\1-\\2", subarea)
  ) %>% 
  filter(!is.na(stock),
         !is.na(prob)) %>% 
  left_join(., stock_key, by = "stock") #%>%
  # # add distance from juan de fuca mouth from creel subarea spatial file
  # left_join(., creel_spatial %>% select(creelsub, dist_123i),
  #           by = "creelsub")


# check for missing regional assignments
long_rec %>%
  filter(
    is.na(region1name)
  ) %>%
  select(
    stock, region
  ) %>%
  arrange(stock) %>%
  distinct()


# saveRDS(
#   long_rec %>%
#     select(
#       stock, region
#     ) %>%
#     arrange(stock) %>%
#     distinct(),
#   here::here("data", "rec", "southcoast_updated_stocks.rds")
# )

## export
saveRDS(long_rec, here::here("data", "rec", "rec_gsi.rds"))


## CLEAN SIZE ------------------------------------------------------------------

wide_rec4 <- readRDS(here::here("data", "rec", "wide_rec.rds"))

wide_size <- wide_rec4 %>% 
  select(id = biokey, date, week_n, month_n, year, cap_region, area, area_n,
         fishing_site = new_location, subarea = subarea_new, strata,
         lat, lon, rkw_habitat, subarea_inc_rkw,
         whale_samples_time, legal, fl, sex, ad = adipose_fin_clipped,
         resolved_stock_source, resolved_stock_region) %>%
  #remove missing and non-sensical size_classes
  filter(
    !is.na(fl)
  ) %>% 
  mutate(
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 551, 701, 851, Inf), 
      labels = c("<55", "55-70", "70-85", ">85")
    ),
    size_bin2 = factor(
      size_bin, labels = c("sublegal", "small", "medium", "large")
    ),
    month_n = lubridate::month(date)
  ) 

# visualize sample coverage through space and time
size_n <- wide_size %>% 
  group_by(month_n, year, cap_region) %>% 
  tally() 

ggplot(size_n) +
  geom_raster(aes(x = month_n, y = year, fill = n)) +
  facet_wrap(~cap_region)

wide_size %>%
  group_by(size_bin, month_n, cap_region) %>%
  summarize(n = length(unique(id))) %>%
  group_by(month_n, cap_region) %>% 
  mutate(total_n = sum(n),
         ppn_obs = n / total_n) %>% 
  ggplot(., aes(x = as.factor(month_n), y = ppn_obs, fill = size_bin)) +
  geom_bar(position="stack", stat="identity") +
  ggsidekick::theme_sleek() +
  facet_wrap(~cap_region)

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


## ADD LOCATION DATA -----------------------------------------------------------

# add georeferenced locations for stocks along main migratory corridor
fish_locs <- read.csv(here::here("data", "spatial", "sc_rec_locations.csv"),
                      na.strings=c("", "#N/A", "-")) %>% 
  select(aa_biodata_location, new_name:new_lat) %>% 
  filter(!(grepl("None", aa_biodata_location)),
         new_pfma %in% c(#"123", "23",
                         "121", "21", "22", "20", "19", "18")) %>% 
  distinct()


sp_rec <- wide_rec3 %>% 
  filter(area %in% c(#"123", "23", 
    "121", "21", "22", "20", "19", "18")) %>%
  select(biokey, area, aa_biodata_location = fishing_location) %>% 
  left_join(., fish_locs, by = "aa_biodata_location") %>% 
  glimpse()
