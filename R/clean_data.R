## Data Clean
# Dec. 1, 2021
# Updated: March 12, 2024

library(tidyverse)
library(sf)
library(ggplot2)


# IMPORT STOCK KEY AND PREP ----------------------------------------------------

# stock key
stock_key <- readRDS(here::here("data", "rec", "finalStockList_Nov2024.rds")) %>%
  janitor::clean_names() %>% 
  mutate(
    stock_group = case_when(
      # move Juan de Fuca populations and Boundary Bay populations to Puget
      region1name == "Juan_de_Fuca" | grepl("NICOMEK", stock) | 
        grepl("SERPEN", stock) | 
        (grepl("CAMPB", stock) & region1name == "Fraser_Fall") ~ "PSD",
      pst_agg %in% c("CR-lower_fa", "CR-upper_su/fa") | 
        region1name == "Willamette_R" ~ 
        "Col_Summer_Fall",
      pst_agg %in% c("CR-lower_sp", "CR-upper_sp")  ~ "Col_Spring",
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN", "NEVI") ~
        "ECVI_SOMN",
      pst_agg %in% c("CA_ORCST", "WACST", "Russia", "NBC_SEAK", "Yukon") ~
        "other",
      grepl("Fraser", region1name) ~ region1name,
      TRUE ~ pst_agg
    ) %>% 
      factor(
        .,
        levels = c("other", "Col_Spring", "Col_Summer_Fall", "PSD",  
                   "WCVI", "ECVI_SOMN", "Fraser_Spring_4.2",
                   "Fraser_Spring_5.2", "Fraser_Summer_5.2", "Fraser_Summer_4.1",
                   "Fraser_Fall"),
        labels = c("other", "Col_Spr", "Col_Sum/Fall", "PSD", "WCVI", 
                   "ECVI_SOMN", "FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2", 
                   "FR_Sum_4.1", "FR_Fall")
      )
  )


# CRITICAL HABITAT CLEAN -------------------------------------------------------

if (Sys.info()['sysname'] == "Windows") {
  shp_path <- "G:/My Drive/spatial"
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
# ggplot() +
#   geom_sf(data = coast_albers, colour = "black", fill = "white") + 
#   geom_sf(data = hab_sf, colour = "red") +
#   # geom_sf(data = pfma_areas %>% filter(mgnt_area %in% c("19", "18")), 
#   #         aes(fill = as.factor(label))) +
#   geom_sf(data = pfma_areas, aes(colour = rkw_overlap), fill = NA) +
#   ggsidekick::theme_sleek()
# 
# saveRDS(
#   hab_sf ,
#   here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
# )
# saveRDS(
#   pfma_areas ,
#   here::here("data", "spatial", "pfma_subareas_sBC.rds")
# )


# INDIVIDUAL DATA CLEAN --------------------------------------------------------

# clean recreational composition data through 2022 (clean to match rec_raw)
rec_raw_new <- read_csv(
  here::here("data", "rec", "raw_data", "sc_biodata_aug_23.csv"),
  na = c("","NA"),
  # remove header
  skip = 6
) %>% 
  janitor::clean_names(.) %>% 
  # remove blanks
  filter(!is.na(biokey)) %>% 
  mutate(
    # convert lat/lon to numeric
    lat = as.numeric(lat),
    lon = as.numeric(long) %>%
      ifelse(. > 0, -1 * ., .)
  )


## clean spatial data
site_ll <- rec_raw_new %>% 
  select(fishing_location, area, lat, lon) %>% 
  distinct() %>% 
  filter(
    !is.na(lon) 
  ) %>% 
  #remove duplicates that have very similar lat and lon
  mutate(
    dd = paste(fishing_location, area, sep = "-")
  ) %>% 
  group_by(dd) %>% 
  slice_head(n = 1) %>% 
  ungroup() %>% 
  select(-dd)

# add distance to coast estimates 
# coast <- readRDS(here::here("data", "spatial", "coast_sf.RDS"))
# coast_dist <- geosphere::dist2Line(
#   p = site_ll %>% select(lon, lat), 
#   line = as(coast, 'Spatial')
# )
# site_ll$shore_dist <- coast_dist[, "distance"]

rec_raw_new <- left_join(
  rec_raw_new, 
  site_ll %>% select(fishing_location, area
                     ),
  by = c("fishing_location", "area")
)

# most recent push from AP, 2022 and 2023 data 
rec_raw_new24 <- read_csv(
  here::here("data", "rec", "raw_data", "sc_biodata_jun_24.csv"),
  na = c("","NA")
) %>% 
  janitor::clean_names(.) %>% 
  # remove blanks, samples already present above, and non sport caught fish
  filter(!is.na(biokey),
         !biokey %in% rec_raw_new$biokey,
         sample_type == "Sport") %>% 
  left_join(., site_ll, by = c("area", "fishing_location")) %>% 
  select(-c(scale_submission_number))


# remove duplicates, id numbers w/ multiple in columns, check unique vals are
# correct, merge with id2
wide_rec <- rec_raw_new %>% 
  select(colnames(rec_raw_new24)) %>% 
  rbind(., rec_raw_new24) %>% 
  mutate(
    date = as.POSIXct(collection_date, format="%d-%b-%Y"),
    month_n = lubridate::month(date),
    week_n = lubridate::week(date),
    # # temporary key for identifying weird sizes below
    temp_key = paste(biokey, date, sep = "_"),
    pbt_brood_year_n = case_when(
      pbt_brood_year %in% c("GSI 0000", "Not Loaded", "0") |
        is.na(pbt_brood_year) ~  NaN,
      TRUE ~ as.numeric(pbt_brood_year)
    ),
    # in some cases resolved age deviates from cwt/pbt brood year; replace w/
    # corrected version
    age = case_when(
      !is.na(cwt_brood_year) ~ year - cwt_brood_year,
      !is.na(pbt_brood_year_n) ~ year - pbt_brood_year_n,
      TRUE ~ resolved_age
    ),
    age_source = case_when(
      !is.na(cwt_brood_year) ~ "cwt",
      !is.na(pbt_brood_year_n) ~ "pbt",
      is.na(pbt_brood_year_n) & is.na(cwt_brood_year) & !is.na(resolved_age) ~ 
        "scale"
    ),
    pbt = ifelse(
      is.na(pbt_brood_year_n),
      FALSE,
      TRUE
    )
  ) 


# export data for age classification error analysis
# rec_for_aging <- wide_rec %>% 
#   filter(
#     (!is.na(cwt_brood_year) & !is.na(age_gr)) |
#       (!is.na(pbt_brood_year_n) & !is.na(age_gr)),
#     # remove uncertain age assignments
#     !grepl("M", age_gr),
#     !grepl("F", age_gr),
#     !grepl("S", age_gr),
#     !grepl("R", age_gr),
#     # remove missing stokcs
#     !(is.na(resolved_stock_origin) | resolved_stock_rollup == "SUS (assumed)")
#   ) %>% 
#   mutate(
#     stock = case_when(
#       resolved_stock_source == "DNA" ~ dna_results_stock_1,
#       resolved_stock_source == "CWT" ~ cwt_result,
#       resolved_stock_source == "Otolith Stock" ~ oto_stock,
#     ) %>% 
#       toupper(),
#     pbt_yes = ifelse(!is.na(pbt_brood_year_n), "yes", "no"),
#     cwt_yes = ifelse(!is.na(cwt_brood_year), "yes", "no")
#   ) %>% 
#   left_join(., stock_key %>% select(stock, stock_group), by = "stock") %>% 
#   select(
#     biokey, year, stock, resolved_stock_rollup, stock_group, age, age_gr, 
#     age_source
#   )
# saveRDS(rec_for_aging, here::here("data", "rec", "aging_data.rds"))


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
# weird_sizes <- wide_rec %>%
#   filter(length_mm < 150 | length_mm > 1500) %>%
#   select(temp_key, length_mm, disposition, contains("size"))

corrected_sizes <- read.csv(
  here::here("data", "rec", "southcoast_size_errors_corrected.csv"),
  colClasses = "character"
) %>% 
  mutate(
    new_length_mm = ifelse(is.na(new_length_mm), "remove", new_length_mm)
  ) %>% 
  select(temp_key, new_length_mm)

# export then paste/add corrections into southcoast_size_errors_corrected.csv
# write.csv(weird_sizes %>%
#             filter(!temp_key %in% corrected_sizes$temp_key)
#           ,
#           here::here("data", "rec", "southcoast_size_errors.csv"),
#           row.names = FALSE)


# define marine age based on total age readings
# estimates seem suspicious (e.g. lots of 4_1s for for Spring 4_2s and fall run 
# stocks); DONT USE
# sw_age_key <- wide_rec %>% 
#   select(temp_key, age_gr)  %>% 
#   filter(!is.na(age_gr)) %>% 
#   mutate(
#     total_age = purrr::map(
#       age_gr, ~ strsplit(.x, "")[[1]][1]
#     ) %>% 
#       as.numeric(),
#     fw_age = purrr::map(
#       age_gr, ~ strsplit(.x, "")[[1]][2] 
#     ) %>% 
#       as.numeric(),
#     sw_age = case_when(
#       grepl("M", age_gr) ~ total_age,
#       is.numeric(total_age) & is.numeric(fw_age) ~ (total_age - fw_age),
#       TRUE ~ NaN
#     )
#   ) %>% 
#   select(temp_key, sw_age)


wide_all <- full_join(wide_rec, corrected_sizes, by = "temp_key") %>%
  # left_join(., sw_age_key, by = "temp_key") %>% 
  mutate(
    fl = ifelse(is.na(new_length_mm), length_mm, new_length_mm) %>% 
      as.numeric(),
    fl = ifelse((fl < 150 | fl > 1500), NaN, fl),
    # shift locations to overlap with subareas
    lat = ifelse(fishing_location == "Cullite", 48.505, lat),
    lon = ifelse(fishing_location == "Moresby I.", -123.287, lon),
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
    id = biokey, date, week_n, month_n, year, area, 
    fishing_site = fishing_location,
    subarea, lat, lon, #shore_dist,
    fl, ad = adipose_fin_clipped, disposition,
    pbt, pbt_brood_year_n, cwt_brood_year_n = cwt_brood_year,
    age_gr, resolved_stock_source, 
    stock_1, stock_2 = dna_stock_2, stock_3 = dna_stock_3,
    stock_4 = dna_stock_4, stock_5 = dna_stock_5,
    starts_with("prob")
  ) %>% 
  mutate(
    area_n = as.numeric(area),
    legal_lim = case_when(
      area_n < 20 & area_n > 11 ~ 620,
      area_n %in% c("28, 29") ~ 620,
      TRUE ~ 450
    ),
    legal = case_when(
      fl >= legal_lim ~ "legal",
      fl < legal_lim ~ "sublegal"
    )
  )


rec_sf <- wide_all %>%
  filter(!is.na(lat), !is.na(lon)) %>% 
  st_as_sf(
    ., 
    coords = c("lon", "lat"), 
    crs = sp::CRS("+proj=longlat +datum=WGS84")
  ) 

## identify whether samples caught inside/outside critical habitat
# NO LONGER NECESSARY WITH SPATIALLY EXPLICIT MODEL
# hab_sf <- readRDS(
#   here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
# )
# rec_sf_sub <- st_intersection(rec_sf, hab_sf)

# assign areas with creel names
pfma_areas <- readRDS(
  here::here("data", "spatial", "pfma_subareas_sBC.rds")
)
rec_sf_areas <- st_intersection(rec_sf, pfma_areas)
rec_sf_areas_trim <- rec_sf_areas %>%
  sf::st_drop_geometry() %>%
  select(id, #new_location,
         subarea_new = label)


wide_rec4 <- wide_all %>% 
  left_join(., rec_sf_areas_trim, by = "id") %>%
  mutate(
    # redefine strata manually
    strata = case_when(
      lat > 48.8 ~ "other",
      (lat < 48.65 & lon < -124.9) | (lat < 48.55 & lon < -124.7) |
        (lon < -125.3) | (lat < 48.7 & lon < -125) ~ "swiftsure",
      lat > 48.61 & lon < -124.78 ~ "swiftsure_nearshore",
      lon > -124.8 & lon < -124.2 ~ "renfrew",
      lat < 48.52 ~ "vic",
      subarea %in% c("18A", "19A") ~ "saanich", 
      lat > 48.5 & lat < 48.8 & lon < -122.9 ~ "haro",
      TRUE ~ "other"
    ),
    strata2 = case_when(
      strata == "vic" & lon > -123.55 ~ "haro",
      TRUE ~ strata
    )
  )

# fix issue with specific samples collected in same location but assigned to 
# different strata
wide_rec4$strata <- ifelse(wide_rec4$fishing_site == "Cape Kepple",
                           "saanich",
                           wide_rec4$strata)

#check disposition
# wide_rec4 %>%
#   filter(!strata == "other") %>%
#   group_by(year, strata2, week_n, month_n) %>%
#   mutate(count = length(unique(id)),
#          released = ifelse(disposition == "Released", 1, 0)) %>%
#   group_by(year, strata2, week_n, month_n, released, count) %>%
#   summarize(
#     ppn_rel = sum(released) / count
#   ) %>%
#   distinct() %>%
#   ggplot(.) +
#   geom_jitter(aes(x = as.factor(month_n), ppn_rel)) +
#   facet_wrap(~strata2)
#   

#checks for mapping locations 
# coast <- readRDS(
#   here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>%
#   # sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>%
#   sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.75, ymax = 48.8)
# ggplot() +
#   geom_sf(data = coast, color = "black", fill = NA) +
#   # geom_sf(data = pfma_areas, aes(colour = rkw_overlap), fill = NA) +
#   ggsidekick::theme_sleek() +
#   scale_x_continuous(expand = c(0, 0)) +
#   scale_y_continuous(expand = c(0, 0)) +
#   geom_point(
#     data = wide_rec4  %>%
#       filter(!strata2 == "other") %>%
#       group_by(
#         fishing_site, lat, lon, strata2
#       ) %>%
#       summarize(
#         n = length(unique(id))
#       ),
#     aes(x = lon, y = lat, size = n, fill = strata2),
#     alpha = 0.7,
#     shape = 21
#   )


# sample sizes
nrow(wide_rec4) #119k samples
nrow(wide_rec4 %>% 
       filter(!is.na(stock_1))) #80k stock
nrow(wide_rec4 %>% 
       filter(!is.na(stock_1) & !is.na(lat))) #70k stock/loc 


saveRDS(wide_rec4 %>% select(-fishing_site),
        here::here("data", "rec", "wide_rec.rds"))


# GSI CLEAN --------------------------------------------------------------------

# import PBT estimates to estimate coverage 
# pbt_rate <- readRDS(
#   here::here("data", "sep", "mean_pbt_rate.rds")
# ) %>% 
#   mutate(
#     pbt_stock = ifelse(
#       collection_extract == "SHUSWAP_RIVER-LOWER (includes Kingfisher)",
#       "SHUSWAP_RIVER_LOWER",
#       gsub("-", "_", collection_extract) %>% 
#         toupper()
#     )
#   ) %>% 
#   select(
#     pbt_stock, brood_year = year, tag_rate
#   ) 
# 
# # ID stocks that have had > 80% in more than 5 years 
# high_rate <- pbt_rate %>% 
#   filter(tag_rate > 0.8) %>% 
#   group_by(pbt_stock) %>% 
#   mutate(n_year = length(unique(brood_year))) %>% 
#   filter(n_year > 5) %>% 
#   pull(pbt_stock) %>% 
#   unique()
# pbt_rate$high <- ifelse(pbt_rate$pbt_stock %in% high_rate, TRUE, FALSE)
# saveRDS(pbt_rate, here::here("data", "sep", "cleaned_pbt.rds"))

pbt_rate <- readRDS(here::here("data", "sep", "cleaned_pbt.rds"))

wide_rec4_trim <- readRDS(here::here("data", "rec", "wide_rec.rds")) %>% 
  filter(
    !is.na(stock_1)
  )


#pivot to long (probs and stock IDs separately) and join 
probs <- wide_rec4_trim %>% 
  # replace with updated probabilities from above
  pivot_longer(., cols = starts_with("prob"), names_to = "rank", 
               names_pattern = "prob_(.+)",
               values_to = "prob") %>%
  select(id, rank, prob)

long_rec <- wide_rec4_trim %>% 
  select(-starts_with("prob")) %>% 
  pivot_longer(., cols = starts_with("stock"), names_to = "rank", 
               names_pattern = "stock_(.+)",
               values_to = "stock") %>% 
  left_join(., probs, by = c("id", "rank")) %>% 
  arrange(desc(date), id, desc(prob)) %>% 
  mutate(
    stock = toupper(stock),
    creelsub = gsub("(\\d+)([[:alpha:]])", "\\1-\\2", subarea_new)
  ) %>% 
  filter(!is.na(stock),
         !is.na(prob)) %>% 
  left_join(., stock_key, by = "stock") %>% 
  # add initial PBT id
  mutate(
    pbt_stock = gsub("-", "_", stock) %>% 
      gsub(" H", "", .) %>% 
      gsub(" ", "_", .) %>% 
      toupper(),
    # replace stocks with mismatched names
    pbt_stock = case_when(
      grepl("CAPILANO", pbt_stock) ~ "CAPILANO_RIVER",
      stock %in% c("CHILLIWACK_RIVER_SUMMER") ~ "CHILLIWACK_RIVER_SUMMER",
      grepl("CHILLIWAC", pbt_stock) ~ "CHILLIWACK_RIVER_FALL",
      grepl("HARRISON", pbt_stock) ~ "HARRISON_RIVER",
      grepl("HORSEFLY", pbt_stock) ~ "HORSEFLY_RIVER",
      grepl("COLDWATER", pbt_stock) ~ "COLWDWATER_RIVER",
      grepl("INDIANPOINT", pbt_stock) ~ "INDIANPOINT_CREEK",
      pbt_stock %in% c("SHUSWAP_RIVER,_MIDDLE,") ~ "SHUSWAP_RIVER_MIDDLE",
      grepl("CHEAKAMUS", pbt_stock) ~ "CHEAKAMUS_RIVER",
      grepl("CHEHALIS_RIVER", pbt_stock) ~ "CHEHALIS_RIVER_SUMMER",
      pbt_stock %in% c("L_QUALICUM") ~ "LITTLE_QUALICUM_RIVER",
      grepl("QUALICUM", pbt_stock) ~ "QUALICUM_RIVER",
      grepl("ATNARKO", pbt_stock) ~ "ATNARKO_RIVER",
      grepl("BURMAN", pbt_stock) ~ "BURMAN_RIVER",
      grepl("MCGREGOR", pbt_stock) ~ "MCGREGOR_RIVER",
      grepl("MORKILL", pbt_stock) ~ "MORKILL_RIVER",
      grepl("PHILLIPS", pbt_stock) ~ "PHILLIPS_RIVER",
      grepl("NAZKO", pbt_stock) ~ "NAZKO_RIVER",
      grepl("CONUMA", pbt_stock) ~ "CONUMA_RIVER",
      grepl("NECHAKO", pbt_stock) ~ "NECHAKO_RIVER",
      grepl("NICOLA", pbt_stock) ~ "NICOLA_RIVER",
      grepl("PORTAGE", pbt_stock) ~ "PORTAGE_CREEK",
      grepl("COWICHAN", pbt_stock) ~ "COWICHAN_RIVER",
      grepl("QUESNEL", pbt_stock) ~ "QUESNEL_RIVER",
      grepl("SERPENTINE", pbt_stock) ~ "SERPENTINE_RIVER",
      grepl("TORPY", pbt_stock) ~ "TORPY_RIVER",
      grepl("SLIM_C", pbt_stock) ~ "SLIM_CREEK",
      grepl("SPIUS", pbt_stock) ~ "SPIUS_CREEK",
      grepl("SWIFT", pbt_stock) ~ "SWIFT_CREEK",
      pbt_stock == "SALMON_RIVER(THOM)" ~ "SALMON_RIVER_SOTH",
      pbt_stock %in% 
        c("S_FIRST_LK/GSVI", "S_NANAIMO_R", "NANAIMO_RIVER_FALL", "NANAIMO",
          "NANAIMO_RIVER_UPPER", "NANAIMO_RIVER_F", "NANAIMO_F", 
          "NANAIMOUPPER") ~ 
        "NANAIMO_RIVER_FALL", 
      pbt_stock %in% c("NANAIMO_RIVER_S", "NANAIMO_SU", "NANAIMOSU_BACKX") ~ 
        "NANAIMO_RIVER_SUMMER",       
      grepl("ASHLU", pbt_stock) ~ "ASHLU_CREEK",
      grepl("NAHMINT", pbt_stock) ~ "NAHMINT_RIVER",
      grepl("THORNTON", pbt_stock) ~ "THORNTON_CREEK",
      grepl("GOLD", pbt_stock) ~ "GOLD_RIVER",
      grepl("SOOKE", pbt_stock) ~ "SOOKE_RIVER",
      grepl("NITINAT", pbt_stock) ~ "NITINAT_RIVER",
      grepl("MAMQUAM", pbt_stock) ~ "MAMQUAM_RIVER",
      grepl("SHOVEL", pbt_stock) ~ "SHOVELNOSE_CREEK",
      grepl("ROBERTSON", pbt_stock) ~ "ROBERTSON_CREEK",
      grepl("TAHSIS", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("TAHSISH", pbt_stock) ~ "TAHSIS_RIVER",
      grepl("CHEMAINUS", pbt_stock) ~ "CHEMAINUS_RIVER",
      grepl("PUNTLEDGE", pbt_stock) ~ "PUNTLEDGE_RIVER_FALL",
      pbt_stock == "S_SALMON_R/JNST" ~ "SALMON_RIVER_JNST",
      grepl("SARITA", pbt_stock) ~ "SARITA_RIVER",
      grepl("TLUPANA", pbt_stock) ~ "TLUPANA_RIVER",
      grepl("COLDWATER", stock) ~ "COLDWATER_RIVER",
      grepl("DEADMAN", stock) ~ "DEADMAN_RIVER",
      pbt_stock == "COTTONWOOD_RIVER_LOWER" ~ "COTTONWOOD_RIVER_(LOWER)",
      grepl("CARIBOO", pbt_stock) ~ "CARIBOO",
      grepl("ENDAKO", pbt_stock) ~ "ENDAKO_RIVER",
      grepl("LEINER", pbt_stock) ~ "LEINER_RIVER",
      grepl("WILLOW", pbt_stock) ~ "WILLOW_RIVER",
      grepl("MARBLE", pbt_stock) ~ "MARBLE_RIVER",
      grepl("QUINSAM", pbt_stock) ~ "QUINSAM_RIVER",
      grepl("WOSS", pbt_stock) ~ "WOSS_RIVER",
      grepl("KENNEDY", pbt_stock) ~ "KENNEDY_RIVER_LOWER",
      grepl("NIMPKISH", pbt_stock) ~ "NIMPKISH_RIVER",
      grepl("TRANQUIL", pbt_stock) ~ "TRANQUIL_CREEK",
      grepl("TOQUART", pbt_stock) ~ "TOQUART_RIVER",
      grepl("BEDWELL", pbt_stock) ~ "BEDWELL_RIVER",
      grepl("SAN_JUAN", pbt_stock) ~ "SAN_JUAN_RIVER",
      grepl("BOWRON", pbt_stock) ~ "BOWRON_RIVER",
      grepl("BRIDGE", pbt_stock) ~ "BRIDGE_RIVER",
      stock %in% c("CARIBOO_RIVER-UPPER", "U_CARIBOO") ~ 
        "CARIBOO_RIVER_UPPER",
      grepl("CHILAKO", pbt_stock) ~ "CHILAKO_RIVER",
      grepl("CHILCOTIN", pbt_stock) ~ "CHILCOTIN_RIVER_UPPER",
      grepl("CHILKO", pbt_stock) ~ "CHILKO_RIVER",
      TRUE ~ pbt_stock
    ),
    # some estimated ages don't line up with pbt brood; assume latter correct
    pbt_age = ifelse(!is.na(pbt_brood_year_n), year - pbt_brood_year_n, NaN),
    cwt_age = ifelse(!is.na(cwt_brood_year_n), year - cwt_brood_year_n, NaN),
    true_age = ifelse(!is.na(cwt_age), cwt_age, pbt_age),
    est_age = case_when(
      age_gr %in% c("2M", "3M", "4M", "1M", "0F", "44", "S1", "5M") ~ NA,
      is.na(age_gr) ~ NA,
      !is.na(age_gr) ~ substr(age_gr, 1, 1)
    ) %>% 
      as.numeric(.),
    age = ifelse(is.na(true_age), est_age, true_age), 
    sw_age = case_when(
      # pull first digit if M included in scale read age (indicating only 
      # saltwater age readable)
      grepl("M", age_gr) ~ map_dbl(age_gr, ~ str_split(.x, "(?<=\\d)(?=\\D)") %>%
                                     unlist() %>%
                                     .[[1]] %>%
                                     as.numeric()),
      # young 2.1s likely 1.2s
      (stock_group %in% c("FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2") |
         pst_agg %in% c("NBC_SEAK")) & age_gr == "21" ~ 1,
      stock_group %in% c("FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2") |
        pst_agg %in% c("NBC_SEAK") ~ age - 2,
      # if stock group has variable life history use identified yearlings
      (pst_agg %in% c("CR-upper_sp", "CR-upper_su/fa", "CR-lower_sp",
                      "CA_ORCST", "WACST","CR-lower_fa", "PSD") |
         stock_group %in% c("Fraser_Sum_4.1")) &
        age_gr %in% c("32", "42", "52", "62") ~ age - 2,
      TRUE ~ age - 1
    ),
    brood_year = year - age
  ) %>% 
  # left_join(
  #   ., pbt_rate, by = c("brood_year", "pbt_stock")
  # ) %>%
  mutate(
    #define hatchery status
    # origin = case_when(
    #   pbt == TRUE |
    #     resolved_stock_source %in% c("CWT", "Otolith Stock") ~ "hatchery",
    #   ad == "Y" ~ "hatchery",
    #   # if PBT rate for a stock is uniformly high and caught after ~2013 BY
    #   stock_group %in% c("FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2",
    #                      "FR_Sum_4.1", "FR_Fall", "ECVI_SOMN", "WCVI") &
    #     resolved_stock_source == "DNA" & year > 2016 & high == TRUE ~ "wild",
    #   # if PBT rate varied or caught before widespread adoption
    #   stock_group %in% c("FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2",
    #                      "FR_Sum_4.1", "FR_Fall", "ECVI_SOMN", "WCVI") &
    #     resolved_stock_source == "DNA" & tag_rate > 0.8 ~ "wild",
    #   # stock soruce not in PBT database
    #   stock_group %in% c("FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2",
    #                      "FR_Sum_4.1", "FR_Fall", "ECVI_SOMN", "WCVI") &
    #     resolved_stock_source == "DNA" & !(pbt_stock %in% pbt_rate$pbt_stock) ~
    #     "wild",
    #   grepl("HANFORD", stock) ~ "unknown",
    #   pst_agg %in% c("PSD", "WACST") | grepl("CR", pst_agg) & ad == "N" ~
    #     "wild",
    #   TRUE ~ "unknown"
    # ) %>%
    #   factor(., levels = c("hatchery", "unknown", "wild")),
    # nation = ifelse(
    #   pst_agg %in% c("FR-early", "FR-late", "SOG", "Yukon", "WCVI") |
    #     region4name == "NBC",
    #   "Can",
    #   "USA"
    # ),
    # origin2 = paste(origin, nation, sep = "_") %>%
    #   factor(.,
    #          levels = c("hatchery_Can",  "hatchery_USA", "unknown_Can",
    #                     "unknown_USA", "wild_Can", "wild_USA"),
    #          labels = c("Can Ha", "USA Ha", "Can Un", "USA Un", "Can Wi",
    #                     "USA Wi")),
    # adjust Capilano based on changes in brood source (ECVI pre 2013 brood);
    # assume fish caught in 2016 onward from Fraser Fall stock
    stock_group = ifelse(
      year > 2015 & grepl("CAPI", stock),
      "FR_Fall",
      as.character(stock_group)
      ) %>%
      factor(
        .,
        levels = c("other", "Col_Spr", "Col_Sum/Fall", "PSD", "WCVI",
                   "ECVI_SOMN", "FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2",
                   "FR_Sum_4.1", "FR_Fall")
      )
  )


# check for missing regional assignments
long_rec %>%
  filter(
    is.na(region1name)
  ) %>%
  select(
    stock
  ) %>%
  arrange(stock) %>%
  distinct()

## export
saveRDS(long_rec, here::here("data", "rec", "rec_gsi.rds"))



## CLEAN SIZE ------------------------------------------------------------------

wide_rec4 <- readRDS(here::here("data", "rec", "wide_rec.rds"))

wide_size <- wide_rec4 %>% 
  #remove missing and non-sensical size_classes
  filter(
    !is.na(fl)
  ) %>% 
  mutate(
    size_bin = cut(
      fl, 
      breaks = c(-Inf, 651, 751, 851, Inf), 
      labels = c("55-65", "65-75", "75-85", ">85")
    ),
    # add factor accounting for slot limits that went into place in different
    # years depending on whether west of 20-4/-5 line
    slot_limit = ifelse(
      ((lon < -124 | strata == "saanich") & year < 2019) | 
        ((lon < -124 | strata == "saanich")  & week_n > 28), "no", "yes" 
    )
  ) 


# visualize sample coverage through space and time
# size_n <- wide_size %>% 
#   group_by(month_n, year, strata) %>% 
#   tally() 

# ggplot(size_n) +
#   geom_raster(aes(x = month_n, y = year, fill = n)) +
#   facet_wrap(~strata)
# 
# wide_size %>%
#   group_by(size_bin, month_n, strata) %>%
#   summarize(n = length(unique(id))) %>%
#   group_by(month_n, strata) %>% 
#   mutate(total_n = sum(n),
#          ppn_obs = n / total_n) %>% 
#   ggplot(., aes(x = as.factor(month_n), y = ppn_obs, fill = size_bin)) +
#   geom_bar(position="stack", stat="identity") +
#   ggsidekick::theme_sleek() +
#   facet_wrap(~strata)
# 
# wide_size %>%
#   filter(strata %in% c("renfrew", "swiftsure", "swiftsure_nearshore")) %>% 
#   group_by(size_bin, month_n, strata, slot_limit) %>%
#   summarize(n = length(unique(id))) %>%
#   group_by(month_n, strata, slot_limit) %>% 
#   mutate(total_n = sum(n),
#          ppn_obs = n / total_n) %>% 
#   ggplot(., aes(x = as.factor(month_n), y = ppn_obs, fill = size_bin)) +
#   geom_bar(position="stack", stat="identity") +
#   ggsidekick::theme_sleek() +
#   facet_wrap(slot_limit~strata)

saveRDS(wide_size, here::here("data", "rec", "rec_size.rds"))


## CLEAN RKW DIET DATA ---------------------------------------------------------


raw_dat <- readRDS(
  here::here(
    "data", "rkw_diet", "raw_data", "RKW predation_chin samples_long_filtered.RDS"
  )
) %>% 
  filter(
    population == "SR"
  )


rkw_dat <- raw_dat %>% 
  # correct weird stock 
  mutate(
    stock = ifelse(stock == "BIGQUL@LANG", "BIG_QUALICUM", stock)
  ) %>% 
  left_join(., stock_key, by = "stock") %>%
  mutate(
    strata = ifelse(is.na(strata), "swiftsure", as.character(strata)) %>%
      factor(
        .,
        levels = c("swiftsure", "swiftsure_nearshore", "renfrew", "cJDF",
                   "sooke", "vic"),
        labels = c("Swiftsure", "Nitinat", "Renfrew", "cJDF",
                   "Sooke/\nVictoria", "San Juan\nIslands")
      ),
    month = lubridate::month(date),
    yday = lubridate::yday(date),
    # in-fill missing ages based on dominant life history strategy
    total_year = case_when(
      grepl("M", gr_age) & grepl(".2", smu) ~ sw_year + 2,
      grepl("M", gr_age) & !grepl(".2", smu) ~ sw_year + 1,
      TRUE ~ total_year
    ) %>% 
      as.factor(),
    era = ifelse(year < 2015, "early", "current") %>% 
      fct_relevel(., "current", after = Inf),
    # sampling event = all samples collected in a given strata-year-week
    sample_id = paste(year, week, strata, sep = "_"),
    sw_age = as.factor(sw_year),
    age_stock_group = case_when(
      grepl("Fraser", smu) ~ smu,
      stock == "CAPILANO" ~ "Fraser_Fall",
      agg == "SOG" ~ "ECVI_SOMN",
      TRUE ~ agg
    ) %>% 
      as.factor()
  ) %>% 
  # since weeks may span multiple months, calculate median month for each week
  group_by(sample_id) %>% 
  mutate(month = median(month)) %>% 
  ungroup() %>% 
  mutate(sample_id_pooled = paste(era, month, strata, sep = "_")) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("longitude", "latitude"), ll_crs = 4326, units = "km"
  ) %>% 
  select(
    id, sample_id, sample_id_pooled,
    fw_year, sw_age, total_year, age_f = total_year,
    era, year, month, week_n = week, yday, strata, utm_y = Y, utm_x = X,
    stock, stock_prob, stock_group, age_stock_group, stock_id_method,
    lat = latitude, lon = longitude
  )

saveRDS(rkw_dat, here::here("data", "rkw_diet", "cleaned_diet_samples.rds"))


## plot spatio temporal distribution of sampling events
rkw_dat %>% 
  group_by(year, yday) %>% 
  summarize(nn = length(unique(id))) %>% 
  ggplot(.) +
  geom_point(aes(y = as.factor(year), x = yday, size = nn))
