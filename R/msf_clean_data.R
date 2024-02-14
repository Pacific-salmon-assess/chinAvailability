## Clean and compare MSF data to historical average
# Feb 14, 2024


library(tidyverse)

# stock key with aggs
stock_key <- readRDS(here::here("data", "rec", "finalStockList_Jan2024.rds")) %>%
  janitor::clean_names()

rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region)


site_ll <- rec_raw %>% 
  select(fishing_location = fishing_site, lat, lon) %>% 
  distinct()


## CLEAN -----------------------------------------------------------------------


# MSF data from 2024
msf_wide <- read_csv(
  here::here("data", "rec", "rec_msf_feb_24.csv"),
  na = c("","NA")
) %>% 
  janitor::clean_names(.) %>% 
  # add lat/lons
  left_join(
    ., site_ll, by = "fishing_location"
  ) %>% 
  mutate(
    date = as.POSIXct(collection_date, format="%d-%b-%Y"),
    month_n = lubridate::month(date),
    week_n = lubridate::week(date)
  )


# trim for GSI purposes (all DNA so no need to combine with CWT/otolith)
msf_wide_trim <- msf_wide %>% 
  # focus on taggable size (> 60 cm)
  filter(
    !is.na(resolved_stock_source),
    length_mm > 600
  ) %>%
  select(
    id = biokey, date, week_n, month_n, year, area, 
    fishing_site = fishing_location, subarea, lat, lon, 
    fl = length_mm, ad = adipose_fin_clipped,
    resolved_stock_source, 
    stock_1 = dna_results_stock_1, stock_2 = dna_stock_2, stock_3 = dna_stock_3,
    stock_4 = dna_stock_4, stock_5 = dna_stock_5,
    starts_with("prob")
  ) 



#pivot to long (probs and stock IDs separately) and join 
probs <- msf_wide_trim %>% 
  # replace with updated probabilities from above
  pivot_longer(., cols = starts_with("prob"), names_to = "rank", 
               names_pattern = "prob_(.+)",
               values_to = "prob") %>%
  select(id, rank, prob)

long_msf <- msf_wide_trim %>% 
  select(-starts_with("prob"), -starts_with("region")) %>% 
  pivot_longer(., cols = starts_with("stock"), names_to = "rank", 
               names_pattern = "stock_(.+)",
               values_to = "stock") %>% 
  left_join(., probs, by = c("id", "rank")) %>% 
  arrange(desc(date), id, desc(prob)) %>% 
  select(-rank) %>% 
  mutate(
    stock = toupper(stock)
  ) %>% 
  filter(!is.na(stock),
         !is.na(prob)) %>% 
  left_join(., stock_key, by = "stock") %>% 
  # add SMU-centric stock groups
  mutate(
    sample_id = paste(subarea, week_n, year, sep = "_"),
    stock_group = case_when(
      grepl("CR", pst_agg) ~ "Columbia",
      pst_agg %in% c("CA_ORCST", "WACST", "Russia", "NBC_SEAK") ~ "other",
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN") ~ "ECVI_SOMN",
      grepl("Fraser", region1name) ~ region1name,
     TRUE ~ pst_agg
    )
  )


comp_in <- long_msf  %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)) %>% as.numeric#,
    # # add abbreviated key for model fitting
    # stock_group2 = tolower(stock_group) %>% 
    #   str_replace_all(., " ", "_") %>% 
    #   paste("stock", ., sep = "-")
  ) %>% 
  ungroup() %>% 
  group_by(sample_id, 
           subarea, area, week_n, month_n, year, nn,
           stock_group) %>% 
  summarize(prob = sum(prob), .groups = "drop") 



## EXPLORE ---------------------------------------------------------------------

# sampling coverage
comp_in %>% 
  select(-c(stock_group, prob)) %>% 
  distinct() %>% 
  ggplot(.) +
  geom_jitter(aes(x = week_n, y = subarea, size = nn, colour = area),
              alpha = 0.4
  ) +
  facet_wrap(~ area, scales = "free_y") +
  ggsidekick::theme_sleek()


# size distribution
ggplot(msf_wide) + 
  geom_histogram(
    aes(x = length_mm)
  ) +
  facet_wrap(
    ~ area
  )


# number of samples by date in areas of interest
samps <- msf_wide %>% 
  filter(area %in% c(18, 19, 20),
         !subarea %in% c("19A", "18A"),
         length_mm > 59
         ) %>% 
  group_by(
    date, area 
  ) %>% 
  tally() %>% 
  print( 
    n = Inf
    )

comp_in %>% 
  mutate(
    agg_ppn = prob / nn
  ) %>% 
  ggplot(., aes(fill = stock_group, y = agg_ppn, x = week_n)) + 
  geom_bar(position="stack", stat="identity") +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  facet_wrap( ~ subarea) +
  labs(x = "Week", y = "Agg Probability") +
  ggsidekick::theme_sleek()
