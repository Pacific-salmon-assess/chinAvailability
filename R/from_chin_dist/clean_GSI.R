## Clean individual genetic data - Ver 2
# Updated version of gsiIndProbsClean.R that uses new flat file with full
# probabilities (i.e. not trimmed to 5 stocks max) for WCVI commercial troll 
# fishery and trimmed probabilities for JDF, SoG, and JS rec fisheries
# May 26, 2020

library(tidyverse)


## COMMERCIAL GSI --------------------------------------------------------------

dat_raw_comm <- read.csv(here::here("data", "gsiCatchData", "commTroll", 
                                    "hidden", "wcviIndProbsLong_RAW.txt"), 
                   stringsAsFactors = FALSE)

# Big chunk of code to separate ID variable into meaningful individual vectors
id_vec <- dat_raw_comm$szLineInfo %>% 
  as.vector() %>% 
  strsplit(., split = " ") %>% 
  unlist() %>% 
  matrix(., nrow = 5, ncol = length(dat_raw_comm$szLineInfo)) %>%
  t() %>% 
  data.frame() %>% 
  rename("statArea" = X1, "year" = X2, "gear" = X3, "jDay" = X4,
         "fishNum" = X5) %>% 
  mutate(
    jDay = as.numeric(as.character(jDay)),
    statArea = 
      case_when(
         statArea %in% c("Area023", "Area23", "Area_23") ~ "23",
         statArea %in% c("Area123", "Area123SWVI", "Area123Comm") ~ "123",
         statArea %in% c("Area124", "Area124SWVI", "Area124Comm", "Area123-124",
                         "Area124_24") ~ "124",
         statArea %in% c("Area125", "Area125NWVI") ~ "125",
         statArea %in% c("Area126", "Area126NWVI", "Area125-126", 
                         "Area126-127") ~ "126",
         statArea %in% c("Area127", "Area127NWVI") ~ "127",
         statArea %in% c("Area026") ~ "26",
         statArea %in% c("Area24", "Area_24", "Area_24xgill") ~ "24",
         TRUE ~ as.character(statArea)
         ),
    abbYear = sapply(strsplit(as.character(year), '[()]'), 
                function(x) (x)[2]),
    year = paste("20", abbYear, sep = ""),
    #adjust sampling day to correct for errors by genetics lab and 
    jDay = 
      case_when(
        statArea == "126" & year == "2012" & jDay > 48 &
          jDay < 116 ~ 48,
        statArea == "126" & year == "2012" & jDay > 116 &
          jDay < 121 ~ 116,
        TRUE ~ jDay),
    date = as.Date(as.numeric(as.character(jDay)),
                   origin = as.Date(paste(year, "01", "01", sep = "-"))),
    month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
    week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d")),
    weekDay = weekdays(date),
    fishNum = as.numeric(as.character(fishNum))
    ) %>% 
  # adjust sampling week
  # for specific strata based on when samples were landed relative to julian 
  # date of when fishing occurred (based on reviewing FOS database and pers. 
  # comm. Lee Kearey - South Coast)
  mutate(
    adjWeek = case_when(
      statArea == "23" & month == "3" & year == "2013" ~ week - 1,
      statArea == "24" & month == "5" & year == "2013" ~ week - 1,
      statArea == "26" & month == "2" & year == "2" ~ week - 1,
      statArea == "123" & month == "6" & year == "2007" ~ week - 1,
      statArea == "123" & month == "6" & year == "2008" ~ week - 1,
      statArea == "125" & month == "6" & year == "2007" ~ week - 1,
      statArea == "125" & month == "10" & year == "2007" ~ week - 1,
      statArea == "125" & month == "5" & year == "2014"~ week - 1,
      statArea == "126" & month == "3" & year == "2007"~ week - 1,
      statArea == "126" & month == "2" & year == "2012"~ week - 1,
      statArea == "126" & month == "3" & year == "2013"~ week - 1,
      statArea == "126" & month == "4" & year == "2013"~ week - 1,
      statArea == "126" & month == "3" & year == "2015"~ week - 1,
      statArea == "127" & month == "6" & year == "2007"~ week - 1,
      statArea == "127" & month == "4" & year == "2014"~ week - 1,
      statArea == "23" & month == "6" & year == "2007"~ week - 1,
      # shift back to account for catch being harvested 3-4 days before landing
      weekDay %in% c("Sunday", "Monday") ~ week - 1,
      TRUE ~ week
      )
    ) %>% 
  #shuffle adjusted
  rename(week = adjWeek, unadjWeek = week)

#Merge id vector with original data frame and trim
dat_comm <- cbind(id_vec, dat_raw_comm) %>% 
  #calculate total summed probability for each sample
  group_by(szLineInfo) %>%
  mutate(stock = toupper(szStock),
         totalProb = sum(dProb)) %>% 
  ungroup() %>% 
  mutate(adj_prob = dProb / totalProb) %>% 
  select(-iRun, -iSample, -iYearMix, id = szLineInfo, stock,  -szStock,
         prob = dProb, adj_prob, -totalProb, -szExclude, -iRegionId, -unadjWeek,
         Region1Name = szRegion, -weekDay) %>% 
  rename(area = statArea, fish_num = fishNum) 


## Export list of stocks to be passed to makeFullStockKey script in
# stockKey repo 
# stks_out <- dat_comm %>%
#   select(stock, Region1Name) %>%
#   distinct()
# saveRDS(stks_out, here::here("data", "stockKeys", "comm_gsi_stocks.rds"))

# stock key generated in stockKey repo
stockKey <- readRDS(here::here("data", "stockKeys", "finalStockList_Nov2020.rds"))

comm_long <- dat_comm %>%
  select(-Region1Name) %>%
  left_join(., stockKey, by = c("stock")) %>% 
  arrange(area, year, month, jDay, fish_num) %>% 
  mutate(
    season_c = case_when(
      month %in% c("12", "1", "2") ~ "w",
      month %in% c("3", "4", "5") ~ "sp",
      month %in% c("6", "7", "8") ~ "su",
      month %in% c("9", "10", "11") ~ "f"
    ),
    season = fct_relevel(season_c, "sp", "su", "f", "w"),
    month_n = as.numeric(month),
    month = as.factor(month_n),
    year =  as.factor(year),
    pres = 1,
    area_n = as.numeric(as.character(area)),
    region = case_when(
      area_n < 125 & area_n > 27 ~ "SWVI",
      area_n < 25 ~ "SWVI",
      TRUE ~ "NWVI"
    ),
    region = as.factor(region), 
    area = as.factor(area),
    temp_strata = paste(month_n, region, sep = "_")
    ) %>% 
  select(id, fish_num, temp_strata, region, area, year, month, week, jDay, date,
         gear, pres, season, month_n, area_n, adj_prob, stock, 
         Region1Name:pst_agg) %>% 
  arrange(year, region, id, desc(adj_prob))

saveRDS(comm_long, here::here("data", "gsiCatchData", "commTroll", 
                              "wcviIndProbsLong.rds"))
# comm_long <- readRDS(here::here("data", "gsiCatchData", "commTroll",
#                            "wcviIndProbsLong.rds"))

## Clean recreation GSI data ---------------------------------------------------

rec_full <- read.csv(here::here("data", "gsiCatchData", "rec", 
                                "hidden", "rec_gsi_may2020.txt"), 
                     stringsAsFactors = F)

# pull stocks to add to stockkey repo
# stk_out <- rec_full %>%
#   filter(!DNA_RESULTS_STOCK_1 == "",
#          !is.na(PROB_1)) %>%
#   select(s1 = DNA_RESULTS_STOCK_1, s2 = DNA_STOCK_2, s3 = DNA_STOCK_3,
#          s4 = DNA_STOCK_4, s5 = DNA_STOCK_5, REGION_1_ROLLUP) %>%
#   pivot_longer(., cols = s1:s5, names_to = "rank", values_to = "stock") %>%
#   select(stock, sc_reg1 = REGION_1_ROLLUP) %>%
#   filter(!stock == "") %>%
#   distinct()
# saveRDS(stk_out, here::here("data", "stockKeys", "rec_gsi_stocks.rds"))


stockKey <- readRDS(here::here("data", "stockKeys", "finalStockList_Nov2020.rds"))

# data frame of probabilities
temp_prob <- rec_full %>% 
  filter(!DNA_RESULTS_STOCK_1 == "",
         RESOLVED_STOCK_SOURCE == "DNA") %>% 
  rename(p1 = PROB_1, p2 = PROB_2, p3 = PROB_3, p4 = PROB_4, p5 = PROB_5) %>% 
  pivot_longer(., cols = c(p1, p2, p3, p4, p5),
               names_to = "rank_prob", values_to = "prob") %>%  
  select(BIOKEY, COLLECTION_DATE, rank_prob, prob)

# rec_full %>% 
#   filter(PFMA == "13") %>% 
#   group_by(SUBAREA) %>% 
#   tally()

# data frame of stock IDs
rec_long <- rec_full %>% 
  filter(!DNA_RESULTS_STOCK_1 == "",
         RESOLVED_STOCK_SOURCE == "DNA") %>% 
  rename(s1 = DNA_RESULTS_STOCK_1, s2 = DNA_STOCK_2, s3 = DNA_STOCK_3, 
         s4 = DNA_STOCK_4, s5 = DNA_STOCK_5,
         p1 = PROB_1, p2 = PROB_2, p3 = PROB_3, p4 = PROB_4, p5 = PROB_5) %>% 
  pivot_longer(., cols = c(s1, s2, s3, s4, s5), 
               names_to = "rank", values_to = "stock") %>% 
  cbind(., temp_prob %>% select(rank_prob, prob)) %>% 
  #add regional roll ups
  left_join(., stockKey, by = "stock") %>% 
  filter(!is.na(prob)) %>%
  mutate(date = as.Date(as.numeric(as.character(DAYOFYEAR - 1)),
                        origin = as.Date(paste(YEAR, "01", "01", sep = "-"))),
         month = lubridate::month(as.POSIXlt(date, format="%Y-%m-%d")),
         week = lubridate::week(as.POSIXlt(date, format="%Y-%m-%d"))
  ) %>% 
  #adjust probabilities 
  group_by(BIOKEY) %>% 
  mutate(total_prob = sum(prob),
         adj_prob = prob / total_prob) %>% 
  ungroup() %>% 
  mutate(
    region = case_when(
      # separate northern areas of 13 (normally in JS) and add to NSoG
      SUBAREA %in% c("13M", "13N") ~ "N. Strait of Georgia",
      PFMA > 124 ~ "NWVI",
      PFMA < 28 & PFMA > 24 ~ "NWVI",
      PFMA %in% c("20", "121", "21") ~ "Juan de Fuca Strait",
      is.na(PFMA) ~ "Juan de Fuca Strait",
      PFMA < 125 & PFMA > 120 ~ "SWVI",
      PFMA < 25 & PFMA > 20 ~ "SWVI",
      PFMA %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      PFMA %in% c("17", "18", "19", "28", "29") ~ "S. Strait of Georgia",
      PFMA %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      PFMA %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
    ),
    season_c = case_when(
      month %in% c("12", "1", "2") ~ "w",
      month %in% c("3", "4", "5") ~ "sp",
      month %in% c("6", "7", "8") ~ "su",
      month %in% c("9", "10", "11") ~ "f"
    ),
    legal_lim = case_when(
      PFMA < 20 & PFMA > 11 ~ 620,
      PFMA %in% c("28, 29") ~ 620,
      TRUE ~ 450
    ),
    legal = case_when(
      LENGTH_MM >= legal_lim ~ "legal",
      LENGTH_MM < legal_lim ~ "sublegal",
      KEPTREL == "Kept" ~ "legal",
      KEPTREL == "Rel" ~ "sublegal"
    ),
    season = fct_relevel(season_c, "sp", "su", "f", "w"),
    month_n = as.numeric(month),
    month = as.factor(month_n),
    year =  as.factor(YEAR),
    pres = 1,
    area_n = as.numeric(as.character(PFMA)),
    region = as.factor(region), 
    area = case_when(
      is.na(PFMA) ~ "20_121_21",
      TRUE ~ as.character(PFMA) 
    ),
    area =  as.factor(area),
    temp_strata = paste(month_n, region, sep = "_"),
    gear = "sport",
    ad_clip = case_when(
      ADIPOSE_FIN_CLIPPED == "Y" ~ "Y",
      TRUE ~ "N"
    )
  ) %>% 
  select(id = BIOKEY, fish_num = FISH_NO, temp_strata, region, area, 
         subarea = SUBAREA, 
         year, month, week, jDay = DAYOFYEAR, date, gear = gear, 
         fl = LENGTH_MM, release = KEPTREL, legal, sex = SEX, 
         ad_clip, pres, season, sampler = SAMPLER_TYPE,
         month_n, area_n, adj_prob, stock, Region1Name:pst_agg) %>% 
  arrange(year, region, id, desc(adj_prob))

saveRDS(rec_long, here::here("data", "gsiCatchData", "rec",
                            "recIndProbsLong.rds"))
