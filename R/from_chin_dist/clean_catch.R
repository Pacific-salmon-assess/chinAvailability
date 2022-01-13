## Clean catch/effort data from southcoast 
# Commercial data provided by BR Aug 19, 2019 as Access database which
# cannot be hosted publicly; imported csv attached here
# Recreational data provided by WL May 27, 2020

library(RODBC); library(tidyverse); library(ggplot2)

 
# ## CLEAN COMMERCIAL ------------------------------------------------------------
#
# # Access database saved locally to work computer
# file_path <- "C:/Users/FRESHWATERC/Documents/chinook/southCoastDatabase/southCoastGSICatch.accdb"
# con <- odbcConnectAccess2007(file_path)
# 
# # Check out table names
# sqlTables(con, tableType = "TABLE")$TABLE_NAME
# 
# ## Import various datatables
# # all gsi samples
# sampInvQry <- "SELECT *
#                FROM [Sample Inventory]"
# sampInv <- sqlQuery(con, sampInvQry) %>% 
#   select(-F14, -F15, -F16) %>% 
#   rename_all(list(~make.names(.))) #get rid of spaces in col names
# trimSampInv <- sampInv %>% 
#   select(Sample.ID, Vial.ID, Year, Fishery, Catch.Region, Catch.Area, 
#          Month.Code, Month.Name, Lab.Reported.N)
# 
# # stock composition data
# stockCompQry <- "SELECT *
#                  FROM [Stock Level Results] 
#                  INNER JOIN [Stocks with Region Codes]
#                  ON [Stock Level Results].[Stock Code] = 
#                   [Stocks with Region Codes].[Stock Code];"
# stockComp <- sqlQuery(con, stockCompQry) %>% 
#   rename_all(list(~make.names(.))) %>% 
#   left_join(., trimSampInv, by = "Sample.ID") %>% 
#   select(Sample.ID:SD, Region1Code:Lab.Reported.N) %>% 
#   rename(catchReg = Catch.Region, month = Month.Code, year = Year) %>% 
#   mutate(catchReg = as.character(catchReg)) %>% 
#   filter(Fishery == "Area G Com") 
# yrs <- unique(stockComp$year)
# 
# # confirm one sample collected per catch region, year, and month
# # dum <- stockComp %>%
# #   select(-c(Stock.Code:Vial.ID)) %>%
# #   distinct() %>%
# #   group_by(Catch.Region, Year, Month.Name) %>%
# #   tally()
# 
# # region key to reference stocks at different roll ups
# reg1Qry <- "SELECT *
#             FROM [Region 1 Stock Names];"
# reg1 <- sqlQuery(con, reg1Qry) %>% 
#   select(-ID)
# reg3Qry <- "SELECT *
#             FROM [Region 3 Stock Names];"
# reg3 <- sqlQuery(con, reg3Qry) %>% 
#   select(-ID)
# # frNameQry <- "SELECT *
# #               FROM [Fraser Stock Grouping Names]"
# reg2Qry <- "SELECT *
#            FROM [Stocks with Region Codes]
#            LEFT JOIN [Region 2 Stock Names]
#            ON [Stocks with Region Codes].Region2Code = 
#               [Region 2 Stock Names].Region2Code;"
# regKey <- sqlQuery(con, reg2Qry) %>%
#   rename_all(list(~make.names(.))) %>%
#   left_join(., reg1, by = "Region1Code") %>% 
#   left_join(., reg3, by = "Region3Code") %>% 
#   mutate(Stock = Stock.Name) %>% 
#   select(Stock.Code, Stock, Region2Name, Region1Name, Region3Name, 
#          Region1Code:FraserGroupCode)
# head(regKey)
# 
# trimRegKey <- regKey %>% 
#   mutate(Stock = as.character(Stock)) %>% 
#   select(Stock, Region1Name, Region2Name, Region3Name)
# 
# # write.csv(regKey, here::here("data", "southcoastStockKey.csv"), row.names = F)
# 
# # all catch/effort data
# catchQry <- "SELECT [Area G FOS Catch Estimates].ESTIMATE_TYPE, 
#                     [Area G FOS Catch Estimates].LICENCE_AREA, 
#                     [Area G FOS Catch Estimates].OPNG_CAT, 
#                     [Area G FOS Catch Estimates].OPNG_DESC, 
#                     [Area G FOS Catch Estimates].TARGETS_CHINOOK, 
#                     [Area G FOS Catch Estimates].STAT_WEEK, 
#                     [Area G FOS Catch Estimates].FISHING_DATE, 
#                     Month([FISHING_DATE]) AS [FISHING MONTH], 
#                     Year([FISHING_DATE]) AS [FISHING YEAR], 
#                     [Area G FOS Catch Estimates].MGMT_AREA, 
#                     [Area G FOS Catch Estimates].AREA_NAME, 
#                     [Management Regions].CATCH_REGION, 
#                     [Area G FOS Catch Estimates].HRS_OPEN, 
#                     [Area G FOS Catch Estimates].VESSELS_OP, 
#                     [Area G FOS Catch Estimates].CHINOOK_KEPT, 
#                     [Area G FOS Catch Estimates].CHINOOK_RELD, 
#                     [Area G FOS Catch Estimates].COMMENTS
#                FROM [Area G FOS Catch Estimates] 
#                INNER JOIN [Management Regions] 
#                ON [Area G FOS Catch Estimates].MGMT_AREA = 
#                   [Management Regions].[Managment Area]
#                WHERE ((([Area G FOS Catch Estimates].TARGETS_CHINOOK)='Yes'));
# "
# areaCatch <- sqlQuery(con, catchQry) %>% 
#   rename_all(list(~make.names(.))) %>%
#   filter(FISHING.YEAR %in% yrs) %>% 
#   # correct anomalously low vessel operating after looking at FOS
#   mutate(VESSELS_OP = 
#            case_when(
#              FISHING_DATE == "2013-05-15" & MGMT_AREA == "123" ~ 39,
#              TRUE ~ VESSELS_OP)
#          )
# write.csv(areaCatch, here::here("data", "gsiCatchData", "commTroll",
#                                 "fosCatch.csv"), 
#           row.names = FALSE)

#summary of how catch/effort data is distributed through time
# fosSumm <- areaCatch %>%
#   group_by(CATCH_REGION, FISHING.MONTH, FISHING.YEAR, MGMT_AREA) %>%
#   tally(name = "daysWithData")
# 
# write.csv(fosSumm, here::here("data", "gsiCatchData", "commTroll", 
#                               "fosSummary.csv"))

## Generate aggregate catch by Julian Day to match individual data from genetics
# lab
areaCatch <- read.csv(here::here("data", "gsiCatchData", "commTroll",
                                 "fosCatch.csv"))

dailyCatch <- areaCatch %>% 
  mutate(jDay = lubridate::yday(as.POSIXlt(FISHING_DATE, 
                                           format="%Y-%m-%d"))) %>% 
  group_by(CATCH_REGION, MGMT_AREA, FISHING.MONTH, jDay, FISHING.YEAR) %>% 
  #normally only KEPT chinook contribute to GSI samples, but certain
  #fisheries include live sampling where CHINOOK_RELD = catch
  summarize(catch =
              case_when(
                grepl("CN DNA Sampling", OPNG_DESC) ~ sum(CHINOOK_RELD),
                TRUE ~ sum(CHINOOK_KEPT)),
            boatDays = sum(VESSELS_OP)) %>% 
  rename(catchReg = CATCH_REGION, area = MGMT_AREA, year = FISHING.YEAR, 
         month = FISHING.MONTH) %>%
  ungroup() %>% 
  mutate(catchReg = as.character(catchReg),
         cpue = catch / boatDays) %>% 
  arrange(area, year, month, jDay)

saveRDS(dailyCatch, here::here("data", "gsiCatchData", "commTroll",
                                 "day_area_commCatch.rds"))


# CLEAN REC DATA ---------------------------------------------------------------

rec_catch <- read.csv(here::here("data", "gsiCatchData", "rec",
                                 "southcoast_rec_cpue_Feb2020.csv")) %>% 
  mutate(strata = paste(PFMA, YEAR, MONTH, DISPOSITION, sep = "_"))

#initial explore
# spatio-temporal data coverage
rec_catch %>% 
  select(PFMA, YEAR, MONTH, DISPOSITION) %>%
  distinct() %>% 
  ggplot() + 
  geom_bar(aes(x = MONTH, fill = DISPOSITION)) +
  ggsidekick::theme_sleek() +
  facet_wrap(~PFMA)

# estimate type breakdown
rec_catch %>% 
  group_by(STATUS) %>% 
  tally()

unpub_strata <- rec_catch %>%
  filter(STATUS != "Published Estimate - Full Month") %>%
  pull(strata)

# rec_catch %>%
#   filter(STATUS == "Published Estimate - Full Month",
#          strata %in% unpub_strata) %>% 
#   pull(strata)
#none seem to be cross referenced so retain both

# clean
rec_catch1 <- rec_catch %>% 
  mutate(
    PFMA = case_when(
      # separate northern areas of 13 (normally in JS) and add to NSoG; 
      # note that 13M = 13A and 13N = 13B (changed in 2011)
      CREEL_SUB_AREA %in% c("13M", "13N", "13A", "13B") ~ "PFMA 14",
      TRUE ~ PFMA
    ),
    pfma_n = as.numeric(str_remove_all(PFMA, "PFMA ")),
    region = case_when(
      pfma_n > 124 ~ "NWVI",
      pfma_n < 28 & pfma_n > 24 ~ "NWVI",
      pfma_n %in% c("20", "121", "21") ~ "Juan de Fuca Strait",
      is.na(pfma_n) ~ "Juan de Fuca Strait",
      pfma_n < 125 & pfma_n > 120 ~ "SWVI",
      pfma_n < 25 & pfma_n > 20 ~ "SWVI",
      pfma_n %in% c("14", "15", "16") ~ "N. Strait of Georgia",
      pfma_n %in% c("17", "18", "19", "28", "29") ~ "S. Strait of Georgia",
      pfma_n %in% c("10", "11", "111") ~ "Queen Charlotte Sound",
      pfma_n %in% c("12", "13") ~ "Queen Charlotte and\nJohnstone Straits"
      ),
    strata = paste(pfma_n, CREEL_SUB_AREA, YEAR, MONTH, sep = "_"),
    month = fct_relevel(MONTH, "January", "February", "March", "April", 
                        "May", "June", "July", "August", "September", 
                        "October", "November", "December"),
    month_n = as.numeric(month)
  ) %>%
  select(strata, month, month_n, year = YEAR, area = pfma_n,
         subarea = CREEL_SUB_AREA, region, 
         DISPOSITION, ADIPOSE_MARK, ESTIMATE, STANDARD_ERROR)

# separate effort and catch data, reformat, then recombine  
rec_eff <- rec_catch1 %>% 
  filter(DISPOSITION == "Effort") %>% 
  select(strata:region, mu_boat_trips = ESTIMATE, 
         se_boat_trips = STANDARD_ERROR)

rec_catch2 <- rec_catch1 %>% 
  filter(!DISPOSITION == "Effort") %>% 
  mutate(
    kept = case_when(
      DISPOSITION == "Kept" ~ "y",
      grepl("Released", DISPOSITION) ~ "n"
    ),
    legal = case_when(
      DISPOSITION == "Released Sub-Legal" ~ "sublegal",
      DISPOSITION %in% c("Released Legal", "Kept") ~ "legal"
      ),
    kept_legal = paste(kept, legal, sep = "_"),
    adipose_clip = case_when(
      ADIPOSE_MARK == "Adipose Marked" ~ "y",
      ADIPOSE_MARK == "Not Adipose Marked" ~ "n",
      TRUE ~ NA_character_)
    ) %>% 
  select(strata:region, legal, kept_legal, adipose_clip, mu_catch = ESTIMATE,
         se_catch = STANDARD_ERROR) %>% 
  distinct() %>% 
  filter(!is.na(mu_catch))

# rec_catch2 %>% 
#   group_by(strata, kept_legal, adipose_clip) %>% 
#   filter(n()>1)

# save list that includes uncertainty in both variables
rec_list <- list("catch" = rec_catch2, "effort" = rec_eff)

# combine into single dataframe that excludes uncertainty
rec_dat_out <- rec_catch2 %>% 
  select(-se_catch) %>%
  left_join(.,
            rec_eff %>% 
              select(strata, mu_boat_trips),
            by = "strata") %>% 
  select(-strata)

saveRDS(rec_list, here::here("data", "gsiCatchData", "rec",
                             "month_subarea_recCatch_list.rds"))
saveRDS(rec_dat_out, here::here("data", "gsiCatchData", "rec",
                             "month_subarea_recCatch.rds"))


# CLEAN DATA FOR MODELING ------------------------------------------------------

# clean function 
clean_catch <- function(dat) {
  dat %>% 
    filter(!is.na(eff),
           !eff == "0") %>% 
    mutate(region_c = as.character(region),
           region = factor(abbreviate(region, minlength = 4)),
           area_n = as.numeric(area),
           area = as.factor(area),
           month = as.factor(month),
           month_n = as.numeric(month),
           month = as.factor(month_n),
           year = as.factor(year),
           eff_z = as.numeric(scale(eff)),
           offset = log(eff)
    ) %>% 
    arrange(region, month) 
}

#commercial catch data
comm_catch_area_month <- readRDS(here::here("data", "gsiCatchData", "commTroll",
                                 "day_area_commCatch.rds")) %>% 
  # calculate monthly catch and effort 
  group_by(catchReg, area, month, year) %>% 
  summarize(catch = sum(catch), 
            eff = sum(boatDays),
            .groups = "drop") %>% 
  rename(region = catchReg) %>% 
  clean_catch(.) %>% 
  # drop inside areas where seasonal catches not available
  filter(!area_n < 100) %>% 
  droplevels() %>% 
  select(region, region_c, area, area_n, month, month_n, year, catch, eff, 
         eff_z, offset)
saveRDS(comm_catch_area_month, here::here("data", "gsiCatchData", "commTroll",
                             "month_area_commCatch.rds"))

#recreational catch data - sampling unit is area-month-year catch estimate
rec_catch_area_month <- readRDS(here::here("data", "gsiCatchData", "rec",
                                "month_subarea_recCatch.RDS")) %>% 
  #drop sublegal fish and regions without genetics
  filter(legal == "legal",
         !region %in% c("NWVI", "SWVI", "Queen Charlotte Sound")) %>% 
  #group by subarea to get rid of adipose and released legal duplicates
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  summarize(subarea_catch = sum(mu_catch, na.rm = T),
            subarea_eff = mean(mu_boat_trips, na.rm = T),
            .groups = "drop") %>% 
  filter(!is.na(subarea_eff)) %>% 
  #group by area to be consistent with commercial data
  group_by(month, month_n, year, area, region, legal) %>% 
  summarize(catch = sum(subarea_catch),
            eff = sum(subarea_eff),
            .groups = "drop") %>% 
  clean_catch(.) 
saveRDS(rec_catch_area_month, here::here("data", "gsiCatchData", "rec",
                                          "month_area_recCatch.rds"))

# export areas to make maps
areas_retained <- c(unique(rec_catch_area_month$area_n), 
                    unique(comm_catch_area_month$area_n))
saveRDS(areas_retained,
        here::here("data", "gsiCatchData", "pfma", "areas_to_plot.RDS"))


# REC C/E CHECK ----------------------------------------------------------------

# check breakdown among subareas to make sure no data gaps
c_list <- rec_dat_out %>% 
  #drop sublegal fish and regions without genetics
  filter(
    legal == "legal",
    !is.na(mu_boat_trips),
    # kept_legal == "y_legal",
    !region %in% c("NWVI", "SWVI", "Queen Charlotte Sound")) %>% 
  #group by subarea to get rid of adipose and released legal duplicates
  group_by(month, month_n, year, area, subarea, region, legal) %>% 
  summarize(subarea_catch = sum(mu_catch, na.rm = T),
            subarea_eff = mean(mu_boat_trips, na.rm = T),
            .groups = "drop") %>% 
  mutate(region_c = as.character(region),
         region = factor(abbreviate(region, minlength = 4)),
         min_m = case_when(
           region == "NSoG" ~ 5,
           region %in% c("QCaJS") ~ 6,
           region == "JdFS" ~ 6,
           region == "SSoG" ~ 5
         ),
         max_m = case_when(
           region %in% c("SSoG") ~ 9,
           region %in% c("JdFS", "NSoG", "QCaJS")  ~ 9
         ),
         date = zoo::as.yearmon(paste(.$year, .$month_n), "%Y %m"),
         subarea_cpue = subarea_catch / subarea_eff
  ) %>%
  group_by(area) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>%
  droplevels() %>% 
  split(., .$region)

pdf(here::here("figs", "hidden", "rec_sample_breakdown.pdf"))
map(c_list, function (x) {
  ggplot(x, aes(x = date, y = subarea, fill = subarea_cpue)) +
    geom_point(shape = 21) +
    scale_fill_viridis_c() +
    ggtitle(label = unique(x$region))
})
dev.off()

