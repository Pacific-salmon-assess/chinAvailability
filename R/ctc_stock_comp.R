### Stock Specific Run Sizes
## Estimates of stock-specific age composition in WCVI AABM troll, WCVI ISBM 
# sport, and JdF ISBM sport fisheries.
# 1) Identify trends through time 
# 2) Generate average estimates for period with diet sample


library(tidyverse)

fishery_lookup <- read.csv(
  here::here("data", "ctc", "2022ERA_FisheryLookup.csv")
) %>% 
  janitor::clean_names()
stock_dat <- read.csv(
  here::here("data", "ctc", "CWTDataforMapping_CTC_big.csv")
) %>% 
  janitor::clean_names()
isbm <- read.csv(
  here::here("data", "ctc", "2304P_fish_ISBM_CCC.csv")
) %>% 
  janitor::clean_names()
aabm <- read.csv(
  here::here("data", "ctc", "2304P_fish_AABM_CCC.csv")
) %>% 
  janitor::clean_names()


## identify relevant fisheries
fishery1 <- fishery_lookup %>% 
  filter(
    fishery_name %in% c("WCVI F/W T", "WCVI SPR T", "WCVI SUM T", "WCVI AABM S",
                        "WCVI ISBM S", "BC JF S", "TBC JF TERM S")
  )

# filter and calculate total catch, by age, across types 
subset_catch <- rbind(isbm, aabm) %>% 
  filter(fishery %in% fishery1$fishery_number) %>% 
  left_join(., 
            fishery_lookup %>% select(fishery = fishery_number, fishery_name),
            by = "fishery") 
sum_catch <- subset_catch %>% 
  select(starts_with("term") | starts_with("preterm")) %>% 
  apply(., 1, sum)
subset_catch$sum_catch <- sum_catch

subset_catch  %>% 
  group_by(year, fishery, fishery_name) %>% 
  summarize(sum_catch = sum(total_catch)) %>% 
  ggplot() +
  geom_jitter(aes(x = year, y =  fishery_name, size = sum_catch),
              alpha = 0.5, width = 0.25) +
  ggsidekick::theme_sleek() 

