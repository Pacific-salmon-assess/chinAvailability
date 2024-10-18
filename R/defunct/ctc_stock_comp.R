### Stock Specific Run Sizes
## Estimates of stock-specific age composition in WCVI AABM troll, WCVI ISBM 
# sport, and JdF ISBM sport fisheries.
# 1) Identify trends through time 
# 2) Generate average estimates for period with diet sample


library(tidyverse)

fishery_lookup <- read.csv(
  here::here("data", "ctc", "fisheryCodes.csv")
) %>% 
  janitor::clean_names()
stock_lookup <- read.csv(
  here::here("data", "ctc", "stockCodes.csv")
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
    fishery_short_name %in% c(
      "WCVI T", "GEO ST T", "WCVI AABM S", "WCVI ISBM S", "GEO ST S", "BC JF S"
    )
  )


# filter and calculate total catch, by age, by fisheries/stock
subset_catch <- rbind(isbm, aabm) %>% 
  filter(fishery %in% fishery1$fishery) %>% 
  left_join(., 
            fishery_lookup %>% select(fishery, fishery_short_name),
            by = "fishery") %>% 
  mutate(age = as.factor(age))
total_catch_at_age <- subset_catch %>% 
  select(starts_with("term") | starts_with("preterm")) %>% 
  apply(., 1, sum)
subset_catch$total_catch_at_age <- total_catch_at_age


# calculate total catch per fishery (regardless of age/stock)
total_catch <- subset_catch %>% 
  group_by(year, fishery, fishery_short_name) %>% 
  summarize(total_catch = sum(total_catch_at_age), .groups = "drop") %>% 
  ungroup()


# calculate total catch at age per fishery (pooling stocks)
subset_agg_catch <- subset_catch %>% 
  group_by(year, age, fishery, fishery_short_name) %>% 
  summarize(total_catch_at_age = sum(total_catch_at_age), .groups = "drop") %>% 
  ungroup() %>% 
  left_join(., total_catch, by = c("year", "fishery", "fishery_short_name")) %>% 
  mutate(ppn_catch_at_age = total_catch_at_age / total_catch) 


# age composition by fishery
age_comp_bar <- ggplot(subset_agg_catch, 
                       aes(fill = age, y = ppn_catch_at_age, x = year)) + 
  geom_bar(position="stack", stat="identity") +
  facet_wrap(~fishery_short_name) +
  ggsidekick::theme_sleek()


# visualize total catch by fisheries
fishery_catch <- ggplot(total_catch) +
  geom_jitter(aes(x = year, y =  fishery_short_name, size = total_catch),
              alpha = 0.5, height = 0.05) +
  ggsidekick::theme_sleek() 


# visualize average since 2014
mean_ppn_bar <- subset_agg_catch %>% 
  filter(year >= 2014, 
         !fishery_short_name == "GEO ST T") %>% 
  group_by(fishery, fishery_short_name, age) %>% 
  summarize(mean_ppn = mean(ppn_catch_at_age)) %>% 
  ggplot(., 
         aes(fill = age, y = mean_ppn, x = fishery_short_name)) + 
  geom_bar(position = "stack", stat = "identity") +
  labs(x =  "Fishery", y = "Proportion of Catch (since 2014)") +
  ggsidekick::theme_sleek()


pdf(here::here("figs", "fishery_age_comp", "age_comp_exp.pdf"))
fishery_catch
age_comp_bar
mean_ppn_bar
dev.off()

