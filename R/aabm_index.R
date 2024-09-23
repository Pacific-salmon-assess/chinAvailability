## Estimates of AABM abundance
# Sep 23, 2024

library(tidyverse)

aabm_dat <- read.csv(
  here::here("data", "run_size", "aabm_fishery_abundance_indices.csv"),
  stringsAsFactors = FALSE
) %>% 
  group_by(aabm_fishery) %>% 
  mutate(
    rolling_mean = zoo::rollmean(abundance_index, k = 5, fill = NA, align = "right")
  )



ggplot(aabm_dat) +
  geom_point(aes(x = year, y = abundance_index, fill = aabm_fishery),
             shape = 21) +
  geom_line(aes(x = year, y = rolling_mean, colour = aabm_fishery)) +
  facet_wrap(~aabm_fishery) +
  ggsidekick::theme_sleek()
