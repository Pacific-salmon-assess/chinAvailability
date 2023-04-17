### Stock Specific Run Sizes
## Generate estimates of seasonal abundance by MU using estimates of run size
## and seasonal run timing (i.e. time series of annual run timing curves)
## Initially just for Fraser, but ultimately extent to WCVI, ECVI, Puget Sound
## and Columbia River


library(tidyverse)


fr_return <- read.csv(here::here("data", "run_size", "fraser_run_size.csv")) %>% 
  janitor::clean_names() %>% 
  pivot_longer(cols = -(year), names_to = "population", values_to = "return") 


fr_timing <- read.csv(
  here::here("data", "run_size", "fraser_ck_run_timing.csv")) %>% 
  mutate(pop_new = tolower(population) %>%
           str_split(., " ") %>% 
           purrr::map(., ~ .x[1]) %>% 
           as.character() %>% 
           c()) 

timing_pops <- fr_timing %>% 
  pull(pop_new) %>% 
  unique()

return_pops <- fr_return %>% 
  pull(population) %>% 
  unique()

pop_data_summary <- data.frame(
  population = c(timing_pops, return_pops) %>% unique()
) %>% 
  mutate(
    timing_data = case_when(
      population %in% c("alouette", "chilliwack", "stave", "coquitlam") ~ "yes", 
      population %in% timing_pops ~ "yes",
      TRUE ~ "no"
      ),
    run_size_data = ifelse(population %in% return_pops, "yes", "no")
  ) %>% 
  arrange(timing_data, run_size_data, population) %>% 
  filter(!population == "total")

write.csv(pop_data_summary, here::here("data", "run_size", "stock_summary.csv"),
          row.names = FALSE)


# calculate density 
focal_stocks <- pop_data_summary %>% 
  filter(timing_data == "yes" & run_size_data == "yes")

ret_tbl <- fr_return %>% 
  filter(population %in% focal_stocks$population) %>% 
  left_join(., 
            fr_timing %>% 
              select(population = pop_new, mean_yday = mean_yday_gsi, sd_gsi),
            by = "population") %>% 
  as_tibble()

ret_list <- vector(mode = "list", length = nrow(ret_tbl))
for (i in seq_along(ret_list)) {
  ret_list[[i]] <- data.frame(
    year = ret_tbl$year[i],
    population = ret_tbl$population[i],
    yday = 75:300
  ) %>% 
    mutate(
      run = dnorm(yday, mean = ret_tbl$mean_yday[i], sd = ret_tbl$sd_gsi[i]) *
        ret_tbl$return[i]
    ) 
}

