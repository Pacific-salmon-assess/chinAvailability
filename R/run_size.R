### Stock Specific Run Sizes
## Generate estimates of seasonal abundance by MU using estimates of run size
## and seasonal run timing (i.e. time series of annual run timing curves)
## Initially just for Fraser, but ultimately extent to WCVI, ECVI, Puget Sound
## and Columbia River


library(tidyverse)

# import and edit population names based on stock_reference_lsit.xlsx from CP
fr_return <- read.csv(here::here("data", "run_size", "fraser_run_size.csv")) %>% 
  janitor::clean_names() %>% 
  pivot_longer(cols = -(year),
               names_to = "population", 
               values_to = "return") %>% 
  mutate(
    pop_new = case_when(
      population == "adams" ~ "lower_adams",
      population == "fraser" ~ "fraser_tj",
      population %in% c("harrison", "chilliwack", "stave") ~ "harrison",
      TRUE ~ population
    )
  )


fr_timing <- read.csv(
  here::here("data", "run_size", "fraser_ck_run_timing.csv")) %>% 
  mutate(
    pop_new1 = tolower(population) %>%
      str_split(., " ") %>% 
      purrr::map(., ~ .x[1]) %>% 
      as.character() %>% 
      c(),
    # change stock group to match grouping
    stock_group = if_else(pop_new1 == "duteau", "spring_42", stock_group),
    pop_new = case_when(
      pop_new1 == "thompson" ~ "lower_thompson",
      pop_new1 %in% c("baezaeko", "west_road", "nazko") ~ "blackwater",
      pop_new1 == "duteau" ~ "bessette",
      pop_new1 %in% c("fontoniko", "herrick", "james") ~ "mc_gregor",
      pop_new1 == "indianpoint" ~ "bowron",
      TRUE ~ pop_new1
    )
  ) 

# make key of common stocks for Chuck
timing_pops <- fr_timing %>% 
  pull(pop_new) %>% 
  unique() %>% sort()

return_pops <- fr_return %>% 
  pull(pop_new) %>% 
  unique() %>% sort()

pop_data_summary <- data.frame(
  population = c(timing_pops, return_pops) %>% unique()
) %>% 
  mutate(
    timing_data = case_when(
      # population %in% c("alouette", "chilliwack", "stave", "coquitlam") ~ "yes", 
      population %in% timing_pops ~ "yes",
      TRUE ~ "no"
      ),
    run_size_data = ifelse(population %in% return_pops, "yes", "no")
  ) %>% 
  arrange(timing_data, run_size_data, population) %>% 
  filter(!population == "total")

write.csv(pop_data_summary, here::here("data", "run_size", "stock_summary.csv"),
          row.names = FALSE)


# POOLED ESTIMATES -------------------------------------------------------------

# generate pooled estimates for stocks groups with multiple run timing estimates
stock_groups <- fr_timing %>% 
  group_by(pop_new) %>% 
  tally() %>% 
  filter(n > 1) %>% 
  pull(pop_new)


group_tbl <- fr_timing %>% 
  filter(pop_new %in% stock_groups) %>% 
  select(mu = stock_group, pop_new, mean_yday_gsi, n_gsi, sd_gsi) %>% 
  group_by(mu, pop_new) %>% 
  group_nest() %>% 
  mutate(
    mean = NA,
    sd = NA,
    n = NA
  )

# use 100 draws from combined normal distribution to generate shared distribution
# weighted by observed sample size
ndraw <- 100
for (i in 1:nrow(group_tbl)) {
  dum <- group_tbl$data[[i]]
  
  dist_mat <- matrix(NA, nrow = sum(dum$n_gsi), ncol = ndraw)
  for(k in seq_len(ndraw)) {
    dist_list <- vector(mode = "list", length = nrow(dum))
    for (j in 1:nrow(dum)) {
      dist_list[[j]] <- rnorm(
        n = dum$n_gsi[j], mean = dum$mean_yday_gsi[j],sd = dum$sd_gsi[j]
      ) %>% 
        round(., digits = 0)
    }
    dist_mat[ , k] <- do.call(c, dist_list)
  }
  group_tbl$mean[i] <- mean(dist_mat)
  group_tbl$sd[i] <- sd(dist_mat)
  group_tbl$n[i] <- sum(dum$n_gsi)
}

# merge into original
fr_timing2 <- fr_timing %>%
  select(mu = stock_group, pop_new, mean = mean_yday_gsi, sd = sd_gsi, 
         n = n_gsi) %>% 
  filter(!pop_new %in% group_tbl$pop_new) %>% 
  rbind(., group_tbl %>% select(-data))


# CALCULATE ANNUAL RUN DISTRIBUTIONS -------------------------------------------

focal_stocks <- pop_data_summary %>% 
  filter(timing_data == "yes" & run_size_data == "yes",
         # only one run timing sample
         !population == "mahood") 

ret_tbl <- fr_return %>% 
  select(-population) %>% 
  filter(pop_new %in% focal_stocks$population) %>% 
  left_join(., 
            fr_timing2 %>% 
              select(mu, pop_new, mean, sd),
            by = "pop_new") %>% 
  as_tibble()

ret_list <- vector(mode = "list", length = nrow(ret_tbl))
for (i in seq_along(ret_list)) {
  ret_list[[i]] <- data.frame(
    year = ret_tbl$year[i],
    population = ret_tbl$pop_new[i],
    mu = ret_tbl$mu[i],
    yday = 75:300
  ) %>% 
    mutate(
      run = dnorm(yday, mean = ret_tbl$mean[i], sd = ret_tbl$sd[i]) *
        ret_tbl$return[i]
    ) 
}

annual_ret <- ret_list %>%
  bind_rows() %>% 
  filter(year %in% c("1979", "1980")) 

annual_ret %>% 
  group_by(mu, yday, year) %>% 
  summarize(mu_run = sum(run), .groups = "drop") %>% 
  ggplot(.) +
  geom_line(aes(x = yday, y = mu_run)) +
  facet_grid(mu~year, scales = "free_y")
