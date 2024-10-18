## Import and clean Fraser escapement
# Dec. 15, 2023

library(tidyverse)


## FRASER ----------------------------------------------------------------------

dat_raw <- read_csv(here::here("data", "run_size", "Escape1979-22_ADJ.csv")) %>% 
  janitor::clean_names() %>% 
  .[-1, ] 

fraser_dat <- dat_raw %>% 
  pivot_longer(cols = starts_with("x"), names_prefix = "x", names_to = "year",
               values_to = "esc") %>% 
  filter(!is.na(stock_name),
         !stock_name %in% c("New Stocks Total", "Total")) %>% 
  mutate(
    year = as.numeric(year),
    management_period = case_when(
      year < 2014 ~ "historical",
      year >= 2014 & year <= 2018  ~ "base",
      year > 2018 ~ "current"
    )
  ) %>% 
  group_by(agg_name, year, management_period) %>% 
  summarize(
    smu_esc = sum(esc)
    ) %>% 
  ungroup()


## ECVI ------------------------------------------------------------------------

ecvi_dat <- read.csv(
  here::here(
    "data", "run_size", "ecvi_esc.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  rename(year = x) %>% 
  pivot_longer(cols = -year, names_to = "stock", values_to = "esc") %>% 
  # remove data prior to 2004
  filter(year > 2003) %>% 
  mutate(
    agg_name = "ECVI",
    management_period = case_when(
      year < 2014 ~ "historical",
      year >= 2014 & year <= 2018  ~ "base",
      year > 2018 ~ "current"
    ) 
  ) %>% 
  group_by(stock, agg_name, management_period) %>% 
  mutate(
    infill_esc = ifelse(is.na(esc), mean(esc, na.rm = T), esc)
  ) %>% 
  group_by(agg_name, year, management_period) %>% 
  summarize(
    smu_esc = sum(infill_esc, na.rm = T)
  ) %>% 
  ungroup()


## WCVI ------------------------------------------------------------------------

wcvi_dat <- read.csv(
  here::here(
    "data", "run_size", "wcvi_esc.csv"
  )
) %>% 
  janitor::clean_names() %>% 
  pivot_longer(cols = starts_with("x"), names_to = "year", values_to = "esc",
               names_prefix = "x") %>% 
  # same yaer as Fraser data
  filter(year > 1979) %>% 
  mutate(
    year = as.numeric(year),
    esc = as.numeric(esc),
    na_flag = ifelse(is.na(esc), 1, 0),
    agg_name = "WCVI",
    management_period = case_when(
      year < 2014 ~ "historical",
      year >= 2014 & year <= 2018  ~ "base",
      year > 2018 ~ "current"
    ) 
  ) %>% 
  group_by(stream_name) %>% 
  mutate(n_na = sum(na_flag)) %>% 
  filter(!n_na > 10) %>% 
  group_by(stream_name, agg_name, management_period) %>% 
  mutate(
    infill_esc = ifelse(is.na(esc), mean(esc, na.rm = T), esc)
  ) %>% 
  group_by(agg_name, year, management_period) %>% 
  summarize(
    smu_esc = sum(infill_esc, na.rm = T)
  ) %>% 
  ungroup()


## POOL DATA -------------------------------------------------------------------

mean_esc <- list(fraser_dat, wcvi_dat, ecvi_dat) %>%
  bind_rows() %>% 
  mutate(
    region = case_when(
      agg_name == "ECVI" ~ "ECVI",
      agg_name == "WCVI" ~ "WCVI",
      TRUE ~ "Fraser"
    ),
    agg_name = as.factor(agg_name),
    agg_name = fct_relevel(agg_name, "Fall", after = Inf) %>% 
      fct_reorder(., as.numeric(as.factor(region)))
  ) %>% 
  group_by(agg_name, region) %>% 
  summarize(
    mean_esc = mean(smu_esc),
    sd_esc = sd(smu_esc)
  )

smu_dat_period <- list(fraser_dat, wcvi_dat, ecvi_dat) %>%
  bind_rows() %>% 
  mutate(
    region = case_when(
      agg_name == "ECVI" ~ "ECVI",
      agg_name == "WCVI" ~ "WCVI",
      TRUE ~ "Fraser"
    )
  ) %>% 
  filter(!management_period == "historical") %>% 
  group_by(management_period, agg_name, region) %>%
  summarize(
    mean_esc = mean(smu_esc),
    up = mean_esc + (1.96 + sd(smu_esc)),
    lo = mean_esc - (1.96 + sd(smu_esc)),
    agg_name = as.factor(agg_name),
    agg_name = fct_relevel(agg_name, "ECVI", "WCVI", "Spring 4.2", 
                           "Spring 5.2", "Summer 5.2", "Summer 4.1", "Fall")
  ) %>%
  distinct()

png(here::here("figs", "avg_escapement.png"),
    height = 4.5, width = 5, units = "in", res = 250)
ggplot() +
  geom_pointrange(
    data = smu_dat_period %>% filter(!agg_name == "WCVI"),
    aes(x = management_period, y = mean_esc, ymax = up, ymin = lo, 
        fill = region),
    shape = 21) +
  geom_hline(data = mean_esc %>% filter(!agg_name == "WCVI"), 
             aes(yintercept = mean_esc), colour = "red") +
  ggsidekick::theme_sleek() +
  facet_wrap(~agg_name, scales = "free_y") +
  labs(x = "Management Period", y = "Mean Escapement") +
  scale_fill_discrete(guide = "none")
dev.off()