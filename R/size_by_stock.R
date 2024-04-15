## Size-Age Models
# Use GSI samples through 2022 to evaluate trends in size-at-age by stock
# Oct 30, 2023


library(tidyverse)


gsi <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  mutate(
    #add management measure effect accounting for slot limits
    mm = ifelse(year %in% c("2019", "2020", "2021", "2022"), "yes", "no"),
    stock_group = case_when(
      stock == "Capilano" | region1name %in% c("ECVI", "SOMN") ~ "ECVI_SOMN",
      grepl("Fraser", region1name) ~ region1name,
      TRUE ~ pst_agg
    )
  ) %>%
  group_by(
    id, strata, stock_group, week_n, month_n,  mm, fl, age, year
  ) %>% 
  summarize(
    sum_prob = sum(prob), .groups = "drop"
  ) %>% 
  ungroup() %>% 
  # focus on: study area; Canadian stocks + Puget Sound + Columbia; May-Oct
  # remove uncertain assignments
  filter(
    sum_prob < 0.75,
    !strata == "other"
  ) %>% 
  mutate(
    strata = factor(
      strata,
      levels = c("swiftsure", "swiftsure_nearshore", "renfrew", "vic",
                 "haro", "saanich"),
      labels = c("Swiftsure", "Nitinat", "Renfrew", "Sooke/\nVictoria",
                 "S. Gulf\nIslands", "Saanich")
    ),
    stock_group = factor(
      stock_group,
      levels = c(
        "CA_ORCST", "CR-upper_sp", "CR-upper_su/fa", "CR-lower_sp",
        "CR-lower_fa", "WACST", "PSD", "WCVI", "ECVI_SOMN",
        "Fraser_Spring_4.2", "Fraser_Spring_5.2", "Fraser_Summer_5.2",
        "Fraser_Summer_4.1", "Fraser_Fall",  "NBC_SEAK"
      )
    ),
    year_f = as.factor(year)
  )
  

## AGE COMPOSITION -------------------------------------------------------------

age_comp <- gsi %>%
  filter(!is.na(age_gr)) %>% 
  group_by(stock_group) %>%
  mutate(age_n = n()) %>%
  ungroup() %>%
  group_by(stock_group, age_gr, age_n) %>%
  tally() %>% 
  mutate(prop = n / age_n)

labs_age_comp <- age_comp %>%
  ungroup() %>% 
  select(stock_group, age_n) %>% 
  distinct()

age_comp_stacked <- ggplot() +
  geom_bar(data = age_comp,
           aes(fill = as.factor(age_gr), y = prop, x = stock_group),
           position="stack", stat="identity") +
  geom_text(data = labs_age_comp, 
            aes(x = stock_group, y = 0.05, label = age_n)) +
  scale_fill_brewer(name = "Age", palette = "Paired", na.value = "grey60" ) +
  
  # scale_fill_viridis_d(name = "Age", na.value = "grey60" ) +
  labs(y = "Proportion Age Composition", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

png(
  here::here("figs", "ms_figs", "comp_bar_fishery_age.png"),
  height = 5, width = 8, units = "in", res = 250
)
age_comp_stacked
dev.off()


## BODY SIZE -------------------------------------------------------------------

size_at_age_box <- gsi %>%
  filter(!is.na(age)#,
         # !age %in% c("2", "6")
         ) %>% 
  ggplot(.) +
  geom_boxplot(aes(x = stock_group, y = fl)) +
  labs(y = "Fork Length", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  ) +
  facet_wrap(~age)

ggplot(gsi %>% filter(fl > 550)) +
  geom_boxplot(aes(x = stock_group, y = fl, fill = mm))  +
  labs(y = "Fork Length", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )



library(mgcv)

fit2 <- gam(
  fl ~ 0 + age + s(week_n, bs = "tp", k = 4, m = 2) +
    s(week_n, bs = "tp", by = stock_group, k = 4, m = 1) +
    mm:stock_group + s(year_f, bs = "re"),
  data = gsi
)


# make predictions, constraining to weeks where stocks present
stock_week <- gsi %>% 
  filter(month_n > 4 & month_n < 11) %>% 
  group_by(week_n, stock_group) %>% 
  tally() 

ggplot(stock_week,
       aes(x = week_n, y = stock_group, size = n, fill = stock_group)) +
  geom_point(shape = 21) +
  scale_size_continuous(trans = "log") +
  ggsidekick::theme_sleek()

obs_weeks <- stock_week %>% 
  group_by(stock_group) %>% 
  summarize(max_obs_week = max(week_n),
            min_obs_week = min(week_n))


# restrict to weeks where at least two individuals sampled
week_month <- data.frame(
  week_n = c(20, 25, 29, 34, 37, 42),
  month = c("May", "Jun", "Jul", "Aug", "Sep", "Oct")
)
new_dat <- expand.grid(
  week_n = unique(week_month$week_n),
  stock_group = unique(gsi$stock_group),
  age = levels(gsi$age)
) %>% 
  left_join(., week_month, by = "week_n") %>% 
  left_join(., obs_weeks, by = "stock_group") %>% 
  # remove weeks where stock not observed
  filter(
    !week_n > max_obs_week,
    !week_n < min_obs_week
  ) %>% 
  mutate(year_f = "2020",
         mm = "no",
         month = fct_reorder(month, week_n))

preds <- predict(fit2, newdata = new_dat,  se.fit = TRUE, 
                 exclude = "s(year_f)")

new_dat2 <- new_dat %>% 
  mutate(
    pred_fl = as.numeric(preds$fit),
    lo_fl = pred_fl + (stats::qnorm(0.025) * as.numeric(preds$se.fit)), 
    up_fl = pred_fl + (stats::qnorm(0.975) * as.numeric(preds$se.fit))
    )

size_month2 <- ggplot(new_dat2) +
  geom_pointrange(
    aes(x = month, y = pred_fl, ymin = lo_fl, ymax = up_fl, 
        fill = age),
    shape = 21
  ) +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )


png(
  here::here("figs", "ms_figs", "temporal_sample_coverage.png"),
  height = 5, width = 7.5, units = "in", res = 250
)



pdf(here::here("figs", "stock_size_age", "summary_figs.pdf"))
age_comp_stacked
size_at_age_box
size_month1
size_month2
dev.off()