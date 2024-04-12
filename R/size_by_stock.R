## Size-Age Models
# Use GSI samples through 2022 to evaluate trends in size-at-age by stock
# Oct 30, 2023


library(tidyverse)


gsi <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  mutate(
    stock_group = case_when(
      grepl("Fraser", region1name) ~ region1name,
      TRUE ~ pst_agg
    ),
    #add management measure effect accounting for slot limits
    mm = ifelse(year %in% c("2019", "2020", "2021", "2022"), "yes", "no")
  ) %>% 
  group_by(
    id, stock_group, week_n, month_n,  mm, fl, age, year
  ) %>% 
  summarize(
    sum_prob = sum(prob), .groups = "drop"
  ) %>% 
  ungroup() %>% 
  # focus on: study area; Canadian stocks + Puget Sound + Columbia; May-Oct
  # remove uncertain assignments
  filter(
    sum_prob < 0.75,
    # !cap_region == "outside",
    month_n > 4 & month_n < 11,
    stock_group %in% c(
      "CR-upper_sp", "CR-lower_sp", "CR-upper_su/fa",
      "CR-lower_fa", "WCVI", "Fraser_Spring_5.2", "Fraser_Spring_4.2",
      "Fraser_Summer_5.2",  "Fraser_Summer_4.1", "Fraser_Fall", "SOG", "PSD"
    )
  ) %>% 
  mutate(
    stock_group = fct_relevel(
      as.factor(stock_group), "CR-upper_sp", "CR-lower_sp", "CR-upper_su/fa",
      "CR-lower_fa", "WCVI", "Fraser_Spring_5.2", "Fraser_Spring_4.2",
      "Fraser_Summer_5.2",  "Fraser_Summer_4.1", "Fraser_Fall", "SOG", "PSD"
    )
  )
  

## AGE COMPOSITION -------------------------------------------------------------

age_comp <- gsi %>%
  filter(!is.na(age)) %>% 
  group_by(stock_group) %>%
  mutate(age_n = n()) %>%
  ungroup() %>%
  group_by(stock_group, age, age_n) %>%
  tally() %>% 
  mutate(prop = n / age_n,
         age = as.factor(age))

labs_age_comp <- age_comp %>%
  ungroup() %>% 
  select(stock_group, age_n) %>% 
  distinct()

age_comp_stacked <- ggplot() +
  geom_bar(data = age_comp,
           aes(fill = age, y = prop, x = stock_group),
           position="stack", stat="identity") +
  geom_text(data = labs_age_comp, 
            aes(x = stock_group, y = 0.05, label = age_n)) +
  scale_fill_viridis_d(name = "Age", na.value = "grey60" ) +
  labs(y = "Proportion Age Composition", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )


## BODY SIZE -------------------------------------------------------------------

size_at_age_box <- gsi %>%
  filter(!is.na(age),
         !age %in% c("2", "6")) %>% 
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


## significant differences among stocks?
size_trim <- gsi %>% 
  filter(fl > 550) %>% 
  mutate(year_f = as.factor(year))


library(mgcv)

fit2 <- gam(
  fl ~ s(week_n, bs = "tp", k = 4, m = 2) +
    s(week_n, bs = "tp", by = stock_group, k = 4, m = 1) +
    mm:stock_group + s(year_f, bs = "re"),
  data = size_trim
)

# make predictions, constraining to weeks where stocks present
stock_week <- gsi %>% 
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
  stock_group = unique(size_trim$stock_group)
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

size_month1 <- ggplot(new_dat2) +
  geom_pointrange(
    aes(x = stock_group, y = pred_fl, ymin = lo_fl, ymax = up_fl, 
        fill = stock_group),
    shape = 21
  ) +
  facet_wrap(~month) +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_blank()
  )

size_month2 <- ggplot(new_dat2) +
  geom_pointrange(
    aes(x = month, y = pred_fl, ymin = lo_fl, ymax = up_fl, 
        fill = stock_group),
    shape = 21
  ) +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank()
  )


pdf(here::here("figs", "stock_size_age", "summary_figs.pdf"))
age_comp_stacked
size_at_age_box
size_month1
size_month2
dev.off()