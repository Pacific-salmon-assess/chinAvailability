## Size-Age Models
# Use GSI samples through 2022 to evaluate trends in size-at-age by stock
# Oct 30, 2023


library(tidyverse)


gsi <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  filter(
    !fl < 551
  ) %>% 
  mutate(
    age_stock_group1 = case_when(
      grepl("CAPI", stock) ~ stock_group,
      stock_group == "PSD" ~ stock_group,
      grepl("Fraser", region1name) ~ stock_group,
      TRUE ~ pst_agg
    ),
    age_stock_group = ifelse(
      age_stock_group1 == "ECVI_SOMN", "SOG", age_stock_group1
      )
  ) %>%
  # calculate the stock group level probability for each individual
  group_by(
    id, strata, stock_group, age_stock_group, week_n, month_n, fl,
    age, sw_age, year
  ) %>% 
  summarize(
    sum_prob = sum(prob), .groups = "drop"
  ) %>% 
  ungroup() %>% 
  # focus on: study area
  # remove uncertain assignments
  filter(
    sum_prob > 0.75,
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
    age_stock_group = factor(
      age_stock_group,
      levels = c(
        "CA_ORCST", "CR-upper_sp", "CR-upper_su/fa", "CR-lower_sp",
        "CR-lower_fa", "WACST", "PSD", "WCVI", "SOG",
        "FR_Spr_4.2", "FR_Spr_5.2", "FR_Sum_5.2",
        "FR_Sum_4.1", "FR_Fall",  "NBC_SEAK"
      )
    ),
    year_f = as.factor(year),
    age_f = as.factor(age),
    sw_age = as.factor(sw_age)
  ) 


# colour palettes
age_pal <- c("grey30", "#a6cee3", "#1f78b4", "#b2df8a", "#33a02c", "black")
names(age_pal) <- c(NA, levels(gsi$sw_age))


## AGE COMPOSITION -------------------------------------------------------------

age_comp <- gsi %>%
  filter(!is.na(sw_age)) %>% 
  group_by(stock_group) %>%
  mutate(age_n = n()) %>%
  ungroup() %>%
  group_by(stock_group, sw_age, age_n) %>%
  tally() %>% 
  mutate(prop = n / age_n)

labs_age_comp <- age_comp %>%
  ungroup() %>% 
  select(stock_group, age_n) %>% 
  distinct()

age_comp_stacked <- ggplot() +
  geom_bar(data = age_comp,
           aes(fill = sw_age, y = prop, x = stock_group),
           position="stack", stat="identity", colour = "black") +
  geom_text(data = labs_age_comp, 
            aes(x = stock_group, y = 0.05, label = age_n)) +
  scale_fill_manual(name = "Marine\nAge", values = age_pal, na.value = "grey60" ) +
  labs(y = "Proportion Age Composition", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )

png(
  here::here("figs", "stock_size_age", "comp_bar_fishery_age_sw.png"),
  height = 5, width = 8, units = "in", res = 250
)
age_comp_stacked
dev.off()


## GAM FIT ---------------------------------------------------------------------

library(mgcv)
library(gratia)


fit <- gam(
  fl ~ #0 + 
    sw_age + age_stock_group + s(week_n, bs = "tp", k = 4, m = 2) +
    s(week_n, bs = "tp", by = sw_age, k = 4, m = 1) +
    s(week_n, bs = "tp", by = age_stock_group, k = 4, m = 1) +
    s(year_f, bs = "re")
  ,
  data = gsi %>% filter(!is.na(sw_age))
)

# make predictions, constraining to weeks where stocks present
stock_week <- gsi %>% 
  filter(month_n > 4 & month_n < 10) %>% 
  group_by(week_n, age_stock_group) %>% 
  tally() 

obs_weeks <- stock_week %>% 
  group_by(age_stock_group) %>% 
  summarize(max_obs_week = max(week_n),
            min_obs_week = min(week_n))

stock_age <- gsi %>% 
  group_by(sw_age, age_stock_group) %>% 
  summarize(
    age_n = length(unique(id))
  ) 

# restrict to weeks where at least five individuals sampled
week_month <- data.frame(
  week_n = c(20, 25, 29, 34, 38, 42),
  month = c("May", "Jun", "Jul", "Aug", "Sep", "Oct"),
  month_n = c(5, 6, 7, 8, 9, 10)
)
new_dat1 <- expand.grid(
  week_n = unique(week_month$week_n),
  age_stock_group = unique(gsi$age_stock_group),
  sw_age = unique(gsi$sw_age) %>% as.factor()
) %>% 
  left_join(., week_month, by = "week_n") %>% 
  left_join(., obs_weeks, by = "age_stock_group") %>% 
  left_join(., stock_age, by = c("age_stock_group", "sw_age")) %>% 
  filter(
    !is.na(sw_age)
  )
new_dat1$year_f <- unique(gsi$year_f)[1]
new_dat1$sw_year <- as.numeric(as.character(new_dat1$sw_age))
  
new_dat <- new_dat1 %>%   
  # remove weeks where stock not observed & rare age classes
  filter(
    !week_n > max_obs_week,
    !week_n < min_obs_week,
    !age_n < 5,
    !is.na(sw_age)
  ) %>% 
  mutate(
    month = fct_reorder(month, week_n)
  )


# generate predictions, integrating out year effects
preds <- predict(fit, newdata = new_dat, se.fit = TRUE, exclude = "s(year_f)")

new_dat2 <- new_dat %>% 
  mutate(
    pred_fl = as.numeric(preds$fit),
    lo_fl = pred_fl + (stats::qnorm(0.025) * as.numeric(preds$se.fit)), 
    up_fl = pred_fl + (stats::qnorm(0.975) * as.numeric(preds$se.fit)),
    age_stock_group = fct_recode(
      age_stock_group, 
      "FR_Spr_4sub2" = "FR_Spr_4.2",
      "FR_Spr_5sub2" = "FR_Spr_5.2",
      "FR_Sum_5sub2" = "FR_Sum_5.2",
      "FR_Sum_4sub1" = "FR_Sum_4.1"
    )
    )

shape_pal <- c(21, 22)
names(shape_pal) <- c("yes", "no")

size_month2 <- ggplot(new_dat2) +
  geom_pointrange(
    aes(x = month, y = pred_fl, ymin = lo_fl, ymax = up_fl, 
        fill = sw_age
        ),
    shape = 21
  ) +
  scale_fill_manual(name = "Marine\nAge", values = age_pal, 
                    na.value = "grey60" ) +
  facet_wrap(~age_stock_group) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    y = "Predicted Fork Length"
  )

png(
  here::here("figs", "stock_size_age", "mean_pred_fishery.png"),
  height = 6.5, width = 8, units = "in", res = 250
)
size_month2
dev.off()


# generate 1000 simulated draws for each combination 
# subset based on data necessary for diet predictions 
# (generated in clean_data.R)
key <- readRDS(here::here("data", "size_at_age_pred_key.rds")) %>% 
  mutate(
    id = paste(age_stock_group, sw_year, month, sep = "_")
  )
new_dat_post_out <- new_dat1 %>% 
  mutate(
    id = paste(age_stock_group, sw_year, month_n, sep = "_")
  ) %>% 
  filter(
    id %in% key$id
  )

set.seed(123)
sims <- simulate(fit, nsim = 1000, data = new_dat_post_out, exclude = "s(year_f)")

size_pred_post <- sims %>% 
  as.data.frame() %>%
  cbind(new_dat_post_out %>% 
          select(month, sw_age, age_stock_group),
        .) %>% 
  pivot_longer(cols = starts_with("V"), names_to = "iter", names_prefix = "V",
               values_to = "fl") %>% 
  mutate(iter = as.numeric(iter))

saveRDS(size_pred_post, here::here("data", "rec", "size_age_post_draws.rds"))


# histogram of stock-specific size distributions
size_pred_hist <- ggplot(
  size_pred_post %>% 
    mutate(
      age_stock_group = fct_recode(
        age_stock_group, 
        "FR_Spr_4sub2" = "FR_Spr_4.2",
        "FR_Spr_5sub2" = "FR_Spr_5.2",
        "FR_Sum_5sub2" = "FR_Sum_5.2",
        "FR_Sum_4sub1" = "FR_Sum_4.1"
      )
    ) %>% 
    filter(month == "Jul")
) +
  geom_histogram(aes(x = fl, fill = sw_age)) +
  scale_fill_manual(name = "Marine\nAge", values = age_pal, 
                    na.value = "grey60" ) +
  facet_wrap(~age_stock_group) +
  ggsidekick::theme_sleek() +
  theme(
    axis.title.x = element_blank(),
    legend.position = "top"
  ) +
  labs(
    y = "Predicted Fork Length"
  )

png(
  here::here("figs", "stock_size_age", "hist_pred_fishery.png"),
  height = 6.5, width = 8, units = "in", res = 250
)
size_pred_hist
dev.off()

# ggplot(gsi %>% filter(!is.na(sw_age))) +
#   geom_histogram(aes(x = fl, fill = sw_age)) +
#   scale_fill_manual(name = "Marine\nAge", values = age_pal, 
#                     na.value = "grey60" ) +
#   facet_wrap(~age_stock_group) +
#   ggsidekick::theme_sleek() +
#   theme(
#     axis.title.x = element_blank(),
#     legend.position = "top"
#   ) +
#   labs(
#     y = "Predicted Fork Length"
#   )
