## Size-Age Models
# Use GSI samples through 2022 to evaluate trends in size-at-age by stock
# Oct 30, 2023


library(tidyverse)


gsi <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  # focus on:
  # study area 
  # Canadian stocks + Puget Sound + Columbia
  # May-Oct
  filter(
    !cap_region == "outside",
    grepl("CR", pst_agg) | grepl("Fraser", stock_group) |
      grepl("Puget", region1name) | pst_agg %in% c("WCVI", "SOG"),
    month_n > 4 & month_n < 11) %>% 
  mutate(
    stock_group = case_when(
      grepl("Fraser", region1name) ~ region1name,
      TRUE ~ pst_agg
    ),
    
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
  geom_boxplot(aes(x = stock_group, y = fl))  +
  labs(y = "Fork Length", x = "Stock") +
  ggsidekick::theme_sleek() +
  theme(
    axis.text.x = element_text(angle = 45, hjust=1)
  )


## significant differences among stocks?
size_trim <- gsi %>% 
  filter(fl > 550) %>% 
  mutate(year_f = as.factor(year))


library(lme4)

fit <- lmer(
  data = size_trim,
  fl ~ stock_group + (1 | year_f)
 )
