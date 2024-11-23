# Calculate Age Classification Errror
# Uses data exported from clean_data.R
# Nov 22, 2024
#1) make age_gr vs age key 
#2) assign scale_age based on age_gr and resolved_age
#3) define true age based on cwt, or pbt
#4) subset to samples that have both a true age and scale age estimate
#5) calculate classification accuracy

library(tidyverse)

age_dat <- readRDS(here::here("data", "rec", "aging_data.rds"))

# make age_key 
age_key <- data.frame(
  age_gr = c("21", "31", "32", "33", "41", "42", "51", "52", "53", "62"),
  est_age = c(2, 3, 3, 3, 4, 4, 5, 5, 5, 6)
)

dum <- age_dat %>% 
  left_join(., age_key, by = "age_gr") %>% 
  mutate(
    bias = case_when(
      age - est_age == 0 ~ "zero",
      age - est_age > 0 ~ "under",
      age - est_age < 0 ~ "over"
    )
  )

samp_size <- dum %>% 
  group_by(stock_group, age) %>% 
  summarize(
    n = length(unique(biokey))
  )

ppn_dat <- dum %>% 
  group_by(stock_group, age) %>% 
  mutate(age_n = n(), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(stock_group, bias, age, age_n) %>% 
  tally() %>% 
  mutate(prop = n / age_n)

ggplot(ppn_dat) +
  geom_bar(aes(x = age, y = prop, fill = bias),
           position="stack", stat="identity") +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  geom_text(
    data = samp_size, aes(x = age, y = -Inf, label = paste(n)),
    vjust = -1.1
  )


