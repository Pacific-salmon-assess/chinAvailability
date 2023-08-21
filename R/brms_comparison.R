## BRMS Comparison
# Explore differences between stockseasonr and brms for estimating composition
# effects


library(tidyverse)
library(stockseasonr)
library(brms)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  # redefine region based on analysis
  mutate(
    cap_region = case_when(
      lat < 48.8 & lon > -125.25 & lon < -124.25 ~ "swiftsure",
      lat < 48.45 & lon < -123.4 & lon > -124.25 ~ "sooke",
      TRUE ~ "outside"
    ),
    whale_samples_time = ifelse(
      (year < 2011 | year > 2017) & month_n %in% c("6", "7", "8") & 
        rkw_habitat == "yes",
      "yes",
      "no"
    ),
    year = as.factor(year),
    yday = lubridate::yday(date),
    day_samp = paste(yday, year)
  )

comp_in <- rec_raw %>% 
  filter(
    !legal == "sublegal",
    !cap_region == "outside"
  ) %>% 
  mutate(
    sample_id = paste(cap_region, week_n,  #rkw_habitat,
                      year, sep = "_")
  ) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)) %>% as.numeric,
    # add abbreviated key for model fitting
    stock_group2 = tolower(stock_group) %>% 
      str_replace_all(., " ", "_") %>% 
      paste("stock", ., sep = "-"),
    # restrict to two stocks
    stock_group2 = ifelse(stock_group2 == "stock-psd", "stock-psd", "stock-other")
  ) %>% 
  ungroup() %>% 
  group_by(sample_id, 
           # rkw_habitat, 
           cap_region, week_n, month_n, year, nn, 
           stock_group2) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  filter(
    month_n >= 5 & month_n <= 10
  )


pred_dat1 <- expand.grid(
  cap_region = unique(comp_in$cap_region)
)

## compare simple model fit with brms and dirichlet
fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + cap_region,
  comp_dat = comp_in_trim,
  pred_dat = pred_dat1,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  # nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)

pred_fit <- clean_pred_foo(fit = fit, preds = pred_dat1)

ggplot(pred_fit$preds) +
  geom_pointrange(aes(x = cap_region, y = pred_prob_est, ymin = pred_prob_low, 
                      ymax = pred_prob_up, colour = stock), 
                  position = position_jitter(width = 0.3))



## brms fit

# extract and convert to matrix of proportions
comp_matrix <- fit$wide_comp_dat %>% 
  select(starts_with("stock-")) %>% 
  as.matrix()
comp_matrix2 <- (comp_matrix / apply(comp_matrix, 1, sum)) %>% 
  cbind()
comp_wide <- comp_wide %>% 
  select(-starts_with("stock-")) %>% 
  mutate(
    y = comp_matrix2
  )

brms_form <- bf(
  y ~ 1 + cap_region,
  family = dirichlet
)

# explore priors
get_prior(brms_form, data = comp_wide)

## Set prior with sd = 1
priors_1 <- c(prior(normal(0,1), class = b),
              prior(normal(0,1), class = Intercept),
              prior(student_t(3, 0, 10), class = phi))
## sample from the priors
ppp_1 <- brm(formula = brms_form, prior = priors_1, 
             data = comp_wide, sample_prior = "only",
             chains = 1, cores = 1)
# plot
pred_1 <- predict(ppp_1, summary = F, ndraws = 100)
plot(density(pred_1[1,,]), 
     col = scales::alpha("blue", 0.2), 
     main = "b ~ Normal(0,1)",
     xlim = c(-0.1, 1.1))
for(i in 2:100) lines(density(pred_1[i,,]), 
                      col = scales::alpha("blue", 0.2))

brms_fit <- brm(
  formula = brms_form, #prior = priors_1, 
  data = comp_wide,
  chains = 4, cores = parallel::detectCores() - 1
)

plot(marginal_effects(brms_fit, categorical = T, effects = "cap_region"),
     plot = T)