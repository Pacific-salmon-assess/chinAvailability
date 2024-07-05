## BRMS Comparison
# Explore differences between stockseasonr and brms for estimating composition
# effects


library(tidyverse)
library(stockseasonr)
library(brms)
library(tidybayes)

options(mc.cores = parallel::detectCores())

source(here::here("R", "utils.R"))


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  # redefine region based on analysis
  mutate(
    cap_region = ifelse(
      lon > -124,
      "east",
      "west"
    ) %>% 
      factor(),
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
    stock_group2 = ifelse(
      stock_group2 %in% c("stock-psd", "stock-fraser_fall"),
      stock_group2,
      "stock-other"
      )
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
  cap_region = unique(comp_in$cap_region),
  month_n = seq(5, 10, length = 50)
)

## compare simple model fit with brms and dirichlet
fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + cap_region + 
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = comp_in,
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
  geom_pointrange(
    aes(x = cap_region, y = pred_prob_est, ymin = pred_prob_low, 
        ymax = pred_prob_up, colour = stock), 
    position = position_dodge(width = 0.35)
  )

ggplot(data = pred_fit$preds, 
           aes(x = month_n, colour = cap_region)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) + 
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = cap_region),
              alpha = 0.3) +
  geom_jitter(data = pred_fit$obs_dat,
              aes(x = month_n, y = obs_ppn, colour = cap_region, 
                  size = samp_nn),
              alpha = 0.2)



## brms fit

# extract and convert to matrix of proportions
comp_matrix <- fit$wide_comp_dat %>% 
  select(starts_with("stock-")) %>% 
  as.matrix()
comp_matrix2 <- (comp_matrix / apply(comp_matrix, 1, sum)) %>% 
  cbind()
comp_wide <- fit$wide_comp_dat %>% 
  select(-starts_with("stock-")) %>% 
  mutate(
    y = comp_matrix2
  )

brms_form <- bf(
  y ~ 1 + cap_region + 
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year)
    ,
  phi ~ nn,
  family = dirichlet
)

# explore priors
get_prior(brms_form, data = comp_wide)

## Set prior with sd = 1 for slopes and intercepts
priors_1 <- c(
  # prior(normal(0, 5), class = Intercept),
  prior(normal(0, 5), class = Intercept, dpar = "mustockother"),
  prior(normal(0, 5), class = Intercept, dpar = "mustockpsd"),
  prior(normal(0, 5), class = b, coef = "cap_regionwest", dpar = "mustockother"),
  prior(normal(0, 5), class = b, coef = "cap_regionwest", dpar = "mustockpsd"),
  prior(normal(0, 5), class = Intercept, dpar = "phi"),
  prior(normal(0, 5), class = b, coef = "nn", dpar = "phi")
)

## sample from the priors
ppp_1 <- brm(formula = brms_form, prior = priors_1,
             data = comp_wide, sample_prior = "only",
             chains = 1, cores = 1)

nd <- data.frame(
  nn = c(10, 100)
)

# plot (ineffective)
pred_1 <- predict(ppp_1, summary = F, ndraws = 100)
pp_draws <- prior_draws(ppp_1, ndraws = 100)
pp_draws <- posterior_predict(ppp_1, newdata = nd)
pp_draws <- simulate(ppp_1, newdata = nd)


pred_2 <- predict(ppp_2, summary = F, ndraws = 100)

# look at par ests
# mean(pred_1[ , , 1])
# sd(pred_2[ , , 1])
# mean(pred_2[ , , 2])
# mean(pred_1[ , , 2])
# sd(pred_1[ , , 3])



brms_fit <- brm(
  formula = brms_form, prior = priors_1,
  data = comp_wide,
  chains = 4, cores = parallel::detectCores() - 1,
  control = list(adapt_delta = 0.95)
)

# plot(marginal_effects(brms_fit, categorical = T, effects = "cap_region"),
#      plot = T)
# plot(marginal_effects(brms_fit2, categorical = T, effects = "cap_region"),
#      plot = T)
# 
# comp_in %>% 
#   mutate(pp = prob / nn) %>% 
#   group_by(cap_region, stock_group2) %>% 
#   summarize(
#     mean_ppn = mean(pp)
#   )

pred_dat1$year <- "2021" 
pred_dat1$nn <- 13 
brms_preds <- epred_draws(brms_fit, newdata = pred_dat1, re_formula = NA)
brms_preds2 <- brms_preds %>% 
  mutate(stock = str_remove(.category, "stock-")) %>% 
  group_by(cap_region, month_n, stock) %>% 
  summarize(
    pred_prob_est = median(.epred),
    pred_prob_low = quantile(.epred, 0.025),
    pred_prob_up = quantile(.epred, 0.975)
  ) %>% 
  mutate(
    model = "brms"
  )


ggplot(data = brms_preds2, 
       aes(x = month_n, colour = cap_region)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) + 
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = cap_region),
              alpha = 0.3) +
  geom_jitter(data = pred_fit$obs_dat,
              aes(x = month_n, y = obs_ppn, colour = cap_region, 
                  size = samp_nn),
              alpha = 0.2)

ss_preds <- pred_fit$preds %>% 
  mutate(
    model = "ss"
  ) %>% 
  select(colnames(brms_preds2)) 

dd <- rbind(ss_preds, brms_preds2)

ggplot(data = dd, 
       aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = model)) + 
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = model),
              alpha = 0.3) +
  geom_jitter(data = pred_fit$obs_dat,
              aes(x = month_n, y = obs_ppn, size = samp_nn),
              alpha = 0.2)
