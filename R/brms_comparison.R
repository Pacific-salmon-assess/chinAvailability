## BRMS Comparison
# Explore differences between stockseasonr and brms for estimating composition
# effects


library(tidyverse)
library(stockseasonr)
library(brms)
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
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) + 
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
  y ~ 1 #+ #cap_region + 
    # s(month_n, bs = "tp", k = 4, m = 2) +
    # (1 | year)
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
  prior(normal(0, 5), class = Intercept, dpar = "phi"),
  prior(normal(0, 5), class = b, coef = "nn", dpar = "phi")
)
# priors_2 <- c(
#   #prior(normal(0, 10), class = b, dpar = ""),
#   #prior(normal(0, 10), class = Intercept, dpar = ""),
#   prior(normal(0, 5), class = Intercept, dpar = ""),
#   prior(normal(0, 5), class = Intercept, dpar = "mustockother"),
#   prior(normal(0, 5), class = Intercept, dpar = "mustockpsd"),
#   prior(normal(0, 5), class = sd, group = "year", dpar = "mustockpsd"),
#   prior(normal(0, 5), class = sd, group = "year", dpar = "mustockother"),
#   prior(student_t(3, 0, 10), class = Intercept, dpar = "phi"),
#   prior(normal(0, 2), class = b, coef = "nn", dpar = "phi")
# )
## sample from the priors
ppp_1 <- brm(formula = brms_form, prior = priors_1,
             data = comp_wide, sample_prior = "only",
             chains = 1, cores = 1)
ppp_2 <- brm(formula = brms_form, prior = priors_2,
             data = comp_wide, sample_prior = "only",
             chains = 1, cores = 1)

pred_1 <- prior_draws(ppp_1, draws = 1000)
pred_2 <- prior_draws(ppp_2, draws = 1000)


nd <- data.frame(
  nn = c(10, 100)
)

# plot
pred_1 <- predict(ppp_1, summary = F, ndraws = 100)
pp_draws <- posterior_predict(ppp_1, newdata = nd)
pp_draws <- simulate(ppp_1, newdata = nd)


pred_2 <- predict(ppp_2, summary = F, ndraws = 100)


# plot(density(pred_1[1,,]), 
#      col = scales::alpha("blue", 0.2), 
#      main = "b ~ Normal(0,1)",
#      xlim = c(-0.1, 1.1))
# for(i in 2:100) lines(density(pred_1[i,,]), 
#                       col = scales::alpha("blue", 0.2))

# look at par ests
# mean(pred_1[ , , 1])
# sd(pred_2[ , , 1])
# mean(pred_2[ , , 2])
# mean(pred_1[ , , 2])
# sd(pred_1[ , , 3])



brms_fit <- brm(
  formula = brms_form, prior = priors_2,
  data = comp_wide,
  chains = 4, cores = parallel::detectCores() - 1
)

# brms_fit2 <- brm(
#   formula = brms_form2, prior = priors_2, 
#   data = comp_wide,
#   chains = 4, cores = parallel::detectCores() - 1
# )

plot(marginal_effects(brms_fit, categorical = T, effects = "cap_region"),
     plot = T)
plot(marginal_effects(brms_fit2, categorical = T, effects = "cap_region"),
     plot = T)

comp_in %>% 
  mutate(pp = prob / nn) %>% 
  group_by(cap_region, stock_group2) %>% 
  summarize(
    mean_ppn = mean(pp)
  )
