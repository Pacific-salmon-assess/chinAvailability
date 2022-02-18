## Fit stock-comp model to generate area-level predictions assuming 
# homogeneous composition within a region
# Feb 17, 2022

library(tidyverse)


## SIMULATE --------------------------------------------------------------------
## simulate data to check indexing 
dum_c <- data.frame(
  strata = c("x", "y", "z"),
  stock_a = c(0.3, 0.2, 0.7)
) %>% 
  mutate(
    stock_b = 1 - a
  )
dum_a <- expand.grid(
  strata = dum_c$strata,
  substrata = c(1, 2, 3)
) %>% 
  mutate(substrata = paste(strata, substrata, sep = "_")) %>% 
  arrange(strata)
dum_a$catch <- exp(rnorm(nrow(dum_a), 0, 1))
key <- as.numeric(dum_a$strata)

dum <- matrix(NA, nrow = nrow(dum_a), ncol = 2)
for (i in 1:nrow(dum)) {
  for(j in 1:ncol(dum)) {
    dum[i, j] <- dum_a$catch[i] * dum_c[key[i], j + 1]
  }
}


## FIT TO JDF ------------------------------------------------------------------


# tmb models
compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal") %>% 
  mutate(reg = abbreviate(cap_region, minlength = 4),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         sample_id = paste(month_n, reg, yday, year, sep = "_")) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)),
    coarse_agg = ifelse(pst_agg %in% c("PSD", "CR-lower_fa", "FR-late", 
                                       "FR-early", "SOG"),
                        "dom", "weak")
  ) %>% 
  ungroup()
stock_comp <- comp1 %>%  
  group_by(sample_id, reg, reg_c = cap_region, month, month_n, year, nn, 
           coarse_agg) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(reg %in% c("JdFS", "SSoG"
                    ),
         month_n > 3 & month_n < 11) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg))


# creel coverage patchy so constrain input data to reasonable months
catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         reg %in% c("JdFS", "SSoG"
                    ),
         month_n > 3 & month_n < 11
         ) %>% 
  mutate(year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort))

# prediction datasets 
pred_dat_comp1 <- group_split(stock_comp, reg) %>%
  map_dfr(., function(x) {
    expand.grid(
      reg = unique(x$reg),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  mutate(
    reg_month_year = paste(reg, month_n, year, sep = "_"),
    key_var = as.factor(reg_month_year)
  )

pred_dat_catch <- group_split(catch, reg) %>%
  map_dfr(., function(x) {
    expand.grid(
      area = unique(x$area),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  # add region IDs back in
  left_join(.,
            catch %>% select(reg, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, reg, sep = "_"),
         offset = mean(catch$offset)) %>%
  arrange(reg) %>%
  # convoluted to ensure ordering is correct for key passed to TMB
  mutate(
    month = as.factor(month_n),
    order = row_number(),
    # year is necessary regardless of whether RIs are included because of year-
    # specific smooths 
    reg_month_year = paste(reg, month_n, year, sep = "_"),
    reg_month_year_f = fct_reorder(factor(reg_month_year), order)
  ) %>%
  select(-order, -strata) %>%
  distinct() %>% 
  # remove years that lack gsi data
  filter(
    year %in% pred_dat_comp1$year,
    area %in% c("121", "21", "20", "19_JdFS", "19_SSoG", "18")
  ) %>% 
  # generate predictions only for swiftsure
  rename(key_var = reg_month_year_f) %>% 
  droplevels()


# subset predicted composition dataset to match pred_dat_catch since fitting 
# data can be more extensive
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  filter(key_var %in% pred_dat_catch$key_var) %>% 
  # filter(reg_month_year %in% pred_dat_catch$reg_month_year) %>% 
  arrange(reg, year, month_n) %>% 
  droplevels()

pred_dat_catch_obs <- catch %>% 
  arrange(reg) %>%
  mutate(
    month = as.factor(month_n),
    order = row_number(),
    reg_month_year = paste(reg, month_n, year, sep = "_"),
    key_var = fct_reorder(factor(reg_month_year), order)
  )

model_inputs <- make_inputs(
  abund_formula = catch ~ 0 + area + 
    s(month_n, bs = "tp", k = 3) +
    s(month_n, by = area, bs = "tp", m = 1, k = 3) +
    s(month_n, by = year, bs = "tp", m = 1, k = 3) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  pred_abund = pred_dat_catch,
  comp_formula = coarse_agg ~ reg + 
    s(month_n, bs = "tp", k = 3, by = reg),
    # s(month_n, bs = "cc", k = 4, by = reg),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_comp = pred_dat_stock_comp,
  model = "integrated"
)

## check 
# names(model_inputs$tmb_data)
# model_inputs$tmb_data$pred_X1_ij %>% dim()
# model_inputs$tmb_data$pred_Zs[[1]] %>% dim
# model_inputs$tmb_data$pred_X2_ij %>% dim()
# model_inputs$tmb_data$pred_rfac_agg %>% unique() %>% length
# model_inputs$tmb_data$pred_rfac_agg_levels %>% length


mod <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  model = "integrated",
  fit_random = FALSE
)

ssdr <- mod$ssdr


## total abundance -------------------------------------------------------------

# Area abundance 
log_abund <- ssdr[rownames(ssdr) == "pred_mu1", ]

log_abund_preds_area <- data.frame(
  link_abund_est = log_abund[ , "Estimate"],
  link_abund_se =  log_abund[ , "Std. Error"]
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 


pred_dat_catch %>%
  cbind(., log_abund_preds_area) %>% 
  ggplot(data = ., aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_abund_est, colour = year)) +
  facet_wrap(~area) #+ 
  # geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
  #           alpha = 0.5)


catch %>% 
  mutate(
    cpue = catch / effort
  ) %>% 
  group_by(year, reg, area, month_n) %>% 
  summarize(
    mean_cpue = mean(cpue)
  ) %>% 
  ggplot(data = ., aes(x = month_n, y = mean_cpue, color = year)) +
  geom_point() +
  facet_wrap(~area)


## composition -----------------------------------------------------------------

logit_pred_ppn <- ssdr[rownames(ssdr) == "logit_pred_Pi_prop", ]

link_preds <- data.frame(
  link_prob_est = logit_pred_ppn[ , "Estimate"],
  link_prob_se =  logit_pred_ppn[ , "Std. Error"]
) %>% 
  mutate(
    pred_prob_est = plogis(link_prob_est),
    pred_prob_low = plogis(link_prob_est + (qnorm(0.025) * link_prob_se)),
    pred_prob_up = plogis(link_prob_est + (qnorm(0.975) * link_prob_se))
  ) 

stock_seq <- colnames(model_inputs$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  split(., .$reg)


# compare to observed composition
obs_comp <- stock_comp %>% 
  group_by(reg, month) %>% 
  mutate(total_samples = sum(prob)) %>% 
  group_by(reg, month_n, total_samples, coarse_agg) %>% 
  summarize(
    agg_prob = sum(prob),
    agg_ppn = agg_prob / total_samples,
    .groups = "drop"
  ) %>%  
  distinct() %>% 
  rename(stock = coarse_agg)

purrr::map2(pred_comp, names(pred_comp), function (x, y) {
  ggplot(data = x, aes(x = month_n)) +
    labs(y = "Predicted Stock Proportion", x = "Month") +
    facet_wrap(~stock) +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_prob_est, colour = year)) +
    labs(title = y) +
    geom_ribbon(data = x,
                aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
                alpha = 0.2) +
    geom_point(data = obs_comp %>% filter(reg == y),
               aes(y = agg_ppn, x = month_n))
})


## abundance -------------------------------------------------------------------

log_abund <- ssdr[rownames(ssdr) == "log_pred_mu1_Pi", ]

log_abund_preds <- data.frame(
  link_abund_est = log_abund[ , "Estimate"],
  link_abund_se =  log_abund[ , "Std. Error"]
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 

stock_seq2 <- colnames(model_inputs$tmb_data$Y2_ik)

pred_abund <- purrr::map(stock_seq2, function (x) {
  dum <- pred_dat_catch
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>% 
  cbind(., log_abund_preds) 


ggplot(data = pred_abund, aes(x = month_n)) +
    labs(y = "Predicted Abundance Index", x = "Month") +
    facet_grid(area~stock, scales = "free_y") +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_abund_est, colour = year))# +
    # geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
    #             alpha = 0.5)
