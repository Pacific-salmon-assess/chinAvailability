### Random Splines Neg Binomial
## Version of integrated_stock_composition but with only negative binomial 
## component of the model and fit to a larger dataset
## Nov. 12, 2021

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


# relevant TMB models
compile(here::here("src", "negbin_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines")))



# DATA CLEAN -------------------------------------------------------------------

catch <- readRDS(here::here("data", "rec", "month_area_recCatch_clean.rds")) %>% 
  filter(region %in% c("SSoG", "JdFS")) %>% 
  droplevels()
# model predictions are sensitive to regions included because information is 
# shared among smooths; remove northern regions


pred_dat_catch <- group_split(catch, region) %>%
  map_dfr(., function(x) {
    expand.grid(
      # region = unique(x$region),
      area = unique(x$area),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  left_join(.,
            catch %>% select(region, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, region, sep = "_"),
         offset = mean(catch$offset)) %>%
  # filter(strata %in% comp_strata) %>%
  arrange(region) %>%
  # convoluted to ensure ordering is correct for key
  mutate(
    month = as.factor(month_n),
    order = row_number(),
    # year is necessary regardless of whether RIs are included because of year-
    # specific smooths 
    reg_month_year = paste(region, month_n, year, sep = "_"),
    reg_month_year_f = fct_reorder(factor(reg_month_year), order)
  ) %>%
  select(-order, -strata) %>%
  distinct() %>% 
  # generate predictions only for swiftsure
  filter(area %in% c("121", "21")) %>% 
  rename(key_var = reg_month_year_f)


## FIT -------------------------------------------------------------------------

tmb_inputs <- make_inputs(
  abund_formula = catch ~ 0 + area + 
    s(month_n, bs = "tp", k = 4, by = region) +
    s(month_n, by = year, bs = "tp", m = 1, k = 4) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  pred_abund = pred_dat_catch,
  model = "negbin"
)

abund_mod <- fit_model(
  tmb_data = tmb_inputs$tmb_data, 
  tmb_pars = tmb_inputs$tmb_pars, 
  tmb_map = tmb_inputs$tmb_map, 
  tmb_random  = tmb_inputs$tmb_random,
  model = "negbin",
  fit_random = TRUE
)

saveRDS(abund_mod$ssdr, 
        here::here("data", "model_fits", "negbin_rsmooths_121_21_only.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- readRDS(
  here::here("data", "model_fits", "negbin_rsmooths_121_21_only.rds"))

unique(rownames(ssdr))


## Aggregate abundance (i.e. total)
log_agg_abund <- ssdr[rownames(ssdr) == "ln_pred_mu1_cumsum", ]

log_abund_preds <- data.frame(
  link_abund_est = log_agg_abund[ , "Estimate"],
  link_abund_se =  log_agg_abund[ , "Std. Error"]
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 

# drop area identifiers and collapse, then add predictions
pred_abund <- pred_dat_catch %>%
  select(-area) %>% 
  distinct() %>% 
  cbind(., log_abund_preds) 

saveRDS(pred_abund, 
        here::here("data", "model_fits", 
                   "negbin_rsmooths_121_21_only_pred_abund.rds"))


ggplot(data = pred_abund, aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_abund_est, colour = year)) 
q + 
  geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
              alpha = 0.5)