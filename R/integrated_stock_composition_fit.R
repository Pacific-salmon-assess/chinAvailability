### Integrated stock composition models 
## Updated with new data pull Feb 10; originals saved in /defunct
## Dec. 8, 2021

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


# relevant TMB models
# compile(here::here("src", "negbin_rsplines.cpp"))
# dyn.load(dynlib(here::here("src", "negbin_rsplines")))
# compile(here::here("src", "dirichlet_mvn.cpp"))
# dyn.load(dynlib(here::here("src", "dirichlet_mvn")))
compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))


# DATA CLEAN -------------------------------------------------------------------

# old data using in chin dist analysis
# stock_comp <- readRDS(here::here("data", "rec", "old", "coarse_rec_comp.rds")) %>% 
#   # filter(region %in% c("SSoG", "JdFS"),
#   #        month %in% c("6", "8"),
#   #        year %in% c("2016", "2018")) %>% 
#   droplevels() %>% 
#   rename(prob = agg_prob)
# 
# catch <- readRDS(here::here("data", "rec", "old", "month_area_recCatch_clean.rds")) %>% 
#   # filter(region %in% c("SSoG", "JdFS"),
#   #        month %in% c("6", "8"),
#   #        year %in% c("2016", "2018")) %>% 
#   droplevels()
# model predictions are sensitive to regions included because information is 
# shared among smooths; remove northern regions


# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal") %>% 
  mutate(reg = abbreviate(cap_region, minlength = 4),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         sample_id = paste(month_n, reg, yday, year, sep = "_")) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id))
  ) %>% 
  ungroup()
stock_comp <- comp1 %>%  
  group_by(sample_id, reg, reg_c = cap_region, month, month_n, year, nn, 
           pst_agg) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels()

# filter to SSoG and JdF


catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal") %>% 
  mutate(offset = log(effort))




# prediction datasets 
pred_dat_comp1 <- group_split(stock_comp, region) %>%
  map_dfr(., function(x) {
    expand.grid(
      region = unique(x$region),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  mutate(
    reg_month_year = paste(region, month_n, year, sep = "_"),
    key_var = as.factor(reg_month_year)
  )

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
  # remove years that lack gsi data
  filter(
    year %in% pred_dat_comp1$year
  ) %>% 
  # generate predictions only for swiftsure
  filter(area %in% c("121", "21")) %>% 
  rename(key_var = reg_month_year_f) %>% 
  droplevels()


# subset predicted composition dataset to match pred_dat_catch since fitting 
# data can be more extensive
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  filter(key_var %in% pred_dat_catch$key_var) %>% 
  # filter(reg_month_year %in% pred_dat_catch$reg_month_year) %>% 
  arrange(region, year, month_n) %>% 
  droplevels()


## FIT -------------------------------------------------------------------------

comp_inputs <- make_inputs(
  abund_formula = catch ~ 0 + area + 
    s(month_n, bs = "tp", k = 3, by = region) +
    s(month_n, by = year, bs = "tp", m = 1, k = 3) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  pred_abund = pred_dat_catch,
  comp_formula = agg ~ region + s(month_n, bs = "cc", k = 4, by = region),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_comp = pred_dat_stock_comp,
  model = "integrated"
)

stock_mod <- fit_model(
  tmb_data = comp_inputs$tmb_data, 
  tmb_pars = comp_inputs$tmb_pars, 
  tmb_map = comp_inputs$tmb_map, 
  tmb_random  = comp_inputs$tmb_random,
  model = "integrated",
  fit_random = TRUE
)

saveRDS(stock_mod$ssdr, 
        here::here("data", "model_fits", "combined_mvn_121_21_only.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- readRDS(
  here::here("data", "model_fits", "combined_mvn_121_21_only.rds"))

unique(rownames(ssdr))


## Stock Composition
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

stock_seq <- colnames(comp_inputs$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) 

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = year)) 

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
              alpha = 0.2)


# generate observed proportions
# number of samples in an event
long_dat <- comp_inputs$wide_comp_dat %>%
  mutate(
    samp_nn = apply(comp_inputs$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(PSD:`FR-early`), names_to = "stock", 
               values_to = "obs_count") %>% 
  mutate(obs_ppn = obs_count / samp_nn,
         sample_size_bin = cut(
           samp_nn, 
           breaks = c(-Inf, 5, 15, 25, Inf), 
           labels=c("<6", "6-15", "16-25", ">25")
         )
         ) %>% 
  filter(region == "JdFS") 
# mean_long_dat <- long_dat %>% 
#   group_by(month_n, year, region, stock) %>% 
#   summarize(obs_ppn = mean(obs_ppn), .groups = "drop")

p_obs <- p + 
  geom_jitter(data = long_dat, aes(x = month_n, y = obs_ppn, colour = year,
                                  alpha = sample_size_bin)) +
  scale_alpha_discrete()


## Predicted Abundance

## aggregate (i.e. total)
log_agg_abund <- ssdr[rownames(ssdr) == "ln_pred_mu1_cumsum", ]
log_abund <- ssdr[rownames(ssdr) == "log_pred_mu1_Pi", ]

log_abund_preds <- data.frame(
  link_abund_est = c(log_agg_abund[ , "Estimate"], log_abund[ , "Estimate"]),
  link_abund_se =  c(log_agg_abund[ , "Std. Error"], log_abund[ , "Std. Error"])
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 

stock_seq2 <- c("total", colnames(comp_inputs$tmb_data$Y2_ik))

pred_abund <- purrr::map(stock_seq2, function (x) {
  dum <- pred_dat_stock_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., log_abund_preds) %>% 
  mutate(stock = fct_reorder(factor(stock), -pred_abund_est))

saveRDS(pred_abund, 
        here::here("data", "model_fits", 
                   "combined_mvn_121_21_only_pred_abund.rds"))

q <- ggplot(data = pred_abund, aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  facet_wrap(~stock, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_abund_est, colour = year)) 
q_ribbon <- q + 
  geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
              alpha = 0.5)



catch %>% 
  group_by(region, month_n, year) %>% 
  summarize(sum_catch = sum(catch),
            sum_eff = sum(eff),
            cpue = sum_catch / sum_eff,
            .groups = "drop") %>% 
  ggplot(., aes(x = month_n)) +
  labs(y = "CPUE", x = "Month") +
  facet_grid(region~year, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_bar(aes(y = cpue), stat = "identity") 



## EXPORT FIGS -----------------------------------------------------------------

pdf(here::here("figs", "jdf_area_preds", "stock_comp_2regions.pdf"))
p
p_ribbon
p_obs
dev.off()

pdf(here::here("figs", "jdf_area_preds", "abundance_2regions.pdf"))
q
q_ribbon
dev.off()