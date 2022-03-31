### Integrated stock composition models 
## Updated with new data pull Feb 10; originals saved in /defunct
## NOTE: this model structure is more complex because predictions are scaled
## differently for abundance and composition component; this functionality has 
## been removed from default fit utility functions (saved as fit_EXPANDED.R)
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
compile(here::here("src", "negbin_rsplines_dirichlet_mvn_EXPANDED.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn_EXPANDED")))


# DATA CLEAN -------------------------------------------------------------------

# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal") %>% 
  mutate(
    reg = abbreviate(cap_region, minlength = 4),
    yday = lubridate::yday(date),
    month_n = lubridate::month(date),
    sample_id = paste(month_n, reg, yday, year, sep = "_"),
    pst_agg = case_when(
      pst_agg %in% c("CA_ORCST", "CR-lower_sp", "CR-upper_sp", "CR-upper_su/fa",
                     "NBC_SEAK", "WACST") ~ "other",
      TRUE ~ pst_agg
    )
  ) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id))
  ) %>% 
  ungroup()
stock_comp <- comp1 %>%  
  group_by(sample_id, area, reg, reg_c = cap_region, month, month_n, year, nn, 
           pst_agg) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(reg %in% c("JdFS", "SSoG")) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg),
         area = as.factor(area))


# creel coverage patchy so constrain input data to reasonable months
catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         reg %in% c("JdFS", "SSoG"),
         month_n > 4 & month_n < 10) %>% 
  mutate(year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort))


# prediction datasets 
area_key <- stock_comp %>% 
  select(area, reg) %>% 
  distinct()

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

# ggplot(pred_dat_catch) +
#   geom_point(aes(x = month_n, y = area)) +
#   facet_wrap(~reg)

# subset predicted composition dataset to match pred_dat_catch since fitting 
# data can be more extensive
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  filter(key_var %in% pred_dat_catch$key_var) %>% 
  # filter(reg_month_year %in% pred_dat_catch$reg_month_year) %>% 
  arrange(reg, year, month_n) %>% 
  droplevels()

stock_comp %>% 
  group_by(reg) %>% 
  summarize(
    min_m = min(month_n),
    max_m = max(month_n)
  )



## RAW OBSERVATIONS ------------------------------------------------------------

catch %>% 
  mutate(cpue = catch / effort) %>% 
  ggplot(.) +
  geom_point(aes(x = month, y = cpue, fill = region), shape = 21) +
  facet_wrap(~area, scales = "free_y")


## FIT -------------------------------------------------------------------------

model_inputs <- make_inputs(
  abund_formula = catch ~ 0 + area + 
    s(month_n, bs = "tp", k = 3, by = area) +
    s(month_n, by = year, bs = "tp", m = 1, k = 3) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  pred_abund = pred_dat_catch,
  comp_formula = pst_agg ~ area + 
    s(month_n, bs = "cc", k = 4, by = reg, m = 2),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_comp = pred_dat_stock_comp,
  model = "integrated"
)

stock_mod <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  model = "integrated",
  fit_random = FALSE
)

saveRDS(stock_mod$ssdr, 
        here::here("data", "model_fits", "combined_mvn_mig_corridor.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- stock_mod$ssdr
ssdr <- readRDS(
  here::here("data", "model_fits", "combined_mvn_mig_corridor.rds"))

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

stock_seq <- colnames(model_inputs$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  split(., .$reg)

comp_plots <- purrr::map2(pred_comp, names(pred_comp), function (x, y) {
  p <- ggplot(data = x, aes(x = month_n)) +
    labs(y = "Predicted Stock Proportion", x = "Month") +
    facet_wrap(~stock) +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_prob_est, colour = year)) +
    labs(title = y)
  
  p_ribbon <- p +
    geom_ribbon(data = x,
                aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
                alpha = 0.2)
  
  list(p, p_ribbon)
})


##### TO BE CORRECTED #####

# generate observed proportions
# number of samples in an event
long_dat <- model_inputs$wide_comp_dat %>%
  mutate(
    samp_nn = apply(model_inputs$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(PSD:CA_ORCST), names_to = "stock", 
               values_to = "obs_count") %>% 
  mutate(obs_ppn = obs_count / samp_nn) 
mean_long_dat <- long_dat %>%
  group_by(month_n, year, region, stock) %>%
  summarize(obs_ppn = mean(obs_ppn), .groups = "drop")

p_obs <- p + 
  geom_jitter(data = long_dat, aes(x = month_n, y = obs_ppn, colour = year,
                                  alpha = sample_size_bin)) +
  scale_alpha_discrete()


## Predicted Abundance

## aggregate (i.e. total)
# log_agg_abund <- ssdr[rownames(ssdr) == "ln_pred_mu1_cumsum", ]
log_abund <- ssdr[rownames(ssdr) == "log_pred_mu1_Pi", ]

log_abund_preds <- data.frame(
  link_abund_est = c(#log_agg_abund[ , "Estimate"], 
    log_abund[ , "Estimate"]),
  link_abund_se =  c(#log_agg_abund[ , "Std. Error"], 
    log_abund[ , "Std. Error"])
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 

stock_seq2 <- c(#"total", 
  colnames(model_inputs$tmb_data$Y2_ik))

pred_abund <- purrr::map(stock_seq2, function (x) {
  dum <- pred_dat_catch
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., log_abund_preds) %>% 
  mutate(stock = fct_reorder(factor(stock), -pred_abund_est)) %>% 
  split(., .$reg)

saveRDS(pred_abund, 
        here::here("data", "model_fits", 
                   "combined_mvn_mig_corridor_pred_abund.rds"))


abund_plots <- purrr::map2(pred_abund, names(pred_abund), function (x, y) {
  q <- ggplot(data = x, aes(x = month_n)) +
    labs(y = "Predicted Abundance Index", x = "Month") +
    facet_grid(area~stock, scales = "free_y") +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_abund_est, colour = year)) +
    labs(title = y)
  q_ribbon <- q + 
    geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
                alpha = 0.5)

  list(q, q_ribbon)  
})


## EXPORT FIGS -----------------------------------------------------------------

pdf(here::here("figs", "jdf_area_preds", "stock_comp_2regions_fe.pdf"))
comp_plots
dev.off()

pdf(here::here("figs", "jdf_area_preds", "abundance_2regions_fe.pdf"))
abund_plots
dev.off()