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
source(here::here("R", "functions", "fit_new.R"))


# relevant TMB models
# compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
# dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))
compile(here::here("src", "negbin_rsplines_dirichlet_ri.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_ri")))


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
      pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa") ~ "CR-fa",
      pst_agg %in% c("CA_ORCST", "CR-lower_sp", "CR-upper_sp", 
                     "NBC_SEAK", "WACST") ~ "other",
      TRUE ~ pst_agg
    ),
    #rename specific areas to match catch data
    area = case_when(
      area %in% c("20E", "20W") ~ "20",
      area == "19JDF" ~ "19_JdFS",
      area == "19GST" ~ "19_SSoG",
      TRUE ~ area
    )
  ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id))) %>% 
  ungroup()
stock_comp <- comp1 %>%  
  group_by(sample_id, area, reg, reg_c = cap_region, month, month_n, year, nn, 
           pst_agg) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(reg %in% c("JdFS", "SSoG"),
         month_n < 10 & month_n > 1.9
         ) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg),
         area = as.factor(area),
  )

# creel coverage patchy so constrain input data to reasonable months
catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         reg %in% c("JdFS", "SSoG"),
         month_n > 4 & month_n < 10) %>% 
  mutate(year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort))


# prediction dataset
area_key <- stock_comp %>% 
  select(area, reg) %>% 
  distinct()

# generate predictive dataframe constrained to variables common to both
# datasets
pred_dat_integrated <- expand.grid(
  area = unique(stock_comp$area[stock_comp$area %in% catch$area]),
  month_n = seq(5, 9, by = 0.1)
) %>% 
  left_join(., area_key, by = "area")
pred_dat_integrated$offset <- mean(catch$offset)


## RAW OBSERVATIONS ------------------------------------------------------------

catch %>% 
  mutate(cpue = catch / effort) %>% 
  ggplot(.) +
  geom_point(aes(x = month, y = cpue, fill = region), shape = 21) +
  facet_wrap(~area, scales = "free_y")

stock_comp %>% 
  select(sample_id, year, reg, area, month_n, nn) %>% 
  distinct() %>% 
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg),
              alpha = 0.3, width = 0.25) +
  facet_wrap(~area)

unique(catch$area)
unique(stock_comp$area)


## FIT -------------------------------------------------------------------------

model_inputs <- make_inputs(
  abund_formula = catch ~ area +
    s(month_n, bs = "tp", k = 4, m = 2) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  comp_formula = pst_agg ~ area + s(month_n, bs = "tp", k = 4, m = 2),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_dat = pred_dat_integrated,
  model = "integrated",
  include_re_preds = FALSE,
  random_walk = TRUE
)

stock_mod <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  fit_random = TRUE,
  ignore_fix = FALSE,
  model_specs = model_inputs$model_specs
)

saveRDS(stock_mod$ssdr, 
        here::here("data", "model_fits", "integrated_area_ri_mig_corridor.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- stock_mod$ssdr
ssdr <- readRDS(
  here::here("data", "model_fits", "combined_mvn_mig_corridor.rds"))

unique(rownames(ssdr))


## Stock Composition (REMOVED FOR ADREPORT)
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
  dum <- pred_dat_integrated
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) 

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = reg))

p_ribbon <- p +
  geom_ribbon(data = x,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
              alpha = 0.2)



## Predicted Abundance

## aggregate (i.e. total)
log_abund <- ssdr[rownames(ssdr) == "log_pred_mu1_Pi", ]

log_abund_preds <- data.frame(
  link_abund_est = log_abund[ , "Estimate"],
  link_abund_se = log_abund[ , "Std. Error"]
) %>% 
  mutate(
    pred_abund_est = exp(link_abund_est),
    pred_abund_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
    pred_abund_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
  ) 

stock_seq2 <- c(colnames(model_inputs$tmb_data$Y2_ik))

pred_abund <- purrr::map(stock_seq2, function (x) {
  dum <- pred_dat_integrated
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
  p <- ggplot(data = x, aes(x = month_n)) +
    labs(y = "Predicted Abundance Index", x = "Month") +
    facet_wrap(area~stock, scales = "free_y") +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_abund_est, colour = year))
  
  p_ribbon <- p +
    geom_ribbon(data = x,
                aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
                alpha = 0.2)

  list(p, p_ribbon)
})







## EXPORT FIGS -----------------------------------------------------------------

pdf(here::here("figs", "jdf_area_preds", "stock_comp_2regions_fe.pdf"))
comp_plots
dev.off()

pdf(here::here("figs", "jdf_area_preds", "abundance_2regions_fe.pdf"))
abund_plots
dev.off()