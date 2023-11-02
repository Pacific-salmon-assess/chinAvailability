### Random Intercepts Integrated
## Adapt stockseasonr model with MVN intercepts
## For now test with area level predictions
## Aug 9, 2021

library(tidyverse)
library(TMB)


# tmb models - use MVN if time-varying predictions are required, use RI if 
# generating predictions for "average" year
compile(here::here("src", "dirichlet_ri_sdmTMB.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_ri_sdmTMB")))
compile(here::here("src", "negbin_rsplines_sdmTMB.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_sdmTMB")))
compile(here::here("src", "integrated_ri_sdmTMB.cpp"))
dyn.load(dynlib(here::here("src", "integrated_ri_sdmTMB")))


# data prep and model fitting functions
source(here::here("R", "functions", "fit_new.R"))


## composition inputs
# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal") %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date),
         reg = case_when(
           # subarea == "123I" ~ "SWVI", # qu
           #correction for subarea 21A IDd as SSoG
           subarea == "21A" ~ "SWVI",
           area %in% c("121", "21") ~ "SWVI",
           # subarea == "19C" ~ "SSoG",
           area %in% c("20W", "20E", "19JDF") ~ "JdFS",
           area %in% c("18", "19GST") ~ "SSoG",
           # area  ~ "SSoG",
           TRUE ~ "out"
         ),
         area = case_when(
           area %in% c("20E", "20W") ~ "20",
           area == "19JDF" ~ "19_JdFS",
           area == "19GST" ~ "19_SSoG",
           TRUE ~ area
         ),
         reg = as.factor(reg),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         sample_id = paste(month_n, area, week, year, sep = "_"),
         can_reg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                          "CR-lower_sp", "CR-upper_sp", "PSD", 
                          "NBC_SEAK", "WACST") ~ "Other",
           Region1Name == "SOMN" ~ "ECVI",
           TRUE ~ Region1Name
         ),
         pst_agg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                          "CR-lower_sp", "CR-upper_sp", 
                          "NBC_SEAK", "WACST") ~ "other",
           TRUE ~ pst_agg
         )
  ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id)) %>% as.numeric) %>% 
  ungroup()

stock_comp <- comp1 %>%  
  group_by(sample_id, 
           #subarea, subarea_original, 
           area, reg, reg_c = region, 
           week, month, month_n, year, nn, pst_agg
           # , core_area
           ) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(!reg == "out",
         # week > 22 & week < 38
         month_n < 10.1 & month_n > 1.9
         ) %>% 
  mutate(year = as.factor(year),
         reg = factor(reg, levels = c("SWVI", "JdFS", "SSoG")),
         #subarea = as.factor(subarea),
         area = as.factor(area)
         )


## catch inputs
catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         area %in% unique(stock_comp$area),
         month_n > 4 & month_n < 10) %>% 
  mutate(year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort))


## data coverage
stock_comp %>%
  select(sample_id, year, reg, area, month_n, nn) %>%
  distinct() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(area, as.numeric(reg))) +
  ggsidekick::theme_sleek()

ggplot(catch) +
  geom_point(aes(x = month_n, y = year, colour = reg),
              alpha = 0.5) +
  facet_wrap(~fct_reorder(area, as.numeric(reg))) +
  ggsidekick::theme_sleek()


# prediction datasets 
pred_dat1 <- expand.grid(
      reg = unique(catch$reg),
      month_n = seq(min(catch$month_n),
                    max(catch$month_n),
                    by = 0.1
      ))

# add areas to composition dataset
area_key <- stock_comp %>% 
  select(area, reg) %>% 
  distinct()

# subset predicted composition dataset
pred_dat <- pred_dat1 %>% 
  left_join(., area_key, by = "reg") %>%
  filter(month_n < 9 & month_n > 5
         )


## FIT MODEL -------------------------------------------------------------------

source(here::here("R", "functions", "fit_new.R"))

# no rand predictions
fit1 <- fit_stockseasonr(
  abund_formula = catch ~ 1 + area +
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  abund_dat = catch, 
  abund_offset = catch$offset,
  comp_formula = pst_agg ~ 1 + area +
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = stock_comp,
  pred_dat = pred_dat,
  model = "integrated",
  fit = TRUE,
  nlminb_loops = 2
)

ssdr <- fit1$ssdr 

beta_mat <- ssdr[rownames(ssdr) == "B2_jk", 2] %>% 
  matrix(., 
            nrow = ncol(fit1$tmb_data$X2_ij),
            ncol = ncol(fit1$tmb_data$Y2_ik))
rownames(beta_mat) <- colnames(fit1$tmb_data$X2_ij)
colnames(beta_mat) <- colnames(fit1$tmb_data$Y2_ik)
beta_mat



## EVALUATE MODEL PREDS --------------------------------------------------------


## stock composition
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

stock_seq <- colnames(fit1$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) 


ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = reg)) +
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)


## stock-specific abundance
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

stock_seq2 <- c(colnames(fit1$tmb_data$Y2_ik))

pred_abund <- purrr::map(stock_seq2, function (x) {
  dum <- pred_dat
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., log_abund_preds) %>% 
  mutate(stock = fct_reorder(factor(stock), -pred_abund_est))

p <- ggplot(data = pred_abund, aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  facet_grid(area~stock, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_abund_est))

p_ribbon <- p +
  geom_ribbon(data = x,
              aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
              alpha = 0.2)