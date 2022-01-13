### Integrated size composition models 
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

size_dat <- readRDS(here::here("data", "rec", "rec_size_ppn.rds")) %>% 
  filter(region %in% c("SSoG", "JdFS")) %>% 
  droplevels() %>% 
  rename(prob = sum_count)

catch2 <- readRDS(
  here::here("data", "rec", "month_area_recCatch_clean.rds")
  ) %>% 
  filter(region %in% c("SSoG", "JdFS")) %>% 
  droplevels()
# model predictions are sensitive to regions included because information is 
# shared among smooths; remove northern regions


# visualize sampling coverage



# PREDICTIVE DATASETS ----------------------------------------------------------

# prediction datasets 
pred_dat_comp1 <- group_split(size_dat, region) %>%
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


pred_dat_catch2 <- group_split(catch2, region) %>%
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
  left_join(.,
            catch2 %>% select(region, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, region, sep = "_"),
         offset = mean(catch2$offset)) %>%
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
pred_dat_size_comp <- pred_dat_comp1 %>% 
  filter(key_var %in% pred_dat_catch2$key_var) %>% 
  # filter(reg_month_year %in% pred_dat_catch$reg_month_year) %>% 
  arrange(region, year, month_n) %>% 
  droplevels()


## FIT -------------------------------------------------------------------------

size_comp_inputs <- make_inputs(
  abund_formula = catch ~ 0 + area +
    s(month_n, bs = "tp", k = 3, by = region) +
    s(month_n, by = year, bs = "tp", m = 1, k = 3) +
    offset,
  abund_dat = catch2,
  abund_rint = "year",
  pred_abund = pred_dat_catch2,
  comp_formula = size_bin ~ region + s(month_n, bs = "cc", k = 4, by = region),
  comp_dat = size_dat,
  comp_rint = "year",
  pred_comp = pred_dat_size_comp,
  model = "integrated"
  # model = "dirichlet"
)

size_mod <- fit_model(
  tmb_data = size_comp_inputs$tmb_data, 
  tmb_pars = size_comp_inputs$tmb_pars, 
  tmb_map = size_comp_inputs$tmb_map, 
  tmb_random  = size_comp_inputs$tmb_random,
  # model = "dirichlet",
  model = "integrated",
  fit_random = TRUE
)

saveRDS(size_mod$ssdr, 
        here::here("data", "model_fits", "size",
                   "size_combined_mvn_121_21_only.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- readRDS(
  here::here("data", "model_fits", "size", "size_combined_mvn_121_21_only.rds"))

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

size_seq <- colnames(size_comp_inputs$tmb_data$Y2_ik)
pred_comp <- purrr::map(size_seq, function (x) {
  dum <- pred_dat_size_comp
  dum$size_bin <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  mutate(
    size_bin = fct_relevel(size_bin, "<45", "45-60", "60-75", "75-90", ">90"))

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Size Proportion", x = "Month") +
  facet_wrap(~size_bin) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = year)) 

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
              alpha = 0.2)

# generate observed proportions
# number of samples in an event
long_dat <- size_comp_inputs$wide_comp_dat %>%
  mutate(
    samp_nn = apply(size_comp_inputs$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(`<45`:`>90`), names_to = "size_bin", 
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

size_seq2 <- c("total", colnames(size_comp_inputs$tmb_data$Y2_ik))

pred_abund <- purrr::map(size_seq2, function (x) {
  dum <- pred_dat_size_comp
  dum$size_bin <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., log_abund_preds) %>% 
  mutate(
    size_bin = fct_relevel(size_bin, 
                           "total", "<45", "45-60", "60-75", "75-90", ">90"))

saveRDS(pred_abund, 
        here::here("data", "model_fits", "size",
                   "combined_mvn_121_21_only_pred_abund.rds"))


q <- ggplot(data = pred_abund, aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  facet_wrap(~size_bin, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_abund_est, colour = year))  +
  scale_colour_viridis_d()
q_ribbon <- q + 
  geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
              alpha = 0.5) +
  scale_fill_viridis_d() 


catch2 %>% 
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

pdf(here::here("figs", "jdf_area_preds", "size_comp_3knots.pdf"))
p
p_ribbon
p_obs
dev.off()

pdf(here::here("figs", "jdf_area_preds", "abundance_2regions.pdf"))
q
q_ribbon
dev.off()