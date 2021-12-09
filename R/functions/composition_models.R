### Compare size and stock composition models 
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
compile(here::here("src", "negbin_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines")))
compile(here::here("src", "dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_mvn")))
compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))


# DATA CLEAN -------------------------------------------------------------------

comp <- readRDS(here::here("data", "rec", "coarse_rec_comp.rds")) %>% 
  filter(!region == "NSoG") %>% 
  droplevels() %>% 
  rename(prob = agg_prob)
catch <- readRDS(here::here("data", "rec", "month_area_recCatch_clean.rds")) %>% 
  filter(!region == "NSoG") %>% 
  droplevels()

# prediction datasets 
pred_dat_comp1 <- group_split(comp, region) %>%
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
    reg_month_year = paste(region, month_n, year, sep = "_")
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
  rename(key_var = reg_month_year_f)


# subset predicted composition dataset to match pred_dat_catch since fitting 
# data can be more extensive
pred_dat_comp <- pred_dat_comp1 %>% 
  filter(reg_month_year %in% pred_dat_catch$reg_month_year) 


## FIT -------------------------------------------------------------------------

comp_inputs <- make_inputs(
  abund_formula = catch ~ area + s(month_n, bs = "tp", k = 4, by = region) +
    s(month_n, by = year, bs = "tp", m = 1, k = 4) +
    offset,
  abund_dat = catch,
  abund_rint = "year",
  pred_abund = pred_dat_catch,
  comp_formula = agg ~ region + s(month_n, bs = "cc", k = 4, 
                            by = region),
  comp_dat = comp,
  comp_rint = "year",
  pred_comp = pred_dat_comp,
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
        here::here("data", "model_fits", "combined_mvn_121_21_onlyB.rds"))







## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr_old <- readRDS(here::here("data", "model_fits", "combined_mvn_121_21_only.rds"))


sdr <- stock_mod$sdr
ssdr <- stock_mod$ssdr

unique(rownames(ssdr))

## predicted stock composition
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

stock_seq <- colnames(obs_comp)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) 

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = year))# +
# geom_ribbon(data = pred_comp,
#             aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
#             alpha = 0.5)
plot(p)


# generate observed proportions
# number of samples in an event
stock_seq <- colnames(obs_comp)
long_dat <- comp_dat %>% 
  mutate(samp_nn = apply(obs_comp, 1, sum), each = length(stock_seq)) %>% 
  pivot_longer(cols = c(PSD:`FR-early`), names_to = "stock", 
               values_to = "obs_count") %>% 
  mutate(obs_ppn = obs_count / samp_nn)
mean_long_dat <- long_dat %>% 
  group_by(month_n, year, region, stock) %>% 
  summarize(obs_ppn = mean(obs_ppn), .groups = "drop")

p + 
  geom_point(data = mean_long_dat, aes(x = month_n, y = obs_ppn, colour = year))


## predicted abundance
log_pred_abund <- ssdr[rownames(ssdr) == "log_pred_abund", ]

link_abund_preds <- data.frame(
  link_fit = log_pred_abund[ , "Estimate"],
  link_se_fit =  log_pred_abund[ , "Std. Error"]
) %>% 
  mutate(
    link_lo = link_fit + (qnorm(0.025) * link_se_fit),
    link_up = link_fit + (qnorm(0.975) * link_se_fit),
    fit = exp(link_fit),
    fit_lo = exp(link_lo),
    fit_up = exp(link_up)
  )

pred_abund <- purrr::map(colnames(obs_comp), function (x) {
  dum <- pred_dat_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_abund_preds) %>% 
  #scale abundance for heat plot 
  group_by(stock, region) %>% 
  mutate(
    fit_z = (fit - mean(fit)) / sd(fit)
  ) %>% 
  ungroup() 
saveRDS(pred_abund, 
        here::here("data", "model_fits", 
                   "combined_mvn_121_21_only_pred_abund.rds"))


p <- ggplot(data = pred_abund, aes(x = month_n)) +
  labs(y = "Predicted Abundance Index", x = "Month") +
  facet_grid(region~stock, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = fit, colour = year)) +
  geom_ribbon(aes(ymin = fit_lo, ymax = fit_up, fill = year),
              alpha = 0.5)
plot(p)



ggplot(pred_abund, 
       aes(x = month_n, y = year)) +
  geom_raster(aes(fill = log(fit))) +
  scale_fill_viridis_c(name = "Predicted\nAbundance\nAnomalies") +
  labs(x = "Month", y = "Year") +
  facet_grid(region~stock) +
  ggsidekick::theme_sleek()


ggplot(data = pred_abund %>% filter(region == "JdFS"),  
       aes(x = month_n, y = year, fill = year)) +
  geom_ridgeline(aes(height = fit), alpha = 0.4, scale = 0.0004) +
  scale_fill_viridis_d() +
  labs(x = "Month", y = "Year") +
  facet_wrap(~ stock) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(#breaks = seq(2, 12, by = 2), limits = c(1, 12),
    #labels = month_labs, 
    expand = c(0, 0)) 
