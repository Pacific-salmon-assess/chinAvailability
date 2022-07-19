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
source(here::here("R", "functions", "fit_new.R"))


# relevant TMB models
# compile(here::here("src", "negbin_rsplines.cpp"))
# dyn.load(dynlib(here::here("src", "negbin_rsplines")))
compile(here::here("src", "negbin_rsplines_sdmTMB.cpp"))
dyn.load(dynlib('src/negbin_rsplines_sdmTMB'))

# DATA CLEAN -------------------------------------------------------------------

catch <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         # reg %in% c("JdFS"#, "SSoG"
         # ),
         month_n > 4 & month_n < 10) %>% 
  mutate(year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort))
# model predictions are sensitive to regions included because information is 
# shared among smooths; remove northern regions


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
  left_join(.,
            catch %>% select(reg, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, reg, sep = "_"),
         offset = mean(catch$offset)) %>%
  # filter(strata %in% comp_strata) %>%
  arrange(reg) %>%
  # convoluted to ensure ordering is correct for key
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
  # generate predictions only for swiftsure
  filter(area %in% c("121", "21")) %>% 
  rename(key_var = reg_month_year_f)


## COMPARE FIT -----------------------------------------------------------------

tmb_inputs <- make_inputs(
  abund_formula = catch ~ 1 +
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | reg) +
    (1 | year)
  ,
  abund_dat = catch,
  abund_offset = "offset",
  pred_dat = catch,
  model = "negbin",
  include_re_preds = FALSE
)

abund_mod <- fit_model(
  tmb_data = tmb_inputs$tmb_data, 
  tmb_pars = tmb_inputs$tmb_pars, 
  tmb_map = tmb_inputs$tmb_map,
  tmb_random  = tmb_inputs$tmb_random,
  model_specs = tmb_inputs$model_specs
)
abund_mod$sdr


m3 <- sdmTMB(catch ~ 1 +
               s(month_n, bs = "tp", k = 4, m = 2) +
               (1 | reg) +
               (1 | year),
             offset = catch$offset,
             data = catch, 
             spatial = "off", family = sdmTMB::nbinom2())
m3$sd_report


ssdr <- abund_mod$ssdr
pred_mu <- ssdr[rownames(ssdr) == "pred_mu1", "Estimate"]

pred_mu_sdm <- predict(m3, newdata = catch) %>% glimpse()



### IGNORE BELOW ####


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- readRDS(
  here::here("data", "model_fits", "negbin_rsmooths_121_21_only.rds"))

unique(rownames(ssdr))


## Area abundance 
log_abund <- ssdr[rownames(ssdr) == "ln_pred_mu1", ]

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
  facet_wrap(~area)
q + 
  geom_ribbon(aes(ymin = pred_abund_low, ymax = pred_abund_up, fill = year),
              alpha = 0.5)


## Regional abundance (i.e. total)
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


## EXPLORE DHARMa RESIDUALS ----------------------------------------------------

## presently does not support dirichlet or multinomial so testing negative 
# binomial independently

library(DHARMa)
library(mgcViz)

nb_gam <- gam(catch ~ 0 + area + s(month_n, bs = "tp", k = 4, by = region) +
                s(month_n, by = year, bs = "tp", m = 1, k = 4) +
                offset + s(year, bs = "re"),
              data = catch,
              family = mgcv::nb())

gam.check(nb_gam)

# throws errors (simulate directly)
# simulationOutput <- simulateResiduals(fittedModel = nb_gam, plot = T)
# plot(simulationOutput)

#
simulate_nb_gam <- function(n){
  pred = predict(nb_gam)
  nObs = length(pred)
  sim = matrix(nrow = nObs, ncol = n)
  for(i in 1:n) sim[,i] = rnbinom(nObs, size = 2.265, mu = 1 / (exp(pred)))
  return(sim)
}
dum <- simulate_nb_gam(n = 100)

rnbinom(nObs, size = 2.265, mu = 1 / (exp(pred)))


dum <- simulate(nb_gam, nsim = 50)

DHARMaRes <- createDHARMa(
  simulatedResponse = dum, observedResponse = catch$catch, 
  fittedPredictedResponse = as.numeric(predict(nb_gam, type = "response")), 
  integerResponse = F
)
plot(DHARMaRes, quantreg = F)



testData = createData(sampleSize = 200, overdispersion = 0.5, family = poisson())
fittedModel <- glm(observedResponse ~ Environment1, family = "poisson", data = testData)

simulatePoissonGLM <- function(fittedModel, n){
  pred = predict(fittedModel, type = "response")
  nObs = length(pred)
  sim = matrix(nrow = nObs, ncol = n)
  for(i in 1:n) sim[,i] = rpois(nObs, pred)
  return(sim)
}

sim = simulatePoissonGLM(fittedModel, 100)
DHARMaRes = createDHARMa(
  simulatedResponse = sim, observedResponse = testData$observedResponse, 
  fittedPredictedResponse = predict(fittedModel), integerResponse = T)
plot(DHARMaRes, quantreg = F)



## SANDBOX ---------------------------------------------------------------------

## pcod example for fitting equivalent
# pcod$fyear = as.factor(pcod$year)
# pcod$catch = round(pcod$density, digits = 0)
# tmb_inputs <- make_inputs(
#   abund_formula = catch ~ 1 +
#     s(depth) +
#     (1 | fyear),
#   abund_dat = pcod,
#   pred_dat = pcod,
#   model = "negbin",
#   include_re_preds = FALSE
# )
# # 
# m2 <- sdmTMB(catch ~ 1 + s(depth), #+ (1|fyear),
#              data = pcod,
#   spatial = "off", family = sdmTMB::nbinom2())
# m2$sd_report
# 
# obj <- TMB::MakeADFun(
#   data = tmb_inputs$tmb_data,
#   parameters = tmb_inputs$tmb_pars,
#   random = tmb_inputs$tmb_random,
#   map = NULL,
#   DLL = "negbin_rsplines_sdmTMB"
# )
# opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
#   control = list(eval.max = 1e4, iter.max = 1e4))
# opt
# sdr <- sdreport(obj)
# sdr
# 
# m2$sd_report
