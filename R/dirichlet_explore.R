### MVN Random Intercepts Dirichlet
## Adapt stockseasonr model with MVN intercepts
## NOTE: attempted to incorporate random smooths, but unclear how to proceed 
## given parameters are a matrix, not vector
## Nov. 29, 2021
## Updated March 19 to explore composition at area levels

library(tidyverse)
library(TMB)


# tmb models
# compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
# dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))
compile(here::here("src", "dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_mvn")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


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
  mutate(nn = length(unique(id))) %>% 
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
         area = as.factor(area),
         )


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


# add areas to composition dataset
area_key <- stock_comp %>% 
  select(area, reg) %>% 
  distinct()

# subset predicted composition dataset
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  arrange(reg, year, month_n) %>% 
  droplevels() #%>% 
  # left_join(., area_key, by = "reg") %>%
  # filter(area %in% c("121", "21", "20", "19JdF", "19GST", "20W", "18"),
  #        month_n < 9.1 & month_n > 4.9) 

# pred_dat_stock_comp2 <- pred_dat_stock_comp %>% 
#   mutate(key_var = paste(reg, month_n, sep = "_")) %>% 
#   select(-c(year, reg_month_year)) %>% 
#   distinct() %>% 
#   glimpse()

## FIT MODEL -------------------------------------------------------------------

model_inputs <- make_inputs(
  comp_formula = pst_agg ~ reg + 
    s(month_n, bs = "cc", k = 4, by = reg, m = 2),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_comp = pred_dat_stock_comp,
  model = "dirichlet"
)

stock_mod <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  model = "dirichlet",
  fit_random = TRUE,
  ignore_fix = TRUE
)

saveRDS(stock_mod$ssdr, 
        here::here("data", "model_fits", 
                   "dirichlet_area_int_mvn_mig_corridor.rds"))

# fit second model without RE predictions for comparison
# model_inputs2 <- make_inputs(
#   comp_formula = pst_agg ~ reg + 
#     s(month_n, bs = "cc", k = 4, by = reg, m = 2),
#   comp_dat = stock_comp,
#   comp_rint = "year",
#   pred_comp = pred_dat_stock_comp2,
#   model = "dirichlet"
# )
# 
# stock_mod2 <- fit_model(
#   tmb_data = model_inputs2$tmb_data, 
#   tmb_pars = model_inputs2$tmb_pars, 
#   tmb_map = model_inputs2$tmb_map, 
#   tmb_random  = model_inputs2$tmb_random,
#   model = "dirichlet",
#   fit_random = TRUE,
#   ignore_fix = TRUE
# )


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- stock_mod$ssdr
# ssdr <- readRDS(
#   here::here("data", "model_fits", "dirichlet_mvn_mig_corridor.rds"))

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
    facet_grid(area~stock) +
    ggsidekick::theme_sleek() +
    geom_line(aes(y = pred_prob_est, colour = year)) +
    labs(title = y)
  
  p_ribbon <- p +
    geom_ribbon(data = x,
                aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
                alpha = 0.2)
  
  list(p, p_ribbon)
})

pdf(here::here("figs", "jdf_area_preds", "stock_comp_area_preds.pdf"))
comp_plots
dev.off()


## compare to observations
# number of samples in an event
long_dat <- model_inputs$wide_comp_dat %>%
  mutate(samp_nn = apply(model_inputs$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(PSD:WCVI), names_to = "stock", 
               values_to = "obs_count") %>% 
  # filter(area %in% c("121", "21", "20", "19JdF", "19GST", "20W", "18")) %>% 
  mutate(obs_ppn = obs_count / samp_nn) 
mean_long_dat <- long_dat %>%
  group_by(month_n, year, area, reg, stock) %>%
  summarize(obs_ppn = mean(obs_ppn), .groups = "drop") %>%
  split(., .$reg)

ggplot() +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area~stock) +
  ggsidekick::theme_sleek() +
  geom_line(data = pred_comp[[1]], aes(x = month_n, y = pred_prob_est, 
                                       colour = year)) +
  geom_jitter(data = mean_long_dat[[1]], 
              aes(x = month_n, y = obs_ppn, colour = year#,
                                   # alpha = sample_size_bin
                                   )) +
  scale_alpha_discrete()


