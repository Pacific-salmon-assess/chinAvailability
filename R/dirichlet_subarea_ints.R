### MVN Random Intercepts Dirichlet
## Adapt stockseasonr model with MVN intercepts
## NOTE: attempted to incorporate random smooths, but unclear how to proceed 
## given parameters are a matrix, not vector
## Nov. 29, 2021
## Updated March 19 to explore composition at area levels
## Updated April 14 to subarea levels

library(tidyverse)
library(TMB)


# tmb models - use MVN if time-varying predictions are required, use RI if 
# generating predictions for "average" year
compile(here::here("src", "dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_mvn")))
compile(here::here("src", "dirichlet_ri.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_ri")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal",
         #remove subareas 18E, 19A, 18A w/ very small sample sizes 
         # (and outside study area)
         !subarea %in% c("18E", "19A", "18A")) %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date),
         subarea = case_when(
           subarea %in% c("21B") ~ "21A",
           subarea %in% c("US7") ~ "19B",
           TRUE ~ subarea
         ),
         reg = case_when(
           # subarea == "123I" ~ "SWVI", # qu
           #correction for subarea 21A IDd as SSoG
           subarea == "21A" ~ "SWVI",
           area %in% c("121", "21") ~ "SWVI",
           subarea == "19C" ~ "SSoG",
           area %in% c("20W", "20E", "19JDF") ~ "JdFS",
           area %in% c("18", "19GST") ~ "SSoG",
           # area  ~ "SSoG",
           TRUE ~ "out"
         ),
         reg = as.factor(reg),
         core_area = case_when(
           reg %in% c("SWVI") ~ "yes", 
           subarea %in% c("18B", "19B", "19C") ~ "yes",
           TRUE ~ "no"
         ),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         sample_id = paste(month_n, reg, yday, year, sep = "_"),
         pst_agg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa") ~ "CR-fa",
           pst_agg %in% c("CA_ORCST", "CR-lower_sp", "CR-upper_sp", 
                          "NBC_SEAK", "WACST") ~ "other",
           TRUE ~ pst_agg
         )
  ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id))) %>% 
  ungroup()

stock_comp <- comp1 %>%  
  group_by(sample_id, subarea, area, reg, reg_c = region, month, month_n, year, nn, 
           pst_agg, core_area) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(!reg == "out",
         month_n < 9.1 & month_n > 4.9) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg),
         subarea = as.factor(subarea),
         area = as.factor(area)
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
  }) 


# add areas to composition dataset
area_key <- stock_comp %>% 
  select(subarea, area, reg, core_area) %>% 
  distinct()

# subset predicted composition dataset
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  left_join(., area_key, by = "reg") %>%
  filter(core_area == "yes",
         #area %in% c("121", "21", "20", "19_JdFS", "19_SSoG", "18"),
         month_n < 9.1 & month_n > 5.9) %>% 
  select(-core_area)

# remove random effects from predictions to generate estimates for "average"
# year
pred_dat_stock_comp_ri <- pred_dat_stock_comp %>%
  select(-year) %>%
  distinct() %>%
  glimpse()


# look at sample coverage in data passed to model
alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")

png(here::here("figs", "data_coverage", "comp_model_inputs.png"))
stock_comp %>% 
  select(sample_id, year, reg, subarea, month_n, nn, core_area) %>% 
  distinct() %>% 
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg, 
                  shape = core_area),
              alpha = 0.3, width = 0.25) +
  facet_wrap(~subarea)
dev.off()

stock_comp %>% 
  select(sample_id, year, reg, area, month_n, nn, core_area) %>% 
  distinct() %>% 
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg, 
                  shape = core_area),
              alpha = 0.3, width = 0.25) +
  facet_wrap(~area)



## FIT MODEL -------------------------------------------------------------------


# rand predictions
model_inputs <- make_inputs(
  comp_formula = pst_agg ~ subarea + 
    s(month_n, bs = "tp", k = 3, by = reg, m = 2),
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_dat = pred_dat_stock_comp,
  model = "dirichlet",
  include_re_preds = FALSE
)

stock_mod <- fit_model(
  tmb_data = model_inputs$tmb_data, 
  tmb_pars = model_inputs$tmb_pars, 
  tmb_map = model_inputs$tmb_map, 
  tmb_random  = model_inputs$tmb_random,
  fit_random = TRUE,
  ignore_fix = TRUE,
  model_specs = model_inputs$model_specs
  )


saveRDS(stock_mod$ssdr, 
        here::here("data", "model_fits", 
                   "dirichlet_subarea_mvn_mig_corridor.rds"))


# no rand predictions
model_inputs_ri <- make_inputs(
  comp_formula = pst_agg ~ area + 
    s(month_n, bs = "tp", k = 4, by = reg, m = 2) ,
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_comp = pred_dat_stock_comp_ri,
  model = "dirichlet",
  include_re_preds = FALSE
)

stock_mod_ri <- fit_model(
  tmb_data = model_inputs_ri$tmb_data, 
  tmb_pars = model_inputs_ri$tmb_pars, 
  tmb_map = model_inputs_ri$tmb_map, 
  tmb_random  = model_inputs_ri$tmb_random,
  model = "dirichlet",
  fit_random = TRUE,
  ignore_fix = TRUE,
  include_re_preds = FALSE
)

saveRDS(stock_mod_ri$ssdr, 
        here::here("data", "model_fits", 
                   "dirichlet_area_ri_mig_corridor.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------


# random effects predictions 
ssdr <- readRDS(
  here::here("data", "model_fits", "dirichlet_area_mvn_mig_corridor.rds"))

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
  mutate(area_f = fct_reorder(as.factor(area), as.numeric(reg)))

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area_f~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = year)) #+

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = year),
              alpha = 0.2)




# fixed effects predictions 
ssdr_ri <- readRDS(
  here::here("data", "model_fits", "dirichlet_area_ri_mig_corridor.rds"))

logit_pred_ppn_ri <- ssdr_ri[rownames(ssdr_ri) == "logit_pred_Pi_prop", ]

link_preds_ri <- data.frame(
  link_prob_est = logit_pred_ppn_ri[ , "Estimate"],
  link_prob_se =  logit_pred_ppn_ri[ , "Std. Error"]
) %>% 
  mutate(
    pred_prob_est = plogis(link_prob_est),
    pred_prob_low = plogis(link_prob_est + (qnorm(0.025) * link_prob_se)),
    pred_prob_up = plogis(link_prob_est + (qnorm(0.975) * link_prob_se))
  ) 

pred_comp_ri <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp_ri
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds_ri) %>% 
  mutate(area_f = fct_reorder(as.factor(area), as.numeric(reg)))

p2 <- ggplot(data = pred_comp_ri, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area_f~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) +
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.5)


png(here::here("figs", "jdf_area_preds", "stock_comp_area_preds_ri.png"),
    height = 5, width = 7, units = "in", res = 200)
p2
dev.off()


## combined plot
pdf(here::here("figs", "jdf_area_preds", "stock_comp_area_preds.pdf"))
p + 
  geom_line(data = pred_comp_ri, aes(x = month_n, y = pred_prob_est),
            color = "black", size = 1.1)
p_ribbon + 
  geom_line(data = pred_comp_ri, aes(x = month_n, y = pred_prob_est),
            color = "black", size = 1.1)
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


