### MVN Random Intercepts Dirichlet
## Adjust subarea model to estimate changes in composition as a function of 
## migration distance; might allow more complex models to be fit

library(tidyverse)
library(TMB)

compile(here::here("src", "dirichlet_ri.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_ri")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))


# CLEAN ------------------------------------------------------------------------

# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal"#,
         #remove subareas 19A, 18A w/ very small sample sizes 
         # (and outside study area)
         # !subarea %in% c("19A", "18A")
         ) %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date),
         subarea_original = subarea,
        reg = case_when(
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
           subarea %in% c("18B", "18D", "18E", "19B", "19C") ~ "yes",
           TRUE ~ "no"
         ),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         # sample_id = paste(month_n, reg, week, year, sep = "_"),
         sample_id = paste(month_n, subarea, yday, year, sep = "_"),
         can_reg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                          "CR-lower_sp", "CR-upper_sp", "PSD", 
                          "NBC_SEAK", "WACST") ~ "Other",
           Region1Name == "SOMN" ~ "ECVI",
           TRUE ~ Region1Name
         ),
         pst_agg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                          "CR-lower_sp", "CR-upper_sp", "PSD", 
                          "NBC_SEAK", "WACST") ~ "other",
           TRUE ~ pst_agg
         )
  ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id)) %>% as.numeric) %>% 
  ungroup()

stock_comp <- comp1 %>%  
  group_by(sample_id, subarea, creelsub, area, reg, reg_c = region, 
           week, month, month_n, year, nn, can_reg, core_area, dist_123i) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  droplevels() %>% 
  filter(!reg == "out"#,
         # week > 22 & week < 38
         # month_n < 9.1 & month_n > 1.9
  ) %>% 
  mutate(year = as.factor(year),
         reg = factor(reg, levels = c("SWVI", "JdFS", "SSoG")),
         subarea = as.factor(subarea),
         area = as.factor(area)
  )


# look at sample coverage in data passed to model
alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")

png(here::here("figs", "data_coverage", "comp_model_inputs_distance.png"), 
    height = 5, width = 7, units = "in", res = 200)
stock_comp %>% 
  select(sample_id, year, reg, subarea, month_n, nn, core_area) %>% 
  distinct() %>% 
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg, 
                  shape = core_area),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(subarea, as.numeric(reg))) +
  ggsidekick::theme_sleek()
dev.off()

hist(stock_comp$dist_123i)


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
  select(subarea, area, reg, core_area, dist_123i) %>% 
  distinct()

# subset predicted composition dataset
pred_dat_stock_comp <- pred_dat_comp1 %>% 
  left_join(., area_key, by = "reg") %>%
  filter(core_area == "yes",
         month_n < 9.1 & month_n > 5.9#,
         # week < 36 & week > 22
  ) %>% 
  select(-core_area)

# remove random effects from predictions to generate estimates for "average"
# year
pred_dat_stock_comp_ri <- pred_dat_stock_comp %>%
  select(-year) %>%
  distinct()


## FIT MODEL -------------------------------------------------------------------

# no rand predictions
model_inputs_ri <- make_inputs(
  comp_formula = can_reg ~ s(dist_123i, k = 5) #+ 
    # s(week, bs = "tp", k = 3, m = 2) +
    # s(month_n, bs = "cc", k = 4, m = 2, by = reg)
  ,
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_dat = pred_dat_stock_comp_ri %>% select(subarea, dist_123i) %>% distinct(),
  model = "dirichlet",
  include_re_preds = FALSE
)

stock_mod_ri <- fit_model(
  tmb_data = model_inputs_ri$tmb_data, 
  tmb_pars = model_inputs_ri$tmb_pars, 
  tmb_map = model_inputs_ri$tmb_map, 
  tmb_random  = model_inputs_ri$tmb_random,
  fit_random = FALSE,
  ignore_fix = FALSE,
  model_specs = model_inputs_ri$model_specs
)

ssdr <- stock_mod_ri$ssdr 
beta_mat <- ssdr[rownames(ssdr) == "B2_jk", 2] %>% 
  matrix(., 
         nrow = ncol(model_inputs_ri$tmb_data$X2_ij),
         ncol = ncol(model_inputs_ri$tmb_data$Y2_ik))
rownames(beta_mat) <- colnames(model_inputs_ri$tmb_data$X2_ij)
colnames(beta_mat) <- colnames(model_inputs_ri$tmb_data$Y2_ik)
beta_mat


## EXPLORE PREDICTIONS ---------------------------------------------------------

logit_pred_ppn <- ssdr[rownames(ssdr) == "logit_pred_Pi_prop", ]
pred_mu <- ssdr[rownames(ssdr) == "pred_Mu2", ]

link_preds <- data.frame(
  link_prob_est = logit_pred_ppn[ , "Estimate"],
  link_prob_se =  logit_pred_ppn[ , "Std. Error"],
  pred_mu = pred_mu[ , "Estimate"]
) %>% 
  mutate(
    pred_prob_est = plogis(link_prob_est),
    pred_prob_low = plogis(link_prob_est + (qnorm(0.025) * link_prob_se)),
    pred_prob_up = plogis(link_prob_est + (qnorm(0.975) * link_prob_se))
  ) 

stock_seq <- colnames(model_inputs_ri$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_stock_comp_ri %>% 
    select(reg, area, subarea, dist_123i) %>% 
    distinct()
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  mutate(
    area_f = fct_reorder(as.factor(area), as.numeric(reg)),
    stock = fct_relevel(
      stock, "Fraser_Spring_4.2", "Fraser_Spring_5.2", "Fraser_Summer_5.2",
      "Fraser_Summer_4.1", "Fraser_Fall", "ECVI", #"SOMN",
      "WCVI", "Other"
    )
  )

ggplot(data = pred_comp, aes(x = dist_123i)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_wrap(~stock) +
  ggsidekick::theme_sleek() +
  geom_point(aes(y = pred_prob_est)) #+

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)