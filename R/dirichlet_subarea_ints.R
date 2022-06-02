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
         #remove subareas 19A, 18A w/ very small sample sizes 
         # (and outside study area)
         !subarea %in% c("19A", "18A")) %>% 
  rename(stock_region = region, region = cap_region) %>% 
  mutate(month_n = lubridate::month(date),
         subarea_original = subarea,
         subarea = case_when(
           subarea %in% c("18B", "18D", "18E") ~ "18BDE",
           subarea %in% c("21B") ~ "21A",
           subarea %in% c("US7") ~ "19B",
           # consolidate subareas due to small sample sizes and doc issues w/
           # convergence
           subarea %in% c("19C", "19D", "19E") ~ "19CDE",
           subarea %in% c("20C", "20D", "20E") ~ "20CDE",
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
           subarea %in% c("18BDE", "19B", "19C") ~ "yes",
           TRUE ~ "no"
         ),
         week = lubridate::week(date),
         yday = lubridate::yday(date),
         month_n = lubridate::month(date),
         # sample_id = paste(month_n, reg, week, year, sep = "_"),
         sample_id = paste(month_n, reg, yday, year, sep = "_"),
         can_reg = case_when(
           pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa", "CA_ORCST", 
                          "CR-lower_sp", "CR-upper_sp", "PSD", 
                          "NBC_SEAK", "WACST") ~ "Other",
           Region1Name == "SOMN" ~ "ECVI",
           # Region1Name %in% c("Fraser_Summer_5.2", "Fraser_Spring_5.2",
           #                    "Fraser_Spring_4.2") ~ "Fraser_Yearling",
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
  mutate(nn = length(unique(id))) %>% 
  ungroup()

stock_comp <- comp1 %>%  
  group_by(sample_id, subarea, subarea_original, area, reg, reg_c = region, 
           week, month, month_n, year, nn, can_reg, core_area) %>% 
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

png(here::here("figs", "data_coverage", "comp_model_inputs.png"), height = 5,
    width = 5, units = "in", res = 200)
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


# prediction datasets 
pred_dat_comp1 <- group_split(stock_comp, reg) %>%
  map_dfr(., function(x) {
    expand.grid(
      reg = unique(x$reg),
      # week = unique(x$week),
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
# month_key <- stock_comp %>% 
#   select(week, month_n) %>% 
#   distinct()

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


# stock_comp %>% 
#   group_by(subarea, reg, month, can_reg) %>% 
#   summarize(sum_prob = prob) %>% 
#   ggplot(.) +
#   geom_col(aes(x = month, y = sum_prob, fill = reg)) +
#   facet_grid(can_reg~subarea, scales = "free_y")
# 

## FIT MODEL -------------------------------------------------------------------

# rand predictions
# model_inputs <- make_inputs(
#   comp_formula = pst_agg ~ subarea + 
#     s(month_n, bs = "tp", k = 3, by = reg, m = 2),
#   comp_dat = stock_comp,
#   comp_rint = "year",
#   pred_dat = pred_dat_stock_comp,
#   model = "dirichlet",
#   include_re_preds = TRUE
# )
# 
# stock_mod <- fit_model(
#   tmb_data = model_inputs$tmb_data, 
#   tmb_pars = model_inputs$tmb_pars, 
#   tmb_map = model_inputs$tmb_map, 
#   tmb_random  = model_inputs$tmb_random,
#   fit_random = TRUE,
#   ignore_fix = TRUE,
#   model_specs = model_inputs$model_specs
#   )
# 
# 
# saveRDS(stock_mod$ssdr, 
#         here::here("data", "model_fits", "subarea", 
#                    "dirichlet_mvn_mig_corridor.rds"))


# no rand predictions
model_inputs_ri <- make_inputs(
  comp_formula = can_reg ~ subarea + 
    # s(week, bs = "tp", k = 3, m = 2) +
    s(month_n, bs = "cc", k = 4, m = 2, by = reg)
    ,
  comp_dat = stock_comp,
  comp_rint = "year",
  pred_dat = pred_dat_stock_comp_ri,
  model = "dirichlet",
  include_re_preds = FALSE
)

stock_mod_ri <- fit_model(
  tmb_data = model_inputs_ri$tmb_data, 
  tmb_pars = model_inputs_ri$tmb_pars, 
  tmb_map = model_inputs_ri$tmb_map, 
  tmb_random  = model_inputs_ri$tmb_random,
  fit_random = TRUE,
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

dum <- model_inputs_ri$tmb_data$pred_X2_ij %*% beta_mat

saveRDS(stock_mod_ri$ssdr, 
        here::here("data", "model_fits", "subarea",
                   "dirichlet_ri_mig_corridor.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------


# random effects predictions 
# ssdr_ri <- readRDS(
#   here::here("data", "model_fits", "subarea",
#              "dirichlet_mvn_mig_corridor.rds"))
# fixed effects predictions 
ssdr <- readRDS(
  here::here("data", "model_fits", "subarea", 
             "dirichlet_ri_mig_corridor.rds"))


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
  dum <- pred_dat_stock_comp_ri
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

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) #+

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)



## compare to observations
# number of samples in a daily observation
long_dat <- model_inputs_ri$wide_comp_dat %>%
  mutate(samp_nn = apply(model_inputs_ri$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(Other:WCVI), 
               names_to = "stock", 
               values_to = "obs_count") %>%
  mutate(obs_ppn = obs_count / samp_nn) %>% 
  filter(subarea %in% pred_comp$subarea,
         month_n %in% stock_comp$month_n
         ) 

p_obs <- ggplot() +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_jitter(data = long_dat,
              aes(x = month_n, y = obs_ppn, size = samp_nn), alpha = 0.2) +
  geom_line(data = pred_comp, aes(x = month_n, y = pred_prob_est), 
            colour = "red") +
  scale_size_continuous() +
  theme(legend.position = "top")


# stacked bar plots
pred_bar <- pred_comp %>% 
  filter(month_n %in% c("6", "7", "8", "9")) %>% 
  ggplot(.) +
  geom_col(aes(x = month_n, y = pred_prob_est, fill = stock)) +
  facet_wrap(~subarea, scales = "free_y") +
  scale_fill_brewer(type = "div", palette = 9) +
  theme(legend.position = "top")


pdf(here::here("figs", "jdf_subarea_preds", "subarea_month.pdf"))
p
p_ribbon
p_obs
pred_bar
dev.off()
