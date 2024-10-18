### MVN Random Intercepts Dirichlet
## Adapt stockseasonr model with MVN intercepts
## Nov. 29, 2021
## Updated March 19 to explore composition at area levels
## Updated April 14 to subarea levels

library(tidyverse)
library(TMB)
library(stockseasonr)

# tmb models - use MVN if time-varying predictions are required, use RI if 
# generating predictions for "average" year
# compile(here::here("src", "dirichlet_mvn.cpp"))
# dyn.load(dynlib(here::here("src", "dirichlet_mvn")))
# compile(here::here("src", "dirichlet_ri.cpp"))
# dyn.load(dynlib(here::here("src", "dirichlet_ri")))
# compile(here::here("src", "dirichlet_ri_sdmTMB.cpp"))
# dyn.load(dynlib(here::here("src", "dirichlet_ri_sdmTMB")))


# utility functions for prepping smooths 
# source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit_new.R"))


# pre-cleaning: aggregate at PST, remove sublegals
size_in <- readRDS(here::here("data", "rec", "rec_size.rds"))
comp_in <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  rename(stock_region = region) %>% 
  mutate(
    agg_new = case_when(
      grepl("CR", pst_agg) ~ "WA_OR_CA",
      grepl("CST", pst_agg) ~ "WA_OR_CA",
      pst_agg == "NBC_SEAK" ~ "SOG",
      Region1Name %in% c("Fraser_Summer_4.1") ~ "Fraser_S",
      Region1Name %in% c("Fraser_Summer_5.2", "Fraser_Spring_5.2",
                         "Fraser_Spring_4.2") | pst_agg == "FR-early" ~ "Fraser_Yearling",
      Region1Name %in% c("Fraser_Fall") ~ "Fraser_F",
      TRUE ~ pst_agg
    )
  )

# trim size and gsi data
trim_foo <- function(dat_in) {
  dat_in %>%
    filter(!legal == "sublegal"#,
           #remove subareas 19A, 18A w/ very small sample sizes 
           # (and outside study area)
           # !subarea %in% c("19A", "18A")
           ) %>% 
    rename(region = cap_region) %>% 
    mutate(
      month_n = lubridate::month(date),
      subarea_original = subarea,
      subarea = case_when(
        subarea == "121C" ~ "121A", #based on fihsing location
        subarea %in% c("18B", "18D", "18E") ~ "18BDE",
        subarea %in% c("21B") ~ "21A",
        subarea %in% c("US7") ~ "19B",
        # consolidate subareas due to small sample sizes and doc issues w/
        # convergence
        subarea %in% c("19C", "19D", "19E") ~ "19CDE",
        subarea %in% c("20B", "20C") ~ "20BC",
        subarea %in% c("20A", "20E") ~ "20AE",
        subarea %in% c("20DB", "20DI", "20DO") ~ "20D",
        subarea %in% c("29G", "29F") ~ "29FG",
        # subarea == "29E" ~ "29DE",
        grepl("29D", subarea) ~ "29D",
        subarea %in% c("29-11", "29I") ~ "29D", 
        subarea %in% c("29B", "29C") ~ "29BC",
        TRUE ~ subarea
      ),
      reg = case_when(
        subarea == "21A" ~ "SWVI",
        area %in% c("121", "21") ~ "SWVI",
        subarea %in% c("19D", "19E") ~ "JdFS",
        area == "20" ~ "JdFS",
        area %in% c("18", "19", "29") ~ "SSoG",
        TRUE ~ "out"
      ),
      reg = as.factor(reg),
      core_area = case_when(
        subarea == "121B" ~ "no",
        subarea == "29DE" ~ "yes",
        reg %in% c("SWVI") ~ "yes", 
        subarea %in% c("18BDE", "19B", "19C", "20AE") ~ "yes",
        TRUE ~ "no"
      ),
      week = lubridate::week(date),
      yday = lubridate::yday(date),
      month_n = lubridate::month(date),
      sample_id = paste(month_n, subarea, week, year, sep = "_"),
      area_f = as.factor(as.numeric(gsub("([0-9]+).*$", "\\1", area)))
    ) %>% 
    group_by(sample_id) %>% 
    mutate(nn = length(unique(id)) %>% as.numeric) %>% 
    ungroup()
} 

trim_size <- trim_foo(size_in) %>% 
  group_by(sample_id, subarea, subarea_original, area_f, reg, reg_c = region, 
           week, month, month_n, year, nn, size_bin, core_area) %>% 
  summarize(prob = length(unique(id)), .groups = "drop")  %>% 
  ungroup() %>% 
  filter(
    month_n < 10.1 & month_n > 1.9,
    !reg == "out" 
  ) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg),
         subarea = as.factor(subarea),
         prob = as.numeric(prob)
  ) %>% 
  droplevels() 

trim_stock <- comp_in %>%
  trim_foo() %>% 
  group_by(sample_id, subarea, subarea_original, area_f, reg, reg_c = region, 
           week, month, month_n, year, nn, agg_new, core_area) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  filter(
    month_n < 10.1 & month_n > 1.9,
    !reg == "out" 
  ) %>% 
  mutate(year = as.factor(year),
         reg = as.factor(reg),
         subarea = as.factor(subarea)
  ) %>% 
  droplevels() 



# look at sample coverage in data passed to model
alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")

trim_size %>%
  select(sample_id, year, reg, subarea, subarea_original, month_n,
         nn, core_area) %>%
  distinct() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg,
                  shape = core_area),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(subarea, as.numeric(reg))) +
  ggsidekick::theme_sleek()


# add areas to composition dataset
area_key <- trim_stock %>% 
  select(subarea, subarea_original, area_f, reg, core_area) %>% 
  distinct()
# saveRDS(area_key, here::here("data", "rec", "subarea_key.RDS"))

# subset predicted composition dataset
pred_dat_comp <- expand.grid(
    reg = unique(trim_stock$reg),
    month_n = seq(min(trim_stock$month_n),
                  max(trim_stock$month_n),
                  by = 0.1
    )
  ) %>% 
  left_join(., area_key, by = "reg") %>%
  filter(#core_area == "yes",
    month_n < 9.1 & month_n > 4.9
  ) %>% 
  select(-core_area) %>% 
  mutate(week = month_n * 4) 


## FIT MODELS ------------------------------------------------------------------

library(stockseasonr)

# stock_mod_ri <- fit_stockseasonr(
#   comp_formula = agg_new ~ 1 + subarea + 
#     s(month_n, bs = "tp", k = 4, m = 2) + (1 | year),
#   comp_dat = trim_stock,
#   pred_dat = pred_dat_comp,
#   model = "dirichlet",
#   random_walk = TRUE,
#   fit = TRUE,
#   nlminb_loops = 2, newton_loops = 1
# )
# 
# ssdr <- stock_mod_ri$ssdr 
# 
# beta_mat <- ssdr[rownames(ssdr) == "B2_jk", 2] %>% 
#   matrix(., 
#             nrow = ncol(stock_mod_ri$tmb_data$X2_ij),
#             ncol = ncol(stock_mod_ri$tmb_data$Y2_ik))
# rownames(beta_mat) <- colnames(stock_mod_ri$tmb_data$X2_ij)
# colnames(beta_mat) <- colnames(stock_mod_ri$tmb_data$Y2_ik)
# beta_mat
# 
# dum <- model_inputs_ri$tmb_data$pred_X2_ij %*% beta_mat
# 
# saveRDS(stock_mod_ri$ssdr, 
#         here::here("data", "model_fits", "subarea",
#                    "dirichlet_ri_mig_corridor.rds"))


# as above but with week as a covariate
stock_mod_ri2 <- fit_stockseasonr(
  comp_formula = agg_new ~ 1 + subarea + 
    s(week, bs = "tp", k = 4, m = 2) + (1 | year),
  comp_dat = trim_stock,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = TRUE,
  fit = TRUE,
  nlminb_loops = 2, newton_loops = 1
)

ssdr2 <- stock_mod_ri2$ssdr 

beta_mat2 <- ssdr[rownames(ssdr2) == "B2_jk", 2] %>% 
  matrix(., 
         nrow = ncol(stock_mod_ri2$tmb_data$X2_ij),
         ncol = ncol(stock_mod_ri2$tmb_data$Y2_ik))
rownames(beta_mat2) <- colnames(stock_mod_ri2$tmb_data$X2_ij)
colnames(beta_mat2) <- colnames(stock_mod_ri2$tmb_data$Y2_ik)
beta_mat2

saveRDS(stock_mod_ri2$ssdr, 
        here::here("data", "model_fits", "subarea",
                   "dirichlet_ri_mig_corridor_week.rds"))



## EVALUATE MODEL PREDS --------------------------------------------------------

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

stock_seq <- colnames(stock_mod_ri$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  mutate(
    area_f = fct_reorder(area_f, as.numeric(reg)),
    stock = fct_relevel(
      stock, "Fraser_Yearling", "Fraser_S", "Fraser_F",
      "SOG", "PSD", "WCVI", "WA_OR_CA"
    )
  )


p <- ggplot(data = pred_comp_comb, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(subarea~stock) +
  # facet_grid(subarea~stock) +
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



# SIZE COMP MODEL --------------------------------------------------------------

# as above but with week as a covariate
size_mod_ri <- fit_stockseasonr(
  comp_formula = size_bin ~ 1 + subarea + 
    s(week, bs = "tp", k = 4, m = 2) + (1 | year),
  comp_dat = trim_size,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = TRUE,
  fit = TRUE,
  nlminb_loops = 2, newton_loops = 1
)

ssdr_size <- size_mod_ri$ssdr 

beta_mat <- ssdr_size[rownames(ssdr_size) == "B2_jk", 2] %>% 
  matrix(., 
         nrow = ncol(size_mod_ri$tmb_data$X2_ij),
         ncol = ncol(size_mod_ri$tmb_data$Y2_ik))
rownames(beta_mat) <- colnames(size_mod_ri$tmb_data$X2_ij)
colnames(beta_mat) <- colnames(size_mod_ri$tmb_data$Y2_ik)
beta_mat

saveRDS(size_mod_ri$ssdr, 
        here::here("data", "model_fits", "subarea",
                   "size_week.rds"))



logit_pred_ppn <- ssdr_size[rownames(ssdr_size) == "logit_pred_Pi_prop", ]
pred_mu <- ssdr_size[rownames(ssdr_size) == "pred_Mu2", ]

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

size_seq <- colnames(size_mod_ri$tmb_data$Y2_ik)
pred_comp_size <- purrr::map(size_seq, function (x) {
  dum <- pred_dat_comp
  dum$size_bin <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  mutate(
    area_f = fct_reorder(area_f, as.numeric(reg)),
    size_bin = fct_relevel(
      size_bin, "<45", "45-60", "60-75", "75-90", ">90"
    )
  )


ggplot(data = pred_comp_size, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(area_f~size_bin) +
  # facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, colour = subarea))
