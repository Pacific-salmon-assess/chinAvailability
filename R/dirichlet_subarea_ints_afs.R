### Stock Composition Estimates
## Extension of dirichlet_subarea_ints.R for AFS presentation (predictions
# for all subareas, using PST aggregates)
## Aug 4, 2022

library(tidyverse)
library(TMB)


# tmb models - use MVN if time-varying predictions are required, use RI if 
compile(here::here("src", "dirichlet_ri_sdmTMB.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_ri_sdmTMB")))


# utility functions for prepping smooths 
# source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit_new.R"))


# pre-cleaning: aggregate at PST, remove sublegals
comp1 <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>%
  filter(!legal == "sublegal",
         #remove subareas 19A, 18A w/ very small sample sizes 
         # (and outside study area)
         !subarea %in% c("19A", "18A", "29I", "29-11")) %>% 
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
           subarea %in% c("20A", "20E") ~ "20AE",
           subarea %in% c("20C", "20D") ~ "20CD",
           subarea %in% c("29G", "29F") ~ "29FG",
           subarea == "29E" ~ "29DE",
           grepl("29D", subarea) ~ "29DE",
           subarea %in% c("29B", "29C") ~ "29BC",
           TRUE ~ subarea
         ),
         reg = case_when(
           subarea == "21A" ~ "SWVI",
           area %in% c("121", "21") ~ "SWVI",
           subarea == "19C" ~ "SSoG",
           area %in% c("20W", "20E", "19JDF") ~ "JdFS",
           area %in% c("18", "19GST", "29") ~ "SSoG",
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
         sample_id = paste(month_n, reg, yday, year, sep = "_"),
         agg_new = case_when(
           grepl("CR", pst_agg) ~ "WA_OR_CA",
           grepl("CST", pst_agg) ~ "WA_OR_CA",
           # pst_agg %in% c("CR-upper_su/fa", "CR-lower_fa") ~ "Col_SF",
           # pst_agg %in% c("CR-lower_sp", "CR-upper_sp") ~ "Col_Sp",
           pst_agg == "NBC_SEAK" ~ "SOG",
           Region1Name %in% c("Fraser_Summer_4.1") ~ "Fraser_S",
           Region1Name %in% c("Fraser_Summer_5.2", "Fraser_Spring_5.2",
                              "Fraser_Spring_4.2") | pst_agg == "FR-early" ~ "Fraser_Yearling",
           Region1Name %in% c("Fraser_Fall") ~ "Fraser_F",
           TRUE ~ pst_agg
         ),
         # simplify area assignments 
         area_f = as.factor(as.numeric(gsub("([0-9]+).*$", "\\1", area)))
  ) %>% 
  group_by(sample_id) %>% 
  mutate(nn = length(unique(id)) %>% as.numeric) %>% 
  ungroup() %>% 
  filter(!reg == "out")

# calculate mean distance after consolidating subareas
mean_dist <- comp1 %>% 
  select(subarea, dist_123i) %>% 
  distinct() %>% 
  group_by(subarea) %>% 
  summarize(dist_123i = mean(dist_123i)) %>% 
  ungroup()

stock_comp <- comp1 %>%  
  group_by(sample_id, subarea, subarea_original, area_f, reg, 
           reg_c = region, 
           week, month, month_n, year, nn, agg_new, core_area) %>% 
  summarize(prob = sum(prob), 
            .groups = "drop") %>% 
  ungroup() %>% 
  filter(month_n < 10.1 & month_n > 1.9
         ) %>% 
  mutate(year = as.factor(year),
         reg = factor(reg, levels = c("SWVI", "JdFS", "SSoG")),
         subarea = as.factor(subarea)
         ) %>% 
  droplevels() %>% 
  left_join(., mean_dist, by = "subarea")


# look at sample coverage in data passed to model
alpha_scale <- c(0.3, 0.95)
names(alpha_scale) <- c("no", "yes")
stock_comp %>%
  select(sample_id, year, reg, subarea, month_n, nn, core_area) %>%
  distinct() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = nn, colour = reg,
                  shape = core_area),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(subarea, as.numeric(reg))) +
  ggsidekick::theme_sleek()


# # prediction datasets 
# pred_dat_comp1 <- group_split(stock_comp, reg) %>%
#   map_dfr(., function(x) {
#     expand.grid(
#       reg = unique(x$reg),
#       month_n = seq(min(x$month_n),
#                     max(x$month_n),
#                     by = 0.1
#       )
#     ) 
#   }) 

# add areas to composition dataset
area_key <- stock_comp %>% 
  select(subarea, subarea_original, area_f, reg, core_area) %>% 
  distinct()
# saveRDS(area_key, here::here("data", "rec", "subarea_key.RDS"))


# subset predicted composition dataset
pred_dat_comp <- #pred_dat_comp1 %>% 
  expand.grid(
    reg = unique(stock_comp$reg),
    month_n = seq(min(stock_comp$month_n),
                  max(stock_comp$month_n),
                  by = 0.1
    )
  ) %>% 
  left_join(., area_key, by = "reg") %>%
  filter(#core_area == "yes",
         month_n < 9.1 & month_n > 4.9
         ) %>% 
  select(-core_area) %>% 
  mutate(week = month_n * 4) 



## FIT MODEL -------------------------------------------------------------------



fit1 <- fit_stockseasonr(
  comp_formula = agg_new ~ 1 + area_f +
    # s(month_n, bs = "tp", k = 4, m = 2) +
    # (1 | reg) +
    (1 | year),
  comp_dat = stock_comp,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = TRUE,
  fit = FALSE,
  nlminb_loops = 2, newton_loops = 1
)
re2 <- fit1$ssdr[rownames(fit1$ssdr) == "re2", ]

# 
# n2 <- nrow(stock_comp)
# n_re2 <- ncol(fit1$tmb_data$re_index2)
# re2 <- fit1$tmb_pars$re2
# re_index2 <- fit1$tmb_data$re_index2
# nobs_re2 <- fit1$tmb_data$nobs_re2
# 
# mu2_i <- eta_re2_i <- rep(0, length.out = n2)
# for (i in seq_len(n2)) {
#   temp <- 0
#   for (g in seq_len(n_re2)) {
#     if (g == 1) eta_re2_i[i] <- re2[re_index2[i, g]]
#     if (g > 1) {
#       dum <- temp
#       temp <- nobs_re2[g - 1] + dum
#       eta_re2_i[i] <- re2[re_index2[i, g] + temp]
#     } 
#   }
# }


# tt <- sdmTMB_dummy$tmb_data$RE_indexes 


source(here::here("R", "functions", "fit_new.R"))

fit1 <- fit_stockseasonr(
  comp_formula = agg_new ~ 1 + area_f + 
    s(month_n, bs = "tp", k = 4, m = 2) + 
    (1 | year),
  comp_dat = stock_comp,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  fit = TRUE,
  nlminb_loops = 2, newton_loops = 1
)
# fails to converge with area included too
fit2 <- fit_stockseasonr(
  comp_formula = agg_new ~ 1 + #area_f +
    s(dist_123i, bs = "tp", k = 3, m = 2) +
    s(month_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = stock_comp,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  fit = TRUE,
  nlminb_loops = 2, newton_loops = 1
)
fit3 <- fit_stockseasonr(
  comp_formula = agg_new ~ 1 + subarea + 
    s(month_n, bs = "tp", k = 4, m = 2) + 
    (1 | year),
  comp_dat = stock_comp,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  fit = TRUE,
  nlminb_loops = 2, newton_loops = 1
)

ssdr1 <- fit1$ssdr 
ssdr2 <- fit2$ssdr 
ssdr3 <- fit3$ssdr 

beta_mat <- ssdr3[rownames(ssdr3) == "B2_jk", 2] %>% 
  matrix(., 
            nrow = ncol(fit3$tmb_data$X2_ij),
            ncol = ncol(fit3$tmb_data$Y2_ik))
rownames(beta_mat) <- colnames(fit3$tmb_data$X2_ij)
colnames(beta_mat) <- colnames(fit3$tmb_data$Y2_ik)
beta_mat

dum <- model_inputs_ri$tmb_data$pred_X2_ij %*% beta_mat

saveRDS(fit1, 
        here::here("data", "model_fits", "subarea",
                   "afs_fit1.rds"))
saveRDS(fit2, 
        here::here("data", "model_fits", "subarea",
                   "afs_fit2.rds"))
saveRDS(fit3, 
        here::here("data", "model_fits", "subarea",
                   "afs_fit3.rds"))


## EVALUATE MODEL PREDS --------------------------------------------------------

ssdr <- fit3$ssdr

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

stock_seq <- colnames(fit3$tmb_data$Y2_ik)
pred_comp <- purrr::map(stock_seq, function (x) {
  dum <- pred_dat_comp
  dum$stock <- x
  return(dum)
}) %>%
  bind_rows() %>%
  cbind(., link_preds) %>% 
  mutate(
    area_f = fct_reorder(as.factor(area_f), as.numeric(reg)),
    stock = fct_relevel(
      stock, "Fraser_Yearling", "Fraser_S", "Fraser_F",
      "SOG", "PSD", "WCVI", "WA_OR_CA"
    )
  ) %>% 
  select(-subarea_original) %>% 
  distinct()

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  # facet_grid(area~stock) +
  facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) #+

p_ribbon <- p +
  geom_ribbon(data = pred_comp,
              aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)



## compare to observations
# number of samples in a daily observation
long_dat <- fit3$wide_comp_dat %>%
  mutate(samp_nn = apply(fit3$tmb_data$Y2_ik, 1, sum)) %>%
  pivot_longer(cols = c(Fraser_S:Fraser_Yearling), 
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


# stacked ribbon plots
# subset areas to focus on core only
cores <- area_key %>% filter(core_area == "yes") %>% pull(subarea)

stack_comp <- pred_comp %>% filter(subarea %in% cores) %>% 
  mutate(subarea = fct_relevel(
    subarea, "21A", "121A", "20AE", "19B", "18BDE", "29DE"
  )) %>% 
  ggplot(.) +
  geom_area(aes(x = month_n, y = pred_prob_est, colour = stock, fill = stock), 
            stat = "identity") + 
  facet_wrap(~subarea, nrow = 1) +
  scale_colour_brewer(type = "div", palette = 1) +
  scale_fill_brewer(type = "div", palette = 1) +
  ggsidekick::theme_sleek() +
  labs(x = "Month", y = "Predicted Stock Composition") +
  theme(legend.position = "right") 

stack_comp_split <- split(pred_comp, pred_comp$subarea) %>% 
  purrr::map2(., names(.), function (x, y) {
    ggplot(x) +
      geom_area(aes(x = month_n, y = pred_prob_est, colour = stock, fill = stock), 
                stat = "identity") + 
      scale_colour_brewer(type = "div", palette = 9) +
      scale_fill_brewer(type = "div", palette = 9) +
      ggsidekick::theme_sleek() +
      labs(x = "Month", y = "Predicted Stock Composition", title = y) 
  })


pdf(here::here("figs", "afs_subarea_preds", "subarea_month.pdf"))
p
p_ribbon
p_obs
dev.off()

png(here::here("figs", "afs_subarea_preds", "stacked_ribbon_facet.png"),
    height = 2.5, width = 9, units = "in", res = 250)
stack_comp
dev.off()

pdf(here::here("figs", "afs_subarea_preds", "stacked_ribbon.pdf"),
    height = 3.5, width = 4.5)
stack_comp
stack_comp_split
dev.off()
