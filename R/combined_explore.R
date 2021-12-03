### Random Splines Neg Binomial + MVN Random Intercepts Dirichlet
## Combines neg_bin_explore and dirichlet_explore to generate new combined model
## Nov. 29, 2021

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)

# utility functions for prepping smooths 
source(here::here("R", "utils.R"))


# use example data from WCVI for initial fitting
# comp_ex <- stockseasonr::comp_ex %>%
#   #collapse some stock levels
#   mutate(
#     agg2 = case_when(
#       agg %in% c("PSD", "SOG", "FR-early", "FR-late", "ECVI") ~ "salish",
#       grepl("CR", agg) ~ "col",
#       TRUE ~ "other"
#     )
#   ) %>%
#   group_by(
#     sample_id, region, year, month_n, agg2, nn
#   ) %>%
#   summarize(agg_prob2 = sum(agg_prob), .groups = "drop") %>%
#   rename(agg = agg2, agg_prob = agg_prob2) %>%
#   ungroup()
# 
# catch <- stockseasonr::catch_ex %>% 
#   filter(
#     year %in% comp$year,
#     month_n %in% comp$month_n
#   )

comp <- readRDS(here::here("data", "rec", "coarse_rec_comp.rds")) %>% 
  filter(!region == "NSoG") %>% 
  droplevels()
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
  )

# subset predicted composition dataset to match pred_dat_catch since fitting 
# data can be more extensive
pred_dat_comp <- pred_dat_comp1 %>% 
  filter(reg_month_year %in% pred_dat_catch$reg_month_year)


## GENERATE TMB INPUTS ---------------------------------------------------------

# abundance inputs 
f1 <- catch ~ area +
  s(month_n, bs = "tp", k = 4, by = region) +
  s(month_n, by = year, bs = "tp", m = 1, k = 4) +
  offset
formula_no_sm <- remove_s_and_t2(f1)
X_ij <- model.matrix(formula_no_sm, data = catch)
sm <- parse_smoothers(f1, data = catch)
pred_X_ij <- predict(gam(formula_no_sm, data = catch), 
                     pred_dat_catch, type = "lpmatrix")
sm_pred <- parse_smoothers(f1, data = catch, newdata = pred_dat_catch)
rfac1 <- as.numeric(as.factor(as.character(catch$year))) - 1
pred_rfac1 <- as.numeric(pred_dat_catch$year) - 1
pred_rfac_agg <- as.numeric(pred_dat_catch$reg_month_year_f) - 1
pred_rfac_agg_levels <- unique(pred_rfac_agg)


# composition inputs 
comp_dat <- comp %>%
  pivot_wider(., names_from = agg, values_from = agg_prob) %>%
  mutate_if(is.numeric, ~ replace_na(., 0.00001)) %>% 
  mutate(dummy = 0)
obs_comp <- comp_dat %>%
  select(-c(sample_id:nn, dummy)) %>%
  as.matrix()

months2 <- unique(comp_dat$month_n)
spline_type <- if (max(months2) == 12) "cc" else "tp"
n_knots <- if (max(months2) == 12) 4 else 3
# response variable doesn't matter, since not fit
m2 <- mgcv::gam(rep(0, length.out = nrow(comp_dat)) ~
                  region + s(month_n, bs = spline_type, k = n_knots, 
                             by = region),
                data = comp_dat
)
X2_ij <- predict(m2, type = "lpmatrix")
pred_X2_ij <- predict(m2, pred_dat_comp, type = "lpmatrix")

rfac2 <- as.numeric(as.factor(as.character(comp_dat$year))) - 1
n_rint <- length(unique(rfac2))
pred_rfac2 <- as.numeric(pred_dat_comp$year) - 1


# input data
data <- list(
  #abundance
  y1_i = catch$catch,
  X1_ij = X_ij,
  rfac1 = rfac1,
  n_rfac1 = length(unique(rfac1)),
  b_smooth_start = sm$b_smooth_start,
  Zs = sm$Zs, # optional smoother basis function matrices
  Xs = sm$Xs, # optional smoother linear effect matrix
  pred_X1_ij = pred_X_ij,
  pred_Zs = sm_pred$Zs,
  pred_Xs = sm_pred$Xs,
  pred_rfac1 = pred_rfac1,
  pred_rfac_agg = pred_rfac_agg,
  pred_rfac_agg_levels = pred_rfac_agg_levels,
  #composition
  Y2_ik = obs_comp,
  X2_ij = X2_ij,
  rfac2 = rfac2,
  n_rfac2 = n_rint,
  pred_X2_ij = pred_X2_ij,
  pred_rfac2 = pred_rfac2
)

if(nrow(pred_X2_ij) != length(pred_rfac_agg_levels)) {
  warning("Prediction DFs do not match")
}

# input parameters
pars <- list(
  #abundance
  b1_j = rep(0, ncol(X_ij)),
  ln_phi = log(1.5),
  bs = rep(0, ncol(sm$Xs)),
  ln_smooth_sigma = rnorm(length(sm$sm_dims), 0, 0.5), #rep(0, length(sm$sm_dims)),
  b_smooth = rnorm(sum(sm$sm_dims), 0, 0.5),
  a1 = rnorm(length(unique(rfac1)), 0, 0.5), #rep(0, length(unique(rfac1))),
  ln_sigma_a1 = log(0.25),
  #composition
  B2_jk = matrix(0,
                 nrow = ncol(X2_ij),
                 ncol = ncol(obs_comp)
  ),
  # mvn matrix of REs
  A2_hk = matrix(rnorm(n_rint * ncol(obs_comp), 0, 0.5), #0,
                 nrow = n_rint,
                 ncol = ncol(obs_comp)
  ),
  ln_sigma_A2 = log(0.25)
)

# mapped parameters
tmb_map <- list()
offset_pos <- grep("^offset$", colnames(X_ij))
pars$b1_j[offset_pos] <- 1
b1_j_map <- seq_along(pars$b1_j)
b1_j_map[offset_pos] <- NA
tmb_map <- c(tmb_map, list(b1_j = as.factor(b1_j_map)))


## region-stock combinations with 0 observations to map
comp_map <- comp %>%
  group_by(region, agg) %>%
  complete(., region, nesting(agg)) %>%
  summarize(total_obs = sum(agg_prob), .groups = "drop") %>%
  filter(is.na(total_obs)) %>%
  select(agg, region)

if (!is.na(comp_map$agg[1])) {
  temp_betas <- pars$B2_jk
  for (i in seq_len(nrow(comp_map))) {
    offset_stock_pos <- grep(
      paste(comp_map$agg[1], collapse = "|"),
      colnames(obs_comp)
    )
    offset_region_pos <- grep(
      paste(comp_map$region[1], collapse = "|"),
      colnames(X_ij)
    )
    for (j in seq_len(ncol(X_ij))) {
      for (k in seq_len(ncol(obs_comp))) {
        if (j %in% offset_region_pos & k %in% offset_stock_pos) {
          temp_betas[j, k] <- NA
        }
      }
    }
  }
  comp_map_pos <- which(is.na(as.vector(temp_betas)))
  B2_jk_map <- seq_along(pars$B2_jk)
  B2_jk_map[comp_map_pos] <- NA
  tmb_map <- c(tmb_map, list(B2_jk = as.factor(B2_jk_map)))
}


# define random variables
tmb_random <- c("b_smooth", "a1", "A2_hk")


## FIT MODEL  ------------------------------------------------------------------

compile(here::here("src", "negbin_dirichlet_mvn_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_dirichlet_mvn_rsplines")))

# fit first w/out REs
tmb_map1 <- c(
  tmb_map,
  list(
    b_smooth = factor(rep(NA, length(pars$b_smooth))),
    ln_smooth_sigma = factor(rep(NA, length(pars$ln_smooth_sigma))),
    a1 = factor(rep(NA, length(pars$a1))),
    ln_sigma_a1 = as.factor(NA),
    A2_hk = factor(rep(NA, length(pars$A2_hk))),
    ln_sigma_A2 = as.factor(NA)
  )
)
obj1 <- TMB::MakeADFun(
  data = data,
  parameters = pars,
  map = tmb_map1,
  DLL = "negbin_dirichlet_mvn_rsplines"
)
opt1 <- stats::nlminb(obj1$par, obj1$fn, obj1$gr,
                      control = list(eval.max = 1e4, iter.max = 1e4)
)
sdr1 <- sdreport(obj1)

# fit with REs using first run as inits
obj <- TMB::MakeADFun(
  data = data,
  parameters = obj1$env$parList(opt1$par),
  map = tmb_map,
  random = tmb_random,
  DLL = "negbin_dirichlet_mvn_rsplines"
)
opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
                     )
nlminb_loops = 2
for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
  opt <- stats::nlminb(opt$par, obj$fn, obj$gr)
}

sdr <- sdreport(obj)
ssdr <- summary(sdr)

saveRDS(ssdr, here::here("data", "model_fits", "combined_mvn.rds"))


# GENERATE PREDICTIONS ---------------------------------------------------------

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
  ungroup() %>% 
  filter(!region == "QCaJS")

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
  geom_raster(aes(fill = fit_z)) +
  scale_fill_viridis_c(name = "Predicted\nAbundance\nAnomalies") +
  labs(x = "Month", y = "Year") +
  facet_grid(region~stock) +
  ggsidekick::theme_sleek()
