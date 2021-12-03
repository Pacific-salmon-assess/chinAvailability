### Random Splines Neg Binomial
## Version of neg_bin_explore constrained to stat areas 121/21 
## Nov. 12, 2021

library(tidyverse)
library(TMB)
library(sdmTMB)

# utility functions for prepping smooths 
source(here::here("R", "utils.R"))


dat <- readRDS(here::here("data", "rec", "month_area_recCatch.rds")) %>% 
  # identify months to exclude based on minimal catch or comp estimates 
  #(make region specific)
  mutate(
    min_m = case_when(
      region %in% c("NSoG", "SSoG") ~ 5,
      region %in% c("JdFS", "QCaJS") ~ 6
    ),
    max_m = case_when(
      region == "QCaJS" ~ 8, 
      region %in% c("JdFS", "NSoG", "SSoG")  ~ 9
    ),
    region = fct_relevel(region, "QCaJS", "JdFS", "NSoG", "SSoG")
  ) %>% 
  # drop areas with fewer than 10 datapoints (month-year-area observation = 1)
  group_by(area_n) %>% 
  filter(!month_n < min_m,
         !month_n > max_m) %>% 
  add_tally() %>% 
  ungroup() %>% 
  filter(!n < 10) %>%
  droplevels() 


## predicted dataset
new_dat <- group_split(dat, region) %>%
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
            dat %>% select(region, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, region, sep = "_"),
         offset = mean(dat$offset)) %>%
  # filter(strata %in% comp_strata) %>%
  arrange(region) %>%
  # convoluted to ensure ordering is correct for key
  mutate(
    month = as.factor(month_n),
    order = row_number(),
    # year is necessary regardless of whether RIs are inlcuded because of year-
    # specific smooths 
    reg_month_year = paste(region, month_n, year, sep = "_"),
    reg_month_year_f = fct_reorder(factor(reg_month_year), order)
  ) %>%
  select(-order, -reg_month_year, -strata) %>%
  distinct() %>% 
  # generate predictions only for swiftsure
  filter(area %in% c("121", "21"))



## PREP TMB --------------------------------------------------------------------

f1 <- catch ~ area +
  s(month_n, bs = "tp", k = 4, by = region) +
  s(month_n, by = year, bs = "tp", m = 1, k = 4) +
  offset


# generate model matrices
formula_no_sm <- remove_s_and_t2(f1)
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(f1, data = dat)
pred_X_ij <- predict(gam(formula_no_sm, data = dat), 
                     new_dat, type = "lpmatrix")
sm_pred <- parse_smoothers(f1, data = dat, newdata = new_dat)
rand_int_vec <- as.numeric(as.factor(as.character(dat$year))) - 1
pred_rand_int_vec <- as.numeric(new_dat$year) - 1
grouping_vec <- as.numeric(new_dat$reg_month_year_f) - 1
grouping_key <- unique(grouping_vec)


# offset
offset_pos <- grep("^offset$", colnames(X_ij))


# input data
data <- list(
  y1_i = dat$catch,
  X1_ij = X_ij,
  rfac1 = rand_int_vec,
  n_rfac1 = length(unique(rand_int_vec)),
  b_smooth_start = sm$b_smooth_start,
  Zs = sm$Zs, # optional smoother basis function matrices
  Xs = sm$Xs, # optional smoother linear effect matrix
  pred_X1_ij = pred_X_ij,
  pred_Zs = sm_pred$Zs,
  pred_Xs = sm_pred$Xs,
  pred_rfac1 = pred_rand_int_vec,
  pred_rfac_agg = grouping_vec,
  pred_rfac_agg_levels = grouping_key
)

# input parameter initial values
pars <- list(
  b1_j = rep(0, ncol(X_ij)),
  ln_phi = log(1.5),
  bs = rep(0, ncol(sm$Xs)),
  ln_smooth_sigma = rep(0, length(sm$sm_dims)),
  b_smooth = rep(0, sum(sm$sm_dims)),
  a1 = rep(0, length(unique(rand_int_vec))),
  ln_sigma_a1 = log(0.25)
)

# define offset
pars$b1_j[offset_pos] <- 1
b1_j_map <- seq_along(pars$b1_j)
b1_j_map[offset_pos] <- NA
tmb_map <- list(b1_j = as.factor(b1_j_map))

# define random parameters
tmb_random <- c("b_smooth", "a1")

compile(here::here("src", "negbin_rsplines_rint.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_rint")))

obj <- MakeADFun(data, pars, 
                 map = tmb_map,
                 random = tmb_random,
                 DLL = "negbin_rsplines_rint")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

saveRDS(ssdr, here::here("data", "model_fits", 
                         "negbin_rsmooths_121_21_only.rds"))


calc_ribbons <- function(pred_dat, fit_vec, se_vec) {
  pred_dat %>%
    mutate(link_fit = as.numeric(fit_vec),
           link_se_fit = as.numeric(se_vec),
           link_lo = link_fit + (qnorm(0.025) * link_se_fit),
           link_up = link_fit + (qnorm(0.975) * link_se_fit),
           fit = exp(link_fit),
           fit_lo = exp(link_lo),
           fit_up = exp(link_up))
}
tmb_pred_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Estimate"]
tmb_pred_se_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Std. Error"]
new_dat_region <- new_dat %>% 
  select(month_n, region, year, offset, month, reg_month_year_f) %>% 
  distinct()
pred_dat_tmb_region <- calc_ribbons(pred_dat = new_dat_region, 
                                    fit_vec = tmb_pred_region,
                                    se_vec = tmb_pred_se_region)

saveRDS(pred_dat_tmb_region, 
        here::here("data", "model_fits", 
                   "negbin_rsmooths_121_21_only_pred_abund.rds"))
