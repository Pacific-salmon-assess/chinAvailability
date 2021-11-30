### Random Splines Neg Binomial
## Adapt stockseasonr model to account for varying splines among REs 
## Nov. 12, 2021

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)

# utility functions for prepping smooths 
source(here::here("R", "utils.R"))


rec_catch <- readRDS(here::here("data", "rec", "month_area_recCatch.rds")) %>% 
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
  droplevels() #%>% 
  # group by region for now
  # group_by(region, region_c, month, month_n, year) %>%
  # summarize(catch = sum(catch),
  #           eff = sum(eff),
  #           offset = log(eff)) %>%
  # ungroup()


## predicted dataset
pred_dat <- group_split(rec_catch, region) %>%
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
            rec_catch %>% select(region, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, region, sep = "_"),
         offset = mean(rec_catch$offset)) %>%
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
  distinct() 



## MGCV VS TMB -----------------------------------------------------------------


dat <- rec_catch %>% 
  # filter(
  #   region %in% c("JdFS", "SSoG"),
  #   year %in% c("2015", "2016", "2017", "2018", "2019")
  # ) %>%
  droplevels() 
new_dat <- pred_dat %>% 
  filter(
    area %in% dat$area,
    year %in% dat$year
  ) %>% 
  droplevels()


f1a <- catch ~ area +
  s(month_n, bs = "tp", k = 3, by = region) +
  s(month_n, by = year, bs = "tp", m = 1, k = 3) +
  #smooth REs can't be readily passed with Anderson's function; code RIs instead
  s(year, bs = "re") + 
  offset(offset)
f1 <- catch ~ area +
  s(month_n, bs = "tp", k = 3, by = region) +
  s(month_n, by = year, bs = "tp", m = 1, k = 3) +
  offset


# fit mgcv for comparison
m1 <- mgcv::gam(f1a,
                data = dat,
                family = mgcv::nb()
)

preds_m1 <- predict.gam(m1,
                        new_dat,
                        # exclude = "s(year)",
                        se.fit = TRUE)

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

pred_dat_m1 <- calc_ribbons(pred_dat = new_dat, fit_vec = preds_m1$fit,
                            se_vec = preds_m1$se.fit)

ggplot() +
  geom_line(data = pred_dat_m1, aes(x = month_n, y = fit, colour = year)) +
  geom_ribbon(data = pred_dat_m1,
              aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year,
                  colour = year),
              alpha = 0.4) +
  facet_wrap(~area, scales = "free_y") +
  ggsidekick::theme_sleek()


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

sdr
tmb_pred <- ssdr[rownames(ssdr) %in% c("pred_eff"), "Estimate"]
tmb_pred_se <- ssdr[rownames(ssdr) %in% c("pred_eff"), "Std. Error"]


# slight deviation is due to RIs in TMB coded as random walk
# new_dat %>% 
#   mutate(tmb = tmb_pred,
#          gam = as.numeric(preds_m1$fit)
#          ) %>% 
#   ggplot(.) +
#   geom_point(aes(x = tmb, y = gam, col = region)) +
#   geom_abline(slope = 1, intercept = 0)


## look at predictions for areas and regions
# NOTE that these include predictions from random smooths (e.g. by = year) and
# random intercepts
pred_dat_tmb <- calc_ribbons(pred_dat = new_dat, fit_vec = tmb_pred,
                            se_vec = tmb_pred_se)

ggplot() +
  geom_line(data = pred_dat_tmb, aes(x = month_n, y = fit, colour = year)) +
  facet_wrap(~area, scales = "free_y") +
  ggsidekick::theme_sleek()


tmb_pred_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Estimate"]
tmb_pred_se_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Std. Error"]
new_dat_region <- new_dat %>% 
  select(month_n, region, year, offset, month, reg_month_year_f) %>% 
  distinct()
pred_dat_tmb_region <- calc_ribbons(pred_dat = new_dat_region, 
                                    fit_vec = tmb_pred_region,
                                    se_vec = tmb_pred_se_region)

ggplot() +
  geom_line(data = pred_dat_tmb_region, aes(x = month_n, y = fit, colour = year)) +
  geom_ribbon(data = pred_dat_tmb_region,
              aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year,
                  colour = year),
              alpha = 0.2) +
  facet_wrap(~region, scales = "free_y") +
  ggsidekick::theme_sleek()

