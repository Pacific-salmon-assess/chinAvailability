### Random Splines Neg Binomial
## Adapt stockseasonr model to account for varying splines among REs 
## Feb 18, 2022
# Updated to focus on area-specific estimates of abundance

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)

# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))


dat <- readRDS(here::here("data", "rec", "rec_creel_area.rds")) %>% 
  filter(!legal == "sublegal",
         reg %in% c("JdFS"#, "SSoG"
         ),
         month_n > 4 & month_n < 10) %>% 
  mutate(month = as.factor(round(month_n, 0)),
         year = as.factor(year),
         area = as.factor(area),
         reg = as.factor(reg),
         offset = log(effort),)

mean_eff <- dat %>% 
  group_by(area, month) %>% 
  summarize(offset = mean(offset, na.rm = T),
            .groups = "drop")

## predicted dataset
new_dat <- group_split(dat, reg) %>%
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
            dat %>% select(reg, area) %>% distinct(),
            by = "area"
  ) %>%
  mutate(strata = paste(month_n, reg, sep = "_"),
         offset = mean(dat$offset)
         ) %>%
  arrange(reg) %>%
  # convoluted to ensure ordering is correct for key
  mutate(
    month = as.factor(round(month_n, 0)),
    order = row_number(),
    # year is necessary regardless of whether RIs are inlcuded because of year-
    # specific smooths 
    reg_month_year = paste(reg, month_n, year, sep = "_"),
    reg_month_year_f = fct_reorder(factor(reg_month_year), order)
  ) %>%
  # left_join(., mean_eff, by = c("month", "area")) %>%
  select(-order, -reg_month_year, -strata) %>%
  distinct() 


obs_cpue <- dat %>% 
  group_by(reg, month_n, year) %>% 
  summarize(sum_catch = sum(catch),
            sum_eff = sum(effort),
            ln_cpue = log(sum_catch / sum_eff),
            .groups = "drop")


## MGCV VS TMB -----------------------------------------------------------------


f1a <- catch ~ area +
  s(month_n, bs = "tp", k = 3) +
  s(month_n, bs = "tp", k = 3, by = area, m =  1) +
  s(month_n, by = year, bs = "tp", m = 1, k = 3) +
  #smooth REs can't be readily passed with Anderson's function; code RIs instead
  s(year, bs = "re") + 
  offset(offset)
f1 <- catch ~ area +
  s(month_n, bs = "tp", k = 3) +
  s(month_n, bs = "tp", k = 3, by = area, m = 1) +
  s(month_n, by = year, bs = "tp", k = 3, m = 1) +
  offset


# fit mgcv for comparison
# m1 <- mgcv::gam(f1,
#                 data = rec_catch,
#                 family = mgcv::nb()
# )
m1a <- mgcv::gam(f1a,
                data = rec_catch,
                family = mgcv::nb()
)


# first check against observed
preds_m1a <- predict.gam(m1,
                         exclude = "s(year)",
                         se.fit = TRUE)
dat$pred_mean <- preds_m1a$fit %>% as.numeric()
dat$pred_se <- preds_m1a$se.fit %>% as.numeric()
dat$ln_cpue <- log(dat$catch / dat$effort)
  
ggplot() +
  geom_abline() +
  facet_wrap(~area) +
  ggsidekick::theme_sleek() +
  geom_point(data = dat, aes(x = log(catch), y = pred_mean, 
                             colour = year))
ggplot(dat) +
  geom_point(aes(x = month_n, y = log(catch)), colour = "black") +
  geom_point(aes(x = month_n, y = pred_mean), colour = "red") +
  facet_grid(year~area) +
  ggsidekick::theme_sleek()

## MGCV has good performance

preds_m1 <- predict.gam(m1a,
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
# extreme deviations due to standardizing effort
ggplot() +
  geom_line(data = pred_dat_m1, aes(x = month_n, y = link_fit - offset, 
                                    colour = year)) +
  geom_point(data = obs_cpue, aes(x = month_n, y = ln_cpue)) +
  facet_grid(year~area, scales = "free_y") +
  ggsidekick::theme_sleek()


##
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

compile(here::here("src", "negbin_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines")))

obj <- MakeADFun(data, pars, 
                 map = tmb_map,
                 random = tmb_random,
                 DLL = "negbin_rsplines")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

saveRDS(ssdr, here::here("data", "model_fits", "negbin_rsmooths.rds"))


tmb_pred <- ssdr[rownames(ssdr) %in% c("pred_mu1"), "Estimate"]
tmb_pred_se <- ssdr[rownames(ssdr) %in% c("pred_mu1"), "Std. Error"]


# slight deviation is due to RIs in TMB coded as random walk
new_dat %>%
  mutate(tmb = tmb_pred,
         gam = as.numeric(preds_m1$fit)
         ) %>%
  ggplot(.) +
  geom_point(aes(x = tmb, y = gam, col = area)) +
  geom_abline(slope = 1, intercept = 0)


## look at predictions for areas and regions
# NOTE that these include predictions from random smooths (e.g. by = year) and
# random intercepts
pred_dat_tmb <- calc_ribbons(pred_dat = new_dat, fit_vec = tmb_pred,
                            se_vec = tmb_pred_se)

ggplot() +
  geom_line(data = pred_dat_tmb, aes(x = month_n, y = fit, colour = year)) +
  facet_wrap(~area, scales = "free_y") +
  ggsidekick::theme_sleek()


# tmb_pred_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Estimate"]
# tmb_pred_se_region <- ssdr[rownames(ssdr) %in% c("ln_pred_eff_cumsum"), "Std. Error"]
# new_dat_region <- new_dat %>% 
#   select(month_n, region, year, offset, month, reg_month_year_f) %>% 
#   distinct()
# pred_dat_tmb_region <- calc_ribbons(pred_dat = new_dat_region, 
#                                     fit_vec = tmb_pred_region,
#                                     se_vec = tmb_pred_se_region)
# 
# ggplot() +
#   geom_line(data = pred_dat_tmb_region, aes(x = month_n, y = fit, colour = year)) +
#   geom_ribbon(data = pred_dat_tmb_region,
#               aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year,
#                   colour = year),
#               alpha = 0.2) +
#   facet_wrap(~region) +
#   ggsidekick::theme_sleek()
# 
# 
# 
# 
# ggplot() +
#   geom_line(data = pred_dat_tmb_region, 
#             aes(x = month_n, y = link_fit - unique(pred_dat$offset), 
#                 colour = year)) +
#   geom_point(data = obs_cpue,
#               aes(x = month_n, y = ln_cpue, colour = year),
#               alpha = 0.4) +
#   facet_wrap(~region) +
#   ggsidekick::theme_sleek()
# 
# 
# 
# dat %>% 
#   mutate(cpue = catch / eff) %>% 
#   ggplot(., aes(x = month_n)) +
#   labs(y = "CPUE", x = "Month") +
#   facet_grid(area~year, scales = "free_y") +
#   ggsidekick::theme_sleek() +
#   geom_bar(aes(y = cpue, fill = region), stat = "identity") 
