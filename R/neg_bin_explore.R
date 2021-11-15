### Random Splines Neg Binomial
## Adapt stockseasonr model to account for varying splines among REs 
## Nov. 12, 2021

library(tidyverse)
library(mgcv)

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
  droplevels() %>% 
  # group by region for now
  group_by(region, region_c, month, month_n, year) %>%
  summarize(catch = sum(catch),
            eff = sum(eff),
            offset = log(eff)) %>%
  ungroup()

## predicted dataset
pred_dat <- group_split(rec_catch, region) %>%
  map_dfr(., function(x) {
    expand.grid(
      region = unique(x$region),
      # area = unique(x$area),
      month_n = seq(min(x$month_n),
                    max(x$month_n),
                    by = 0.1
      ),
      year = unique(x$year)
    )
  }) %>% 
  # left_join(.,
  #           rec_catch %>% select(region, area) %>% distinct(),
  #           by = "area"
  # ) %>%
  mutate(strata = paste(month_n, region, sep = "_"),
         offset = mean(rec_catch$offset)) %>%
  # filter(strata %in% comp_strata) %>%
  arrange(region) %>%
  # convoluted to ensure ordering is correct for key
  mutate(
    month = as.factor(month_n),
    order = row_number(),
    reg_month = paste(region, month_n, sep = "_"),
    reg_month_f = fct_reorder(factor(reg_month), order)
  ) %>%
  select(-order, -reg_month, -strata) %>%
  distinct() 


## MGCV version ----------------------------------------------------------------

## stockseasonr function contents
yr_vec_catch <- as.numeric(as.factor(as.character(rec_catch$year))) - 1

# generate model matrix based on GAM with spline type a function of months
# retained
months1 <- unique(rec_catch$month_n)
spline_type <- if (max(months1) == 12) "cc" else "tp"
n_knots <- if (max(months1) == 12) 4 else 3


## fit model equivalent to ms
m1 <- mgcv::gam(catch ~ region + s(month_n, bs = spline_type, k = n_knots, 
                                 by = region, m = 2) + 
                  s(year, bs="re") + offset,
                data = rec_catch,
                family = mgcv::nb()
)

preds_m1 <- predict.gam(m1, 
                     pred_dat, 
                     # exclude = "s(year)",
                     se.fit = TRUE)

pred_dat_m1 <- pred_dat %>% 
  # filter(year == "2010") %>%
  mutate(link_fit = as.numeric(preds_m1$fit),
         link_se_fit = as.numeric(preds_m1$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se_fit),
         link_up = link_fit + (qnorm(0.975) * link_se_fit),
         fit = exp(link_fit),
         fit_lo = exp(link_lo),
         fit_up = exp(link_up),
         model = "r_ints")

ggplot(pred_dat_m1) +
  geom_line(aes(x = month_n, y = fit, colour = year)) +
  # geom_ribbon(aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year,
  #                 colour = year),
  #             alpha = 0.4) +
  facet_wrap(~region, scales = "free_y") +
  ggsidekick::theme_sleek()


## as above but with splines varying by year as well
m2 <- mgcv::gam(catch ~ region + s(month_n, bs = spline_type, k = 3, 
                                   by = region, m = 2) + 
                  s(month_n, by = year, bs="tp", m = 1, k = 3) +  
                  s(year, bs="re") + offset,
                data = rec_catch,
                family = mgcv::nb()
)

preds_m2 <- predict.gam(m2, 
                     pred_dat, 
                     exclude = "s(year)",
                     se.fit = TRUE)

pred_dat_m2 <- pred_dat %>% 
  # filter(year == "2010") %>% 
  mutate(link_fit = as.numeric(preds_m2$fit),
         link_se_fit = as.numeric(preds_m2$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se_fit),
         link_up = link_fit + (qnorm(0.975) * link_se_fit),
         fit = exp(link_fit),
         fit_lo = exp(link_lo),
         fit_up = exp(link_up),
         model = "r_slopes")

ggplot(pred_dat_m2) +
  geom_line(aes(x = month_n, y = fit, colour = year)) +
  # geom_ribbon(aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year, 
  #                 colour = year),
  #             alpha = 0.4) +
  facet_wrap(~region, scales = "free_y") +
  ggsidekick::theme_sleek()


## heat plot of seasonal patterns
rbind(pred_dat_m1, pred_dat_m2) %>% 
  ggplot(., 
       aes(x = month_n, y = year)) +
  geom_raster(aes(fill = fit)) +
  scale_fill_viridis_c(name = "Predicted\nAbundance") +
  labs(x = "Month", y = "Year") +
  ggsidekick::theme_sleek() +
  facet_grid(model~region)


# as above but with relative abundance (not very informative...)
rbind(pred_dat_m1, pred_dat_m2) %>% 
  group_by(region, model) %>% 
  mutate(mean_abund = mean(fit),
         sd_abund = sd(fit),
         scaled_abund = (fit - mean_abund) / sd_abund) %>% 
  ungroup() %>% 
  ggplot(., 
         aes(x = month_n, y = year)) +
  geom_raster(aes(fill = scaled_abund)) +
  scale_fill_viridis_c(name = "Predicted\nScaled\nAbundance") +
  labs(x = "Month", y = "Year") +
  ggsidekick::theme_sleek() +
  facet_grid(model~region)


# compare model matrices
pred_mm_catch1 <- predict(m1, pred_dat, type = "lpmatrix")
pred_mm_catch2 <- predict(m2, pred_dat, type = "lpmatrix")


## TMB VERSION -----------------------------------------------------------------


library(TMB)
library(sdmTMB)


trim_rec_catch <- rec_catch %>% 
  filter(
    region %in% c("JdFS", "SSoG"),
    year %in% c("2015", "2016", "2017", "2018", "2019")
  )
trim_pred <- pred_dat %>% 
  filter(
    region %in% trim_rec_catch$region,
    year %in% trim_rec_catch$year
  )


m1 <- mgcv::gam(catch ~ region + 
                  s(month_n, bs = "tp", k = 3, by = region, m = 2) + 
                  offset,
                data = trim_rec_catch,
                family = mgcv::nb()
)


n<-400
dat <- gamSim(1, n=n)
# g <- exp(dat$f/5)

## negative binomial data... 
# dat$y <- rnbinom(g,size=3,mu=g)
# ## known theta fit ...
# b0 <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=negbin(3),data=dat)
# plot(b0,pages=1)
# print(b0)
# 
# ## same with theta estimation...
# b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3),family=nb(link = "log"),data=dat)
# plot(b,pages=1)
# print(b)
# b$family$getTheta(TRUE) ## extract final theta estimate

b <- gam(y~s(x0)+s(x1)+s(x2)+s(x3), data=dat)



f1 <- catch ~ region + 
  s(month_n, bs = "tp", k = 3, by = region, m = 2) + 
  # s(year, bs="re") + 
  offset
f2 <- catch ~ region + 
  s(month_n, bs = "tp", k = 3, by = region, m = 2) + 
  s(year, bs="re") +
  offset
f2 <- catch ~ region + 
  s(month_n, bs = "tp", k = 3, by = region, m = 2) + 
  s(month_n, by = year, bs="tp", m = 1, k = 3) +  
  s(year, bs="re") + offset


# utility functions for prepping smooths 
source(here::here("R", "utils.R"))



## prep data for TMB

# generate model matrices
formula_no_sm <- remove_s_and_t2(f1)
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(f1, data = dat)
# sm2 <- parse_smoothers(f2, data = dat)
# sm_pred <- parse_smoothers(f1, data = dat, newdata = trim_pred)
# yr_vec_catch <- as.numeric(as.factor(as.character(dat$year))) - 1

# input data
data <- list(
  # abundance input data
  # y1_i = trim_rec_catch$catch,
  y1_i = dat$y,
  X1_ij = X_ij,
  # factor1k_i = yr_vec_catch,
  # nk1 = length(unique(yr_vec_catch)),
  b_smooth_start = sm$b_smooth_start,
  Zs = sm$Zs, # optional smoother basis function matrices
  Xs = sm$Xs # optional smoother linear effect matrix
  # X1_pred_ij = pred_mm_catch,
  # pred_factor2k_h = grouping_vec,
  # pred_factor2k_levels = grouping_key,
)

# input parameter initial values
pars <- list(
  b1_j = rnorm(ncol(X_ij), 0, 0.01),
  phi = 1.5,
  # log_phi = log(1.5),
  bs = rep(0, ncol(sm$Xs)),
  ln_smooth_sigma = rep(0, length(sm$sm_dims)),
  b_smooth = rep(0, sum(sm$sm_dims))#,
  # z1_k = rep(0, length(unique(yr_vec_catch))),
  # log_sigma_zk1 = log(0.25)
)


compile(here::here("src", "negbin_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines")))

obj <- MakeADFun(data, pars, DLL = "negbin_rsplines")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

ssdr[rownames(ssdr) %in% c("b1_j", "log_phi", "bs"), ]
