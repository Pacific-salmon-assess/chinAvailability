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

# utility functions for prepping smooths 
source(here::here("R", "utils.R"))

trim_rec_catch <- rec_catch %>% 
  filter(
    region %in% c("JdFS", "SSoG"),
    year %in% c("2015", "2016", "2017", "2018", "2019")
  ) %>% 
  droplevels() %>% 
  mutate(log_eff = offset)
trim_pred <- pred_dat %>% 
  filter(
    region %in% trim_rec_catch$region,
    year %in% "2015"#trim_rec_catch$year
  ) %>% 
  droplevels() %>% 
  mutate(log_eff = offset)



f1a <- catch ~ region +
  s(month_n, bs = "tp", k = 3) +
  log_eff
f1 <- catch ~ region +
  s(month_n, bs = "tp", k = 3) +
  log_eff
f2 <- catch ~ region +
  s(month_n, bs = "tp", k = 3, by = region, m = 2) +
  s(year, bs="re") #+
  # offset
f2 <- catch ~ region +
  s(month_n, bs = "tp", k = 3, by = region, m = 2) +
  s(month_n, by = year, bs="tp", m = 1, k = 3) +
  s(year, bs="re") + offset

m1 <- mgcv::gam(f1a,
                data = trim_rec_catch,
                family = mgcv::nb()
)

preds_m1 <- predict.gam(m1,
                        trim_pred,
                        # exclude = "s(year)",
                        se.fit = TRUE)

pred_dat_m1 <- trim_pred %>%
  # filter(year == "2010") %>%
  mutate(link_fit = as.numeric(preds_m1$fit),
         link_se_fit = as.numeric(preds_m1$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se_fit),
         link_up = link_fit + (qnorm(0.975) * link_se_fit),
         fit = exp(link_fit),
         fit_lo = exp(link_lo),
         fit_up = exp(link_up))

ggplot() +
  geom_point(data = trim_rec_catch, aes(x = month_n, y = catch, colour = year)) +
  geom_line(data = pred_dat_m1, aes(x = month_n, y = fit, colour = year)) +
  # geom_ribbon(aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year,
  #                 colour = year),
  #             alpha = 0.4) +
  facet_wrap(~region, scales = "free_y") +
  ggsidekick::theme_sleek()


## prep data for TMB

# generate model matrices
formula_no_sm <- remove_s_and_t2(f1)
X_ij <- model.matrix(formula_no_sm, data = dat)
sm <- parse_smoothers(f1, data = dat)
pred_X_ij <- predict(gam(formula_no_sm, data = dat), 
                     new_dat, type = "lpmatrix")
sm_pred <- parse_smoothers(f1, data = dat, newdata = new_dat)
# rand_int_vec <- as.numeric(as.factor(as.character(trim_rec_catch$year))) - 1

# offset
offset_pos <- grep("^offset$", colnames(X_ij))


# input data
data <- list(
  # abundance input data
  # y1_i = trim_rec_catch$catch,
  y1_i = dat$y,
  X1_ij = X_ij,
  # factor1k_i = rand_int_vec,
  # nk1 = length(unique(rand_int_vec)),
  b_smooth_start = sm$b_smooth_start,
  Zs = sm$Zs, # optional smoother basis function matrices
  Xs = sm$Xs, # optional smoother linear effect matrix
  pred_X1_ij = pred_X_ij,
  pred_Zs = sm_pred$Zs,
  pred_Xs = sm_pred$Xs
  # pred_factor2k_h = grouping_vec,
  # pred_factor2k_levels = grouping_key,
)

# input parameter initial values
pars <- list(
  b1_j = rep(0, ncol(X_ij)),
  ln_phi = log(1.5),
  bs = rep(0, ncol(sm$Xs)),
  ln_smooth_sigma = rep(0, length(sm$sm_dims)),
  b_smooth = rep(0, sum(sm$sm_dims))#,
  # z1_k = rep(0, length(unique(rand_int_vec))),
  # ln_sigma_zk1 = log(0.25)
)

# define offset
pars$b1_j[offset_pos] <- 1
b1_j_map <- seq_along(pars$b1_j)
b1_j_map[offset_pos] <- NA
tmb_map <- list(b1_j = as.factor(b1_j_map))

# define random parameters
tmb_random <- "b_smooth"

compile(here::here("src", "negbin_rsplines.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines")))

obj <- MakeADFun(data, pars, 
                 tmb_map = tmb_map,
                 random = tmb_random,
                 DLL = "negbin_rsplines")
opt <- nlminb(obj$par, obj$fn, obj$gr)
sdr <- sdreport(obj)
ssdr <- summary(sdr)

ssdr[rownames(ssdr) %in% c("b1_j", "log_phi", "bs"), ]
tmb_pred <- ssdr[rownames(ssdr) %in% c("pred_fe"), "Estimate"]

trim_pred %>% 
  mutate(tmb = tmb_pred,
         gam = as.numeric(preds_m1$fit)) %>% 
  # pivot_longer(cols = c(tmb, gam), names_to = "model", values_to = "fit") %>% 
  ggplot(.) +
  geom_point(aes(x = exp(tmb), y = exp(gam), col = region)) +
  geom_abline(slope = 1, intercept = 0)




# split_formula <- glmmTMB::splitForm(f1)
# RE_names <- barnames(split_formula$reTrmFormulas)
# fct_check <- vapply(RE_names, function(x) check_valid_factor_levels(data[[x]], .name = x), TRUE)
# RE_indexes <- vapply(RE_names, function(x) as.integer(data[[x]]) - 1L, rep(1L, nrow(data)))
# nobs_RE <- unname(apply(RE_indexes, 2L, max)) + 1L
# if (length(nobs_RE) == 0L) nobs_RE <- 0L
# formula <- split_formula$fixedFormula
# ln_tau_G_index <- unlist(lapply(seq_along(nobs_RE), function(i) rep(i, each = nobs_RE[i]))) - 1L


## sdmTMB test
set.seed(19203)
# examples from ?mgcv::gam.models
# continuous by example:
dat <- mgcv::gamSim(3, n = 800)
m_mgcv <- mgcv::gam(y ~ s(x2
                          # , by = x1
                          ), data = dat)
p_mgcv <- predict(m_mgcv)
dat$X <- runif(nrow(dat))
dat$Y <- runif(nrow(dat))
spde <- make_mesh(dat, c("X", "Y"), cutoff = 0.1)
m <- sdmTMB(y ~ s(x2
                  # , by = x1
                  ),
                           data = dat,
                           spde = spde,
            include_spatial = FALSE
)
p <- predict(m, newdata = NULL)
plot(p$est, p_mgcv)
abline(a = 0, b = 1)

m$tmb_data$Zs %>% glimpse()



## sim data
dat <- mgcv::gamSim(6, n = 800)
g <- exp(dat$f / 5)
dat$y <- rnbinom(g, size=3, mu = g)

## negative binomial data... 
b <- gam(y~s(x2, by = fac) + fac, 
         family=nb(link = "log"),
         data=dat)
f1 <- y~s(x2, by = fac) + fac

new_dat <- expand_grid(
  x2 = seq(min(dat$x2), max(dat$x2), length.out = 50),
  fac = unique(dat$fac)
)

preds <- predict.gam(b, newdata = new_dat)
new_dat$fit <- exp(preds)
ggplot() +
  geom_line(data = new_dat, aes(x = x2, y = fit, colour = fac)) +
  geom_point(data = dat, aes(x = x2, y = y, colour = fac))
