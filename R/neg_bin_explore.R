### Random Splines Neg Binomial
## Adapt stockseasonr model to account for varying splines among REs 
## Nov. 12, 2021

library(tidyverse)
library(mgcv)
library(TMB)

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


## stockseasonr function contents
yr_vec_catch <- as.numeric(as.factor(as.character(rec_catch$year))) - 1

# generate model matrix based on GAM with spline type a function of months
# retained
months1 <- unique(rec_catch$month_n)
spline_type <- if (max(months1) == 12) "cc" else "tp"
n_knots <- if (max(months1) == 12) 4 else 3


## fit model equivalent to ms
m1 <- mgcv::gam(catch ~ region + s(month_n, bs = spline_type, k = n_knots, 
                                 by = region) + 
                  s(year, bs="re") + offset,
                data = rec_catch,
                family = mgcv::nb()
)

preds_m1 <- predict.gam(m1, 
                     pred_dat, 
                     exclude = "s(year)",
                     se.fit = TRUE)

pred_dat_m1 <- pred_dat %>% 
  # filter(year == "2010") %>% 
  mutate(link_fit = as.numeric(preds_m1$fit),
         link_se_fit = as.numeric(preds_m1$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se_fit),
         link_up = link_fit + (qnorm(0.975) * link_se_fit),
         fit = exp(link_fit),
         fit_lo = exp(link_lo),
         fit_up = exp(link_up))

ggplot(pred_dat2) +
  geom_line(aes(x = month_n, y = fit, colour = year)) +
  # geom_ribbon(aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year, 
  #                 colour = year),
  #             alpha = 0.4) +
  facet_wrap(~region) +
  ggsidekick::theme_sleek()


## as above but with splines varying by year as well
m2 <- mgcv::gam(catch ~ region + s(month_n, bs = spline_type, k = 4, 
                                   by = region, m = 2) + 
                  s(month_n, by = year, bs="tp", m = 1, k = 4) +  
                  s(year, bs="re") + offset,
                data = rec_catch,
                family = mgcv::nb()
)

preds_m2 <- predict.gam(m2, 
                     pred_dat, 
                     exclude = "s(year)",
                     se.fit = TRUE)

pred_dat_m2 <- pred_dat %>% 
  mutate(link_fit = as.numeric(preds_m2$fit),
         link_se_fit = as.numeric(preds_m2$se.fit),
         link_lo = link_fit + (qnorm(0.025) * link_se_fit),
         link_up = link_fit + (qnorm(0.975) * link_se_fit),
         fit = exp(link_fit),
         fit_lo = exp(link_lo),
         fit_up = exp(link_up))

ggplot(pred_dat_m2) +
  geom_line(aes(x = month_n, y = fit, colour = year)) +
  # geom_ribbon(aes(x = month_n, ymin = fit_lo, ymax = fit_up, fill = year, 
  #                 colour = year),
  #             alpha = 0.4) +
  facet_wrap(~region) +
  ggsidekick::theme_sleek()


# compare model matrices
pred_mm_catch1 <- predict(m1, pred_dat, type = "lpmatrix")
pred_mm_catch2 <- predict(m2, pred_dat, type = "lpmatrix")

