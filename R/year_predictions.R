### Explore Model Preds
## Take model predictions from stockseasonr and used to generate time series of 
## abundance in each region. Ultimately will need to generate predictions in TMB
## to integrate uncertainty in FE and RE.
## 1. Calculate aggregate abundance by year/region by summing area-specific 
## estimates within a year (based on gen_rand_int_pred from chinDist repo).
## 2. Calculate year/region-specific composition predictions.
## 3. Generate predictions for stock-specific catch using product.
## 4. Evaluate time series.
## Oct. 13, 2021


library(tidyverse)


# generated in chinDist directory
fit_tbl <- readRDS(here::here("data", "combined_model_dir.RDS"))
# pred_list <- readRDS(here::here("data", "combined_model_predictions.RDS"))


## AGGREGATE ABUNDANCE ---------------------------------------------------------

# fit with commercial data first
pred_dat_catch <- fit_tbl$pred_dat_catch[[1]]
comp_long <- fit_tbl$comp_long[[1]]
ssdr <- fit_tbl$ssdr[[1]]
catch_dat <- fit_tbl$catch_data[[1]]

# trim fit_tbl down
fit_tbl_trim <- fit_tbl %>% 
  filter(grouping == "pst")
  select(fishery, pred_dat_catch, comp_long, ssdr, catch_dat)

  
## Calculate area-specific values
pred_list <- pmap(
  list(fit_tbl_trim$pred_dat_catch,
       fit_tbl_trim$comp_long,
       fit_tbl_trim$ssdr,
       fit_tbl_trim$catch_data
       ),
  .f = function(pred_dat_catch, comp_long, ssdr, catch_dat) {
    abund_pred <- ssdr[rownames(ssdr) %in% "log_pred_abund", ] 
    pred_dat2 <- pred_dat_catch %>% 
      left_join(., comp_long %>% select(region, region_c) %>% distinct(), 
                by = "region") %>% 
      mutate(region_c = fct_reorder(region_c, as.numeric(region)),
             offset = unique(pred_dat_catch$offset))
    
    # fixed effects predictions
    area_fe_preds <- data.frame(est_link = abund_pred[ , "Estimate"],
               se_link =  abund_pred[ , "Std. Error"]) %>%
      cbind(pred_dat2, .) %>%
      mutate(
        pred_est = exp(est_link),
        pred_low = exp(est_link + (qnorm(0.025) * se_link)),
        pred_up = exp(est_link + (qnorm(0.975) * se_link)),
        pred_est_logcpue = est_link - offset,
        pred_low_logcpue = (est_link + (qnorm(0.025) * se_link)) - offset,
        pred_up_logcpue = (est_link + (qnorm(0.975) * se_link)) - offset
      )
    
    # save random intercept parameters and bind to vector of years based on input
    # data
    rand_int <- ssdr[rownames(ssdr) %in% "z1_k", ] 
    year_ints <- data.frame(year = levels(as.factor(as.character(catch_dat$year))),
                            z1_k = as.numeric(rand_int[, "Estimate"])) 
    
    # generate RE predictions
    year_preds <- expand.grid(year = year_ints$year, 
                              month = unique(pred_dat_catch$month)) %>% 
      left_join(pred_dat_catch %>% select(area, month), ., by = "month") %>% 
      left_join(., year_ints, by = "year")
    
    area_re_preds <- area_fe_preds %>% 
      select(month_n:est_link, se_link) %>% 
      left_join(., year_preds, by = c("area", "month")) %>% 
      mutate(
        year = as.factor(year),
        est_link_yr = est_link + z1_k,
        # pred_est_logcpue = est_link_yr / offset,
        pred_catch_yr = exp(est_link_yr),
        effort_yr = exp(offset)
      ) 
    
    # return both fixed and random effect predictions
    return(list(area_fe_preds = area_fe_preds, 
                area_re_preds = area_re_preds#,
                # reg_re_preds = reg_re_preds
                ))
  }
)

pred_list[[1]]$area_re_preds %>% 
  filter(month_n %in% seq(1, 12, by = 1)) %>% 
  select(area, region, offset, effort_yr) %>% 
  distinct()
  glimpse()

reg_pred_list <- map(pred_list, function (x) {
  x$area_re_preds %>% 
    group_by(month_n, region_c, year) %>% 
    summarize(pred_catch_yr = sum(pred_catch_yr), 
              eff_yr = sum(effort_yr),
              pred_catch_rate_yr = pred_catch_yr / eff_yr,
              .groups = "drop") %>% 
    mutate(year_n = as.numeric(as.character(year)),
           month_c = month.abb[month_n],
           month_f = fct_reorder(as.factor(month_c), month_n))
})


# remove predictions not associated with
reg_re_preds_comm <- reg_pred_list[[1]]  
reg_re_preds_rec <- reg_pred_list[[2]] 


## plotting functions
plot_years <- function(x, pal = "Dark2") {
  ggplot(x) +
    geom_line(aes(x = year_n, y = pred_catch_rate_yr, colour = region_c)) +
    scale_colour_brewer(palette = pal, name = "Region") +
    facet_wrap(~month_f) +
    ggsidekick::theme_sleek() +
    labs(x = "Year", y = "Predicted Catch Rate") 
}

plot_seasons <- function(x, pal = "D") {
  ggplot(x) +
    geom_line(aes(x = month_n, y = pred_catch_rate_yr, colour = year)) +
    scale_colour_viridis_d(option = pal, name = "Year") +
    facet_wrap(~region_c) +
    ggsidekick::theme_sleek() +
    labs(x = "Month", y = "Predicted Catch Rate") 
}


# time series of monthly averages in each region
pdf(here::here("figs", "yearly_preds.pdf"))
reg_re_preds_comm %>% 
  filter(month_n %in% seq(1, 12, by = 1)) %>% 
  plot_years()
reg_re_preds_rec %>% 
  filter(month_n %in% seq(1, 12, by = 1)) %>% 
  plot_years(pal = "Set1")

reg_re_preds_comm %>% 
  plot_seasons(pal = "A")
reg_re_preds_rec %>% 
  plot_seasons(pal = "D")

dev.off()

## COMPOSITION -----------------------------------------------------------------


# predicted effects in link space not generated by current model; calculate
# externally by multiplying composition beta parameters by relevant prediction
# matrix
tmb_data <- fit_tbl$tmb_data[[1]]
n_stocks <- ncol(tmb_data$y2_ig)
comp_pred_mat <- fit_tbl$tmb_data[[1]]$X2_pred_ij
comp_pred_dat <- fit_tbl$pred_dat_comp[[1]]
comp_wide_dat <- fit_tbl$comp_wide[[1]]
comp_dat <- fit_tbl$comp_long[[1]]
years <- levels(as.factor(as.character(comp_dat$year)))
stk_names <- unique(comp_dat$agg)

b2 <- ssdr[rownames(ssdr) %in% "b2_jg", ] 
b2_mat <- matrix(b2[, 1], nrow = 6, ncol = n_stocks)

# generate fixed effects predictions
fe_preds_link <- comp_pred_mat %*% b2_mat


# save random intercept parameters and bind to vector of years based on input
# data
rand_int_comp <- ssdr[rownames(ssdr) %in% "z2_k", ] 
year_ints_comp <- data.frame(year = years,
                             z2_k = as.numeric(rand_int_comp[, "Estimate"])) 


# Add random intercepts to fixed effects 
re_preds_link <- vector(mode = "list", length = length(years))
for (i in 1:length(years)) {
  re_preds_link[[i]] <- fe_preds_link + year_ints_comp$z2_k[i]
}


# Calculate derived quantities for each RE level
pred_pi_list <- pred_ppns_list <- vector(mode = "list", length = length(years))
for (i in 1:length(years)) {
  nn <- nrow(re_preds_link[[i]])
  pred_gamma = exp(re_preds_link[[i]])
  pred_gamma_plus = apply(pred_gamma, 1, sum)
  pred_theta = 1 / (pred_gamma_plus + 1)
  
  pred_pi_prop <- pred_pi <- matrix(NA, nrow = nrow(pred_gamma), 
                                    ncol = ncol(pred_gamma))
  for (j in 1:nn) {
    pred_pi[j, ] <- pred_gamma[j, ] / pred_theta[j]
  }
  
  pred_pi_plus <- apply(pred_pi, 1, sum)
  for (j in 1:nn) {
    pred_pi_prop[j, ] <- pred_pi[j, ] / pred_pi_plus[j]
  }
  
  colnames(pred_pi_prop) <- stk_names
  
  pred_pi_list[[i]] <- pred_pi
  pred_ppns_list[[i]] <- cbind(comp_pred_dat, as.data.frame(pred_pi_prop)) %>% 
    mutate(year = years[i])
}
  

pred_gamma = exp(pred_eff.array());
pred_gamma_plus = pred_gamma.rowwise().sum();
pred_theta = 1 / (pred_gamma_plus + 1);
for (int m = 0; m < n_pred_levels; m++) {
  for(int k = 0; k < n_cat; k++) {
    pred_pi(m, k) = pred_gamma(m, k) / pred_theta(m);
  }
}
pred_n_plus = pred_pi.rowwise().sum();
for (int m = 0; m < n_pred_levels; m++) {
  for (int k = 0; k < n_cat; k++) {
    pred_pi_prop(m, k) = pred_pi(m, k) / pred_n_plus(m);
    inv_logit_pred_pi_prop(m, k) = invlogit(pred_pi_prop(m, k));
  }
}