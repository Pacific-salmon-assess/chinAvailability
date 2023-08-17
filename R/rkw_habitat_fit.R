## Model Fitting
# Use GSI samples through 2020 to explore composition coverage and determine
# spatial scale of models; equivalent process with creel data

library(tidyverse)
library(stockseasonr)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  # redefine region based on analysis
  mutate(
    cap_region = case_when(
      lat < 48.8 & lon > -125.25 & lon < -124.25 ~ "swiftsure",
      lat < 48.45 & lon < -123.4 & lon > -124.25 ~ "sooke",
      TRUE ~ "outside"
    ),
    whale_samples_time = ifelse(
      (year < 2011 | year > 2017) & month_n %in% c("6", "7", "8") & 
        rkw_habitat == "yes",
      "yes",
      "no"
    ),
    year = as.factor(year),
    yday = lubridate::yday(date),
    day_samp = paste(yday, year)
  )

comp_in <- rec_raw %>% 
  filter(
    !legal == "sublegal",
    !cap_region == "outside"
  ) %>% 
  mutate(
    sample_id = paste(cap_region, week_n,  #rkw_habitat,
                      year, sep = "_")
  ) %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)) %>% as.numeric,
    # add abbreviated key for model fitting
    stock_group2 = tolower(stock_group) %>% 
      str_replace(., " ", "_") %>% 
      paste("stock", ., sep = "-")
  ) %>% 
  ungroup() %>% 
  group_by(sample_id, rkw_habitat, cap_region, week_n, month_n, year, nn, 
           stock_group2) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  filter(
    month_n >= 5 & month_n <= 10
  )


# make predictive dataframe
area_key <- comp_in %>% 
  select(cap_region) %>% 
  distinct()

# subset predicted composition dataset
pred_dat_comp <- expand.grid(
  cap_region = unique(comp_in$cap_region),
  month_n = seq(min(comp_in$month_n),
                max(comp_in$month_n),
                by = 0.1)
  )



## stockseasonr fits -----------------------------------------------------------


fit <- fit_stockseasonr(
  comp_formula = stock_group2 ~ 1 + cap_region + 
    s(month_n, bs = "tp", k = 4) + (1 | year),
  comp_dat = comp_in,
  pred_dat = pred_dat_comp,
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  # nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)

fit$ssdr


## function to clean and pars and generate preds
clean_pred_foo <- function(fit, preds) {
  # make predictions
  ssdr <- fit$ssdr
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
  
  stock_seq <- colnames(fit$tmb_data$Y2_ik) %>% 
    str_remove(., "stock-")
  pred_out <- purrr::map(stock_seq, function (x) {
    dum <- preds
    dum$stock <- x
    return(dum)
  }) %>%
    bind_rows() %>%
    cbind(., link_preds) %>% 
    mutate(
      stock = factor(
        stock,
        levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer 4.1", 
                   "fraser_fall", "ecvi_somn", "psd", "other")
      )
    )
  
  # make mean observations
  obs_out <- fit$wide_comp_dat %>%
    mutate(samp_nn = apply(fit$tmb_data$Y2_ik, 1, sum)) %>%
    pivot_longer(cols = starts_with("stock-"), 
                 names_to = "stock", 
                 values_to = "obs_count",
                 names_prefix = "stock-") %>%
    mutate(obs_ppn = obs_count / samp_nn) %>% 
    filter(month_n %in% preds$month_n) %>% 
    mutate(
      stock = factor(
        stock,
        levels = c("nbc_seak", "wcvi", "fraser_yearling", "fraser_summer 4.1", 
                   "fraser_fall", "ecvi_somn", "psd", "other")
      )
    )
  
 list(
    preds = pred_out,
    obs_dat = obs_out
    )
}

pred_fit <- clean_pred_foo(fit = fit, preds = pred_dat_comp)



p <- ggplot(data = pred_fit$preds, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(cap_region~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) 

p_ribbon <- p +
  geom_ribbon(data = pred_fit$preds,
              aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)

p_obs <- p +
  geom_jitter(data = pred_fit$obs_dat,
              aes(x = month_n, y = obs_ppn, size = samp_nn), alpha = 0.2) +
  # geom_point(data = full_pred_list$mean_dat,
  #            aes(x = month_n, y = mean_obs_ppn), alpha = 0.6, color = "blue") +
  scale_size_continuous() +
  theme(legend.position = "top")
