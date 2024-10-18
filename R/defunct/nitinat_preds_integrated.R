### Nitinat Stock Specific Abundance Estimates
## Generate estimates of stock-specific abundance for area near
## Nitinat Lake (21-0 and 121-1). 
## Analogue to nitinat_preds_composition.R which has different seasonal 
## and spatial coverage
## To evaluate data quality and model performance generate:
## 1) Exploratory plots of catch effort
## 2) Observed observed CPUE and observed CPUE * composition 
##    (comp in other script)
## 3) Model predictions for restricted (21-0 and 121-1), intermediate (all 
## 21 and 121) and full (SWVI and JdF) datasets. 
## Dec. 8, 2022

library(tidyverse)
# library(TMB)
library(stockseasonr)


## DATA CLEAN ------------------------------------------------------------------


## clean catch data first since more constrained 
catch_subarea_raw <- readRDS(here::here("data", "rec", "rec_creel_subarea.rds"))
catch_subarea <- catch_subarea_raw %>% 
  mutate(
    subarea_original = subarea,
    subarea = case_when(
      subarea %in% c("21A", "121C") ~ "21A-121C",
      subarea %in% c("20DO", "20D", "20C") ~ "20D-20C",
      subarea %in% c("20A", "20E") ~ "20A-20E",
      TRUE ~ subarea
    ),
    # define core areas 
    zone = case_when(
      subarea %in% c("121A", "21A-121C") ~ "core",
      subarea %in% c("121B", "20A-20E", "20E", "123R"
                     # , "Area 20 (WCVI)", "Area 20 (West)"
      ) ~ "intermediate",
      subarea %in% c("20B", "20D-20C", "23J"
                     # "Area 20 (East)",
      ) ~ "full",
      TRUE ~ "out"
    ),
    zone = factor(zone, levels = c("core", "intermediate", "full", "out")),
    # correct apparently missing effort data
    effort = ifelse(catch == "0" & is.na(effort),
                    0,
                    effort)
  ) %>%
  # remove sublegals then sum catch pooling clip rates and AD marks
  filter(
    !zone == "out",
    legal == "legal"#,
    # subarea_est == "yes"
  ) %>% 
  # first sum by original subarea to combine all catches (effort the same)
  group_by(
    month_n, year, subarea, subarea_original, zone, effort
  ) %>% 
  summarize(
    catch = sum(catch),
    .groups = "drop"
  ) %>% 
  # then sum by new subarea to combine catch AND effort across originally distinct
  # subareas
  group_by(
    month_n, year, subarea, zone
  ) %>% 
  summarize(
    catch = sum(catch),
    effort = sum(effort),
    cpue = catch / effort,
    offset = log(effort),
    year = as.factor(year),
    .groups = "drop"
  ) %>% 
  distinct() %>% 
  # remove 0 CPUE observations (skeptical)
  filter(
    !is.na(cpue)
  )


## clean composition and filter to match catch
comp_in_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) 
comp_in <- comp_in_raw %>% 
  rename(stock_region = region) %>% 
  mutate(
    # define stock groups
    agg_new = case_when(
      grepl("NITINAT", stock) ~ "Nitinat",
      grepl("SWVI", resolved_stock_rollup) ~ "SWVI",
      grepl("NWVI", resolved_stock_rollup) ~ "other",
      pst_agg == "WCVI" ~ "SWVI",
      grepl("CR", pst_agg) ~ "other",
      grepl("CST", pst_agg) ~ "other",
      pst_agg == "NBC_SEAK" ~ "other",
      # grepl("Fraser", Region1Name) | pst_agg == "FR-early" ~ "Fraser",
      Region1Name %in% c("Fraser_Summer_4.1") ~ "Fraser_S",
      Region1Name %in% c("Fraser_Summer_5.2", "Fraser_Spring_5.2",
                         "Fraser_Spring_4.2") | pst_agg == "FR-early" ~
        "Fraser_Yearling",
      Region1Name %in% c("Fraser_Fall") ~ "Fraser_F",
      TRUE ~ pst_agg
    ),
    subarea_original = subarea,
    subarea = case_when(
      subarea %in% c("21A", "121C") ~ "21A-121C",
      subarea %in% c("20DO", "20D", "20C") ~ "20D-20C",
      subarea %in% c("20A", "20E") ~ "20A-20E",
      TRUE ~ subarea
    ),
    # subarea = ifelse(subarea %in% c("20DO", "20D"), "20D", subarea),
    # define core areas 
    zone = case_when(
      subarea %in% c("121A", "21A-121C") ~ "core",
      subarea %in% c("121B", "20A-20E", #"20E", "20A", #"123U",
                     "123R") ~ "intermediate",
      subarea %in% c("20B", "20D-20C", #"20E",
                     "23J") ~ "full",
      TRUE ~ "out"
    ),
    zone = factor(zone, levels = c("core", "intermediate", "full", "out"))
  ) 
  

# filter a bit
trim_stock <- comp_in %>%
    filter(
      !legal == "sublegal",
      !zone == "out",
      subarea %in% catch_subarea$subarea
    ) %>% 
    mutate(
      month = lubridate::month(date, label = TRUE),
      month_n = lubridate::month(date),
      week = lubridate::week(date),
      yday = lubridate::yday(date),
      sample_id = paste(month_n, subarea, #week,
                        year, sep = "_"),
      area_f = as.factor(as.numeric(gsub("([0-9]+).*$", "\\1", area)))
    ) %>% 
    group_by(sample_id) %>% 
    mutate(nn = length(unique(id)) %>% as.numeric) %>% 
    ungroup() %>% 
  group_by(sample_id, subarea, subarea_original, area_f, reg = cap_region, 
           #week, 
           month, month_n, year, nn, agg_new, zone) %>% 
  summarize(prob = sum(prob), .groups = "drop") %>% 
  ungroup() %>% 
  filter(
    month_n < 10.1 & month_n > 1.9,
    !reg == "out" 
  ) %>% 
  mutate(
    year = as.factor(year),
    subarea = as.factor(subarea),
    agg_new = paste("stock", agg_new, sep = "-")
  ) %>% 
  droplevels() 


# make predictive dataframe
area_key <- trim_stock %>% 
  select(subarea, #area_f, reg, 
         zone) %>% 
  distinct()

# subset predicted composition dataset
pred_dat_comp <- expand.grid(
  subarea = unique(trim_stock$subarea),
  month_n = seq(min(trim_stock$month_n),
                max(trim_stock$month_n),
                by = 0.1
  )
  ) %>% 
  left_join(., area_key, by = "subarea") %>%
  filter(
    zone %in% c("core", "intermediate"),
    month_n < 9.1 & month_n > 5.9
  )



## FIT MODEL -------------------------------------------------------------------

int_catch <- catch_subarea %>% 
  filter(
    zone %in% c("core", "intermediate"),
    month_n %in% c("6", "7", "8", "9")
  )
int_stock <- trim_stock %>% 
  filter(
    zone %in% c("core", "intermediate"),
    month_n %in% c("6", "7", "8", "9")
  )

fit_intermediate <- fit_stockseasonr(
  abund_formula = catch ~ 1 +
    s(month_n, bs = "tp", k = 3, m = 2) +
    subarea +
    (1 | year),
  abund_dat = int_catch,
  abund_offset = "offset",
  comp_formula = agg_new ~ 1 + subarea +  
    s(month_n, bs = "tp", k = 4, m = 2)
  + (1 | year),
  comp_dat = int_stock,
  pred_dat = pred_dat_comp,
  model = "integrated",
  random_walk = TRUE,
  fit = TRUE,
  nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)


full_catch <- catch_subarea %>% 
  filter(
    month_n %in% c(#"5", 
      "6", "7", "8", "9")
  )
full_stock <- trim_stock %>% 
  filter(
    month_n %in% c(#"5",
      "6", "7", "8", "9")
  )

fit_full <- fit_stockseasonr(
  abund_formula = catch ~ 1 +
    s(month_n, bs = "tp", k = 3, m = 2) +
    subarea +
    (1 | year),
  abund_dat = full_catch,
  abund_offset = "offset",
  comp_formula = agg_new ~ 1 + subarea + 
    s(month_n, bs = "tp", k = 4, m = 2) + (1 | year),
  comp_dat = full_stock,
  pred_dat = pred_dat_comp,
  model = "integrated",
  random_walk = TRUE,
  fit = TRUE,
  nlminb_loops = 2,
  newton_loops = 1,
  silent = FALSE
)


## MAKE PREDICTIONS ------------------------------------------------------------


clean_int_pred_foo <- function(fit, preds) {
  # make predictions
  ssdr <- fit$ssdr
  # log_pred_mu1 <- ssdr[rownames(ssdr) == "log_pred_mu1", ]
  log_pred_abund <- ssdr[rownames(ssdr) == "log_pred_mu1_Pi", ]
  # logit_pnn <- ssdr[rownames(ssdr) == "logit_pred_Pi_prop", ]

  link_preds <- data.frame(
    link_abund_est = log_pred_abund[ , "Estimate"],
    link_abund_se =  log_pred_abund[ , "Std. Error"]#,
    # link_ppn_est = logit_pnn[ , "Estimate"],
    # link_ppn_se =  logit_pnn[ , "Std. Error"]
  ) %>% 
    mutate(
      pred_ppn_mu_est = exp(link_abund_est),
      pred_ppn_mu_low = exp(link_abund_est + (qnorm(0.025) * link_abund_se)),
      pred_ppn_mu_up = exp(link_abund_est + (qnorm(0.975) * link_abund_se))
    ) 
  
  stock_seq <- colnames(fit$tmb_data$Y2_ik)
  dd <- purrr::map(stock_seq, function (x) {
    dum <- preds
    dum$stock <- x
    return(dum)
  }) %>%
    bind_rows() %>%
    cbind(., link_preds) 
  
}

int_preds <- clean_int_pred_foo(fit_intermediate,
                                 preds = pred_dat_comp)
full_preds <- clean_int_pred_foo(fit_full,
                                 preds = pred_dat_comp)


comb_preds <- purrr::map2(
  list(
    # core_preds,
    int_preds,
    full_preds
  ),
  c(#"core",
    "int", "full"),
  ~ {
    .x %>% 
      mutate(
        model = .y
      )
  }
) %>% 
  bind_rows() %>% 
  mutate(
    stock = str_split(stock, "-") %>% 
      purrr::map(., ~ .x[2]) %>% 
      as.character() %>% 
      c()
  ) %>% 
  mutate(
    subarea = fct_reorder(subarea, as.numeric(zone)),
    stock = factor(
      stock, 
      levels = c("Nitinat", "SWVI", 
                 # "Fraser",
                 "Fraser_Yearling", "Fraser_S", "Fraser_F",
                 "PSD", "SOG", "other")
    )
  ) %>% 
  distinct()

# export for use in RMD
saveRDS(comb_preds,
        here::here(
          "data", "model_fits", "nitinat_only", "integrated_preds.RDS"
        ))


p <- ggplot(data = comb_preds, aes(x = month_n)) +
  labs(y = "Predicted Index of Abundance", x = "Month") +
  facet_grid(subarea~stock, scales = "free_y") +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_ppn_mu_est, color = model))

p_ribbon <- p +
  geom_ribbon(data = comb_preds,
              aes(ymin = pred_ppn_mu_low, ymax = pred_ppn_mu_up, fill = model),
              alpha = 0.2)


png(here::here("figs", "nitinat_preds", "pred_ribbons_abund.png"),
    height = 5, width = 6.5, units = "in", res = 250)
p_ribbon
dev.off()


# single model versions of above
pdf(here::here("figs", "nitinat_preds", "abund_preds_fullmodel.pdf"),
    height = 6, width = 9)
ggplot(data = comb_preds %>% filter(model == "full"), 
       aes(x = month_n)) +
  labs(y = "Predicted Index of Abundance", x = "Month") +
  facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_ppn_mu_est)) +
  geom_ribbon(aes(ymin = pred_ppn_mu_low, ymax = pred_ppn_mu_up),
              alpha = 0.2)
dev.off()


## EXPLORATORY FIGS ------------------------------------------------------------


## Look at available catch data
catch_sampling <- catch_subarea %>%
  group_by(year, subarea, month_n, zone) %>%
  tally() %>%
  ggplot() +
  geom_jitter(aes(x = month_n, y = year, size = n, colour = zone),
              alpha = 0.5, width = 0.25) +
  facet_wrap(~fct_reorder(subarea, as.numeric(as.factor(zone)))) +
  ggsidekick::theme_sleek()

seasonal_cpue <- ggplot(catch_subarea) +
  geom_boxplot(aes(x = as.factor(month_n), y = cpue, fill = zone)) +
  facet_wrap(~fct_reorder(subarea, as.numeric(as.factor(zone)))) +
  ggsidekick::theme_sleek()

seasonal_effort <- ggplot(catch_subarea %>% 
                            filter(subarea %in% c("121A", "21A-121C"))) +
  geom_boxplot(aes(x = as.factor(month_n), y = effort)) +
  facet_wrap(~fct_reorder(subarea, as.numeric(as.factor(zone)))) +
  ggsidekick::theme_sleek() +
  labs(
    x = "Month",
    y = "Summed Vessel Days Per Month"
  )

pdf(here::here("figs", "nitinat_preds", "catch_effort_raw.pdf"))
catch_sampling
seasonal_cpue
seasonal_effort
dev.off()


