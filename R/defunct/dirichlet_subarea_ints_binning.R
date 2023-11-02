### Stock Composition Estimates
## Extension of dirichlet_subarea_ints.R to evaluate effects of binning samples
## at four scales
## Aug 9, 2022

library(tidyverse)
library(TMB)


# tmb models - use MVN if time-varying predictions are required, use RI if 
compile(here::here("src", "dirichlet_ri_sdmTMB.cpp"))
dyn.load(dynlib(here::here("src", "dirichlet_ri_sdmTMB")))


# utility functions for prepping smooths 
# source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit_new.R"))

comp_list <- readRDS( here::here("data", "rec", "comp_dat_list.rds"))
comp_tbl <- tibble(
  name = names(comp_list),
  dat = comp_list
)

# stock_comp_area <- stock_comp
# stock_comp_area_w <- stock_comp
# stock_comp_reg <- stock_comp
# stock_comp_reg_w <- stock_comp
# 
# comp_list <- list(stock_comp_area, stock_comp_area_w, stock_comp_reg, 
#                   stock_comp_reg_w)
# names(comp_list) <- c("subarea_yday", "subarea_week",
#                       "reg_yday", "reg_week")
# saveRDS(comp_list, here::here("data", "rec", "comp_dat_list.rds"))


# add areas to composition dataset
area_key <- readRDS(here::here("data", "rec", "subarea_key.RDS"))


# subset predicted composition dataset
pred_dat_comp <- #pred_dat_comp1 %>% 
  expand.grid(
    reg = unique(comp_tbl$dat[[1]]$reg),
    month_n = seq(min(comp_tbl$dat[[1]]$month_n),
                  max(comp_tbl$dat[[1]]$month_n),
                  by = 0.1
    )
  ) %>% 
  left_join(., area_key, by = "reg") %>%
  filter(#core_area == "yes",
         month_n < 9.1 & month_n > 4.9
         ) %>% 
  select(-core_area) %>% 
  mutate(week = month_n * 4) 



## FIT MODEL -------------------------------------------------------------------

fit_list <- comp_tbl$dat %>% 
  purrr::map(
  function (x) {
    fit_stockseasonr(
      comp_formula = agg_new ~ 1 + subarea + 
        s(month_n, bs = "tp", k = 4, m = 2) + 
        (1 | year),
      comp_dat = x,
      pred_dat = pred_dat_comp,
      model = "dirichlet",
      fit = TRUE,
      nlminb_loops = 2, newton_loops = 1
    )
  }
    )
comp_tbl$tmb_data <- purrr::map(fit_list, function (x) x$tmb_data)
comp_tbl$ssdr <- purrr::map(fit_list, function (x) x$ssdr)



## EVALUATE MODEL PREDS --------------------------------------------------------

pred_foo <- function (ssdr, tmb_data) {
  logit_pred_ppn <- ssdr[rownames(ssdr) == "logit_pred_Pi_prop", ]
  pred_mu <- ssdr[rownames(ssdr) == "pred_Mu2", ]
  
  link_preds <- data.frame(
    link_prob_est = logit_pred_ppn[ , "Estimate"],
    link_prob_se =  logit_pred_ppn[ , "Std. Error"],
    pred_mu = pred_mu[ , "Estimate"]
  ) %>% 
    mutate(
      pred_prob_est = plogis(link_prob_est),
      pred_prob_low = plogis(link_prob_est + (qnorm(0.025) * link_prob_se)),
      pred_prob_up = plogis(link_prob_est + (qnorm(0.975) * link_prob_se))
    ) 
  
  stock_seq <- colnames(tmb_data$Y2_ik)
  pred_comp <- purrr::map(stock_seq, function (x) {
    dum <- pred_dat_comp
    dum$stock <- x
    return(dum)
  }) %>%
    bind_rows() %>%
    cbind(., link_preds) %>% 
    mutate(
      area_f = fct_reorder(as.factor(area_f), as.numeric(reg)),
      stock = fct_relevel(
        stock, "Fraser_Yearling", "Fraser_S", "Fraser_F",
        "SOG", "PSD", "WCVI", "WA_OR_CA"
      )
    ) %>% 
    select(-subarea_original) %>% 
    distinct()
}

pred_comp <- purrr::pmap(list(comp_tbl$ssdr, comp_tbl$tmb_data, comp_tbl$name), 
            function(x, y, name) {
              pred_foo(x, y) %>% 
                mutate(dataset = name)
            }) %>% 
  bind_rows()
  

p <- ggplot(data = pred_comp, aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(subarea~stock) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est, color = dataset)) #+

ggplot(data = pred_comp %>% filter(stock == "Fraser_S", subarea == "29DE"), 
       aes(x = month_n)) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  # facet_grid(subarea~dataset) +
  facet_wrap(~dataset)+
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) +
  geom_ribbon(aes(ymin = pred_prob_low, ymax = pred_prob_up),
              alpha = 0.2)



