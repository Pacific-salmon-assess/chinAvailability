## Utility Functions

# take stockseasonr predictions and clean 

clean_pred_foo <- function(fit, preds) {
  ssdr <- fit$ssdr
  
  # stock names
  stock_seq <- colnames(fit$tmb_data$Y2_ik) %>% 
    str_remove(., "stock-")
  
  # clean parameter estimates
  par_names <- fit$tmb_data$X2_ij %>% colnames()
  par_name_vec <- list(length = length(stock_seq), mode = "vector")
  for (i in seq_along(stock_seq)) {
    par_name_vec[[i]] <- paste(par_names, stock_seq[i], sep = "_")
  }
  
  # par_est <- ssdr[rownames(ssdr) == "B2_jk", ]
  # rownames(par_est) <- unlist(par_name_vec)
  par_dat <- data.frame(
    Parameter = unlist(par_name_vec)
  ) %>% 
    cbind(as.data.frame(ssdr[rownames(ssdr) == "B2_jk", ])) %>% 
    janitor::clean_names()
  
  # clean predicted proportions
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
  
  pred_out <- purrr::map(stock_seq, function (x) {
    dum <- preds
    dum$stock <- x
    return(dum)
  }) %>%
    bind_rows() %>%
    cbind(., link_preds) 
  
  # make mean observations
  obs_out <- fit$wide_comp_dat %>%
    mutate(samp_nn = apply(fit$tmb_data$Y2_ik, 1, sum)) %>%
    pivot_longer(cols = starts_with("stock-"), 
                 names_to = "stock", 
                 values_to = "obs_count",
                 names_prefix = "stock-") %>%
    mutate(obs_ppn = obs_count / samp_nn) %>% 
    filter(week_n %in% preds$week_n) 
  
  list(
    pars = par_dat,
    preds = pred_out,
    obs_dat = obs_out
  )
}


clean_pred_foo_size <- function(fit, preds) {
  
  ssdr <- fit$ssdr
  
  # stock names
  stock_seq <- colnames(fit$tmb_data$Y2_ik) %>% 
    str_remove(., "sizebin-")
  
  # clean parameter estimates
  par_names <- fit$tmb_data$X2_ij %>% colnames()
  par_name_vec <- list(length = length(stock_seq), mode = "vector")
  for (i in seq_along(stock_seq)) {
    par_name_vec[[i]] <- paste(par_names, stock_seq[i], sep = "_")
  }
  
  par_dat <- data.frame(
    Parameter = unlist(par_name_vec)
  ) %>% 
    cbind(as.data.frame(ssdr[rownames(ssdr) == "B2_jk", ])) %>% 
    janitor::clean_names()
  
  # clean predicted proportions
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
  
  pred_out <- purrr::map(stock_seq, function (x) {
    dum <- preds
    dum$size_bin2 <- x
    return(dum)
  }) %>%
    bind_rows() %>%
    cbind(., link_preds) 
  
  # make mean observations
  obs_out <- fit$wide_comp_dat %>%
    mutate(samp_nn = apply(fit$tmb_data$Y2_ik, 1, sum)) %>%
    pivot_longer(cols = starts_with("sizebin-"), 
                 names_to = "size_bin2", 
                 values_to = "obs_count",
                 names_prefix = "sizebin-") %>%
    mutate(obs_ppn = obs_count / samp_nn) %>% 
    filter(week_n %in% preds$week_n) 
  
  list(
    pars = par_dat,
    preds = pred_out,
    obs_dat = obs_out
  )
}
