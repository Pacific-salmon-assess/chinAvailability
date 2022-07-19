## Function to make inputs for composition or integrated model, then fit
# Eventually move to stockseasonr

# Dec. 7, 2021

# Dev TO DO:
# 2) Make RE intercepts flexible (i.e. turn on/off or remove MVN component)
# 3) Map 0 observations in composition component
# 4) Combine cpp scripts into one with conditionals (may not be worthwhile...)

library(tidyverse)
library(mgcv)
library(TMB)
library(sdmTMB)

 
# abund_formula = catch ~ 1 +
#   s(month_n, bs = "tp", k = 4, m = 2) +
#   (1 | reg) +
#   (1 | year)
# # abund_dat = catch;
# # abund_rint = "year";
# # comp_formula = pst_agg ~ area + s(month_n, bs = "tp", k = 4, m = 2);
# # comp_dat = stock_comp
# # comp_rint = NULL #"year"
# # pred_dat = pred_dat_catch
# model = "negbin"
# include_re_preds = FALSE

## MAKE INPUTS  ----------------------------------------------------------------

# helper function to map unused parameters
map_foo <- function(x, tmb_pars) {
  out_list <- vector(mode = "list", length = length(x))
  names(out_list) <- x
  for (i in seq_along(x)) {
    out_list[[i]] <- as.factor(rep(NA, length(tmb_pars[[x[i]]])))
  }
  return(out_list)
}

make_inputs <- function(abund_formula = NULL, comp_formula = NULL, 
                        comp_knots = NULL,
                        abund_dat = NULL, comp_dat = NULL,
                        abund_offset = NULL,
                        # abund_rint = NULL, comp_rint = NULL,
                        pred_dat = NULL,
                        model = c("negbin", "dirichlet", "integrated"),
                        random_walk = FALSE,
                        include_re_preds = FALSE) {
  
  # make sure necessary components are present
  if (model != "negbin") {
    if (is.null(comp_dat)) {
      stop("Missing model inputs to fit integrated model")
    }
    if (is.null(comp_dat$prob)) {
      stop("Composition data not identified. Name vector in comp_dat 'prob' to 
       indicate proportions data.")
    }
  }
  
  # initialize empty lists to fill with data and initial parameters
  tmb_data <- list()
  tmb_pars <- list()
  tmb_map <- list()
  tmb_random <- NULL

  if (model %in% c("integrated", "negbin")) {
    # generate inputs for negbin component
    # random intercepts (use sdmTMB to generate proper structure)
    sdmTMB_dummy <- sdmTMB::sdmTMB(
      abund_formula,
      data = abund_dat,
      spatial = "off",
      do_fit = FALSE
    )
    
    # generate offset separately
    if (is.null(offset)) {
      offset <- rep(0, length(y_i))
    } else {
      offset <- abund_dat[["offset"]]
    }
    
    # smooths present (conditional for determining input structure)
    has_smooths <- as.integer(sdmTMB_dummy$tmb_data$has_smooths)
    
    # predictions present
    has_preds <- ifelse(is.null(pred_dat), 0, 1)
    
    if (has_preds == 1) {
      resp <- attr(terms(abund_formula), which = "variables")[[2]] %>% 
        as.character()
      pred_dat[resp] <- 0
      
      sdmTMB_dummy_p <- sdmTMB::sdmTMB(
        abund_formula,
        data = pred_dat,
        spatial = "off",
        do_fit = FALSE
      )
    } else {
      # copy original template inputs as placeholders (not used in cpp due to
      # conditional)
      sdmTMB_dummy_p <- sdmTMB_dummy
    }
    
    # make abundance tmb inputs
    abund_tmb_data <- list(
      y1_i = sdmTMB_dummy$tmb_data$y_i %>% as.numeric(),
      X1_ij = sdmTMB_dummy$tmb_data$X_ij[[1]],
      re_index1 = sdmTMB_dummy$tmb_data$RE_indexes,
      ln_sigma_re_index1 = sdmTMB_dummy$tmb_data$ln_tau_G_index,
      nobs_re1 = sdmTMB_dummy$tmb_data$nobs_RE, # number of random intercepts
      Zs = sdmTMB_dummy$tmb_data$Zs, # optional smoother basis function matrices
      Xs = sdmTMB_dummy$tmb_data$Xs, # optional smoother linear effect matrix
      offset_i = offset,
      has_smooths = has_smooths,
      b_smooth_start = sdmTMB_dummy$tmb_data$b_smooth_start,
      has_preds = has_preds,
      pred_X1_ij = sdmTMB_dummy_p$tmb_data$X_ij[[1]],
      pred_Zs = sdmTMB_dummy_p$tmb_data$Zs,
      pred_Xs = sdmTMB_dummy_p$tmb_data$Xs
    )
    tmb_data <- c(tmb_data, abund_tmb_data)
    
    # TODO: adjust data when RI predictions being made
    # if (include_re_preds == TRUE) {
    #   #vector of predicted random intercepts
    #   pred_rand_ints <- list(
    #     pred_rfac1 = as.numeric(pred_dat[[abund_rint]]) - 1)
    # 
    #   tmb_data <- c(tmb_data, pred_rand_ints)
    # }
    
    abund_tmb_pars <- list(
      b1_j = rep(0, ncol(abund_tmb_data$X1_ij)),
      ln_phi = log(1.5),
      bs = if (has_smooths) {
        sdmTMB_dummy$tmb_params$bs %>% as.vector() 
      } else {
        0
      },
      ln_smooth_sigma = if (has_smooths) {
        sdmTMB_dummy$tmb_params$ln_smooth_sigma %>% as.vector()
      } else {
        0
      },
      b_smooth = if (has_smooths) {
        rep(0, nrow(sdmTMB_dummy$tmb_params$b_smooth))
      } else {
        0
      },
      #for some reason cannot fix as zeros when ignore_fix = FALSE; use sdmTMB
      #inputs for length but redefine with rnorm
      re1 = rnorm(n = length(sdmTMB_dummy$tmb_params$RE %>% as.vector()),
                  0,
                  0.5),
      ln_sigma_re1 = sdmTMB_dummy$tmb_params$ln_tau_G %>% as.vector()
    )
    tmb_pars <- c(tmb_pars, abund_tmb_pars)
    
    # map parameters unless necessary
    # random parameters
    if (abund_tmb_data$nobs_re1[[1]] > 0) {
      tmb_random <- c(tmb_random, "re1")
    } else {
      re_map <- map_foo(x = c("re1", "ln_sigma_re1"),
                        tmb_pars = tmb_pars)
      tmb_map <- c(tmb_map, 
                   re_map)
    }
    if (abund_tmb_data$has_smooths > 0) {
      tmb_random <- c(tmb_random, "b_smooth")
    } else {
      smooth_map <- map_foo(x = c("b_smooth", "ln_smooth_sigma", "bs"),
                            tmb_pars = tmb_pars)
      tmb_map <- c(tmb_map, 
                   smooth_map)
    }
  } 
  
  if (model %in% c("integrated", "dirichlet")) {
    ## composition component of model
    # adjust composition model formula
    comp_formula_split <- str_split(comp_formula, "~")
    comp_formula_new <- formula(paste("dummy", 
                                      comp_formula_split[[3]], sep = "~"))
    group_var <- comp_formula_split[[2]]
    
    # adjust input data to wide and convert observations to matrix
    comp_wide <- comp_dat %>%
      pivot_wider(names_from = as.name(group_var), 
                  values_from = prob) %>%
      mutate_if(is.numeric, ~ replace_na(., 0.00001)) %>% 
      #add dummy response variable
      mutate(dummy = 0)
    group_names <- unique(comp_dat[[group_var]])
    obs_comp <- comp_wide[ , group_names] %>%
      as.matrix()
    
    # dummy model
    dummy_comp <- mgcv::gam(comp_formula_new, data = comp_wide, 
                            knots = comp_knots)
    X2_ij <- predict(dummy_comp, type = "lpmatrix")
    pred_X2_ij <- predict(dummy_comp, pred_dat, type = "lpmatrix")
    
    # check to make sure predictive dataframes for composition and abundance
    # are same length
    if (model == "integrated"){
      if (nrow(pred_X2_ij) != nrow(pred_X_ij)) {
        stop("Dimensions of abundance and composition predictions are not
         compatible.")
      }}
    
    ## TODO: conditionals don't seem to work (perhaps due to improper format
    # of empty data/parameters passed to tmb)
    # if (!is.null(comp_rint)) {
    rfac2 <- as.numeric(as.factor(as.character(comp_wide[[comp_rint]]))) - 1
    n_rint2 <- length(unique(rfac2))
    # } else {
    #   rfac2 <- matrix(, nrow = nrow(comp_wide), ncol = 0) %>% as.vector()
    #   n_rint2 <- 0
    # }
    
    #make composition tmb inputs 
    comp_tmb_data <- list(
      Y2_ik = obs_comp,
      X2_ij = X2_ij,
      rfac2 = rfac2,
      n_rfac2 = n_rint2,
      pred_X2_ij = pred_X2_ij,
      random_walk = ifelse(random_walk == TRUE, 1, 0) %>% as.integer()
    )
    tmb_data <- c(tmb_data, comp_tmb_data)
    
    comp_tmb_pars <- list(
      B2_jk = matrix(0,
                     nrow = ncol(X2_ij),
                     ncol = ncol(obs_comp)
      ),
      ln_sigma_A2 = log(0.25)
    )
    tmb_pars <- c(tmb_pars, comp_tmb_pars)
    
    # adjust data and parameters when RI predictions being made
    if (include_re_preds == TRUE) {
      #vector of predicted random intercepts
      #only added for dirichlet because generated in neg bin component for 
      #integrated model
      if (model == "dirichlet") {
        pred_rand_ints <- list(
          pred_rfac1 = as.numeric(pred_dat[[comp_rint]]) - 1
        )
        tmb_data <- c(tmb_data, pred_rand_ints)
      }
      # mvn matrix of REs
      mvn_rand_inits <- list(
        A2_hk = matrix(rnorm(n_rint2 * ncol(obs_comp), 0, 0.5), 
                       nrow = n_rint2,
                       ncol = ncol(obs_comp))
      )
      
      tmb_pars <- c(tmb_pars, mvn_rand_inits)
      tmb_random <- c(tmb_random, "A2_hk")
    } else if (include_re_preds == FALSE) {
      # vector of random intercepts
      rand_inits <- list(
        A2_h = rep(0, n_rint2)
      )
      
      tmb_pars <- c(tmb_pars, rand_inits)
      tmb_random <- c(tmb_random, "A2_h")
    }
  }
  
   
  # combine model specs to pass as logicals to fitting function
  model_specs <- list(model = model, include_re_preds = include_re_preds)
  
  out_list <- list(model_specs = model_specs,
                   tmb_data = tmb_data, tmb_pars = tmb_pars, tmb_map = tmb_map, 
                   tmb_random = tmb_random)
  if (model %in% c("dirichlet", "integrated")) {
    out_list <- c(out_list, list(wide_comp_dat = comp_wide))
  }
  
  return(out_list)
}


## FIT MODELS  -----------------------------------------------------------------
# tmb_data = tmb_inputs$tmb_data 
# tmb_pars = tmb_inputs$tmb_pars 
# tmb_map = tmb_inputs$tmb_map
# tmb_random  = tmb_inputs$tmb_random
# fit_random = TRUE;
# ignore_fix = FALSE;
# model_specs = tmb_inputs$model_specs
# tmb_pars$re1 <- rnorm(n = length(tmb_pars$re1), 0, 1)
# model_specs = list(model = "dirichlet",
#                    include_re_preds = FALSE)
# nlminb_loops = 2
# newton_loops = 2


fit_model <- function(tmb_data, tmb_pars, tmb_map = NULL, tmb_random = NULL,
                      # fit_random = TRUE, ignore_fix = FALSE, 
                      model_specs,
                      nlminb_loops = 1L, newton_loops = 0L
                      ) {
  
  if (model_specs$model == "negbin") tmb_model <- "negbin_rsplines_sdmTMB"
  # use MVN model if random effects predictions necessary
  if (model_specs$model == "dirichlet") {
    tmb_model <- ifelse(model_specs$include_re_preds == FALSE,
                        "dirichlet_ri",
                        "dirichlet_mvn")
  }
  if (model_specs$model == "integrated") {
    tmb_model <- ifelse(model_specs$include_re_preds == FALSE,
                        "negbin_rsplines_dirichlet_ri",
                        "negbin_rsplines_dirichlet_mvn")
  }
  
  ## fit fixed effects only 
  # map random effects
  # if (fit_random == FALSE | ignore_fix == FALSE) {
    #TODO: this is nonfunctional with smooths (b_smooth flagged as RI); 
    #for now remove entirely, but could replace with conditional excluding those 
    # pars although may not be worthwhile
  #   new_map_list <- tmb_pars[names(tmb_pars) %in% tmb_random]
  #   tmb_map_random <- c(
  #     tmb_map,
  #     map(new_map_list, function (x) factor(rep(NA, length(x))))
  #   )
  #   # fit
  #   obj <- TMB::MakeADFun(
  #     data = tmb_data,
  #     parameters = tmb_pars,
  #     map = tmb_map_random,
  #     DLL = tmb_model
  #   )
  #   opt <- stats::nlminb(obj$par, obj$fn, obj$gr,
  #                         control = list(eval.max = 1e4, iter.max = 1e4)
  #   )
  # }
  
  ## fit with random effects 
  # if (fit_random) {
    # pass parameter inits from above unless specified otherwise
    # if (ignore_fix == TRUE) {
      # pars_in <- tmb_pars   
    # } else {
    #   pars_in <- obj$env$parList(opt$par)
    # } 
    
    obj <- TMB::MakeADFun(
      data = tmb_data,
      parameters = tmb_pars, #pars_in,
      map = tmb_map,
      random = tmb_random,
      DLL = tmb_model
    )
    opt <- stats::nlminb(obj$par, obj$fn, obj$gr)
  # }
  
  if (nlminb_loops > 1) {
    for (i in seq(2, nlminb_loops, length = max(0, nlminb_loops - 1))) {
      temp <- opt[c("iterations", "evaluations")]
      opt <- stats::nlminb(
        start = opt$par, objective = obj$fn, gradient = obj$gr)
      opt[["iterations"]] <- opt[["iterations"]] + temp[["iterations"]]
      opt[["evaluations"]] <- opt[["evaluations"]] + temp[["evaluations"]]
    }
  }
  
  if (newton_loops > 0) {
    for (i in seq_len(newton_loops)) {
      g <- as.numeric(obj$gr(opt$par))
      h <- stats::optimHess(opt$par, fn = obj$fn, gr = obj$gr)
      opt$par <- opt$par - solve(h, g)
      opt$objective <- obj$fn(opt$par)
    }
  }
  
  sdr <- sdreport(obj)
  
  # derived quantities
  ssdr <- summary(sdr)
  
  return(list(sdr = sdr, ssdr = ssdr))
}
