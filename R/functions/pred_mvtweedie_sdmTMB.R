# modified version predict.mvtweedie that a) removes VAST architecture, 
# b) incorporates "exclude" argument from predict.mgcv, c) allows non-tweedie
# families (should be restricted to Poisson and Negative Binomial)

# x = fit_sdmTMB
# se.fit = TRUE
# category_name = "stock_group"
# origdata = agg_dat
# newdata = newdata
# nsim = 500

pred_mvtweedie_sdmTMB <- function(
    x,
    category_name = "group",
    newdata,
    origdata = x$frame,
    se.fit = FALSE,
    re_form_iid = NULL,
    nsim = 250) {
    
    # Defaults
    if(missing(newdata) || is.null(newdata)) newdata = origdata
    
    # Predict each observation for each class
    if (se.fit == TRUE) {
      pred_ic = array(
        NA, 
        dim = c(nrow(newdata), nsim, nlevels(origdata[,category_name]))
      )
    } else if (se.fit == FALSE) {
      pred_ic = array(
        NA, 
        dim = c(nrow(newdata), nlevels(origdata[,category_name]))
      )
    }
    for(cI in 1:nlevels(origdata[,category_name])) {
      
      # Modify data
      data = newdata
      data[,category_name] = factor(
        levels(origdata[,category_name])[cI], 
        levels=levels(origdata[,category_name]) 
        )
      
      # include hacky version to account for year-group specific RIs
      if ("sg_year" %in% colnames(data)) {
        data$sg_year <- paste(
          factor(levels(origdata[,category_name])[cI], 
                 levels=levels(origdata[,category_name])),
          data$year_n,
          sep = "_"
        ) %>% 
          as.factor()
      }
      
      # simulate to generate predictions from sdmTMB
      if( se.fit == TRUE ){
        pred = predict(x,
                       newdata = data,
                       type="response",
                       nsim = nsim,
                       re_form_iid = re_form_iid)
        pred_ic[ , , cI] = pred
      } else{
        pred = predict(x,
                       newdata = data,
                       type="response",
                       re_form_iid = re_form_iid)
        pred_ic[,cI] = pred$est
      }
    }
    
    
    # Normalize probability for each observation and class
    if (se.fit == FALSE) {
      rowsum_pred_ic = outer( rowSums(pred_ic), rep(1,ncol(pred_ic)) )
      prob_ic = pred_ic / rowsum_pred_ic
      prob_i = prob_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ]
      out = prob_i
    }
    
    if( se.fit==TRUE ){
      #create array to store normalized proportions with same dims as
      # pred_ic
      prp_array <- array(NA, dim = dim(pred_ic))
      for (i in 1:nsim) {
        temp_row_sums <- outer(
          rowSums(pred_ic[ , i, ]), rep(1, ncol(pred_ic[ , i, ]))
          )
        prp_array[ , i, ] <- pred_ic[ , i, ] / temp_row_sums
      }
      prob_ic <- se_ic <- low_ic <- up_ic <- array(
        NA, 
        dim = c(nrow(newdata), nlevels(origdata[,category_name]))
      )
      # calculate summary stats among simulation values
      for (j in 1:nlevels(origdata[,category_name])) {
        prob_ic[ , j] <- apply(prp_array[ , , j], 1, mean)
        se_ic[ , j] <- apply(prp_array[ , , j], 1, sd)
        low_ic[ , j] <- apply(prp_array[ , , j], 1, quantile, 0.025)
        up_ic[ , j] <- apply(prp_array[ , , j], 1, quantile, 0.975)
      }
      # save as list and rearrange into vectors with length nrow(newdata)
      out = list(
        "fit" = prob_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ],
        "se.fit" = se_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ], 
        "low" = low_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ], 
        "up" = up_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ]
        )
    }

    return(out)
  }



