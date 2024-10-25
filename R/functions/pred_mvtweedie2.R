# modified version predict.mvtweedie that a) removes VAST architecture, 
# b) incorporates "exclude" argument from predict.mgcv, c) allows non-tweedie
# families (should be restricted to Poisson and Negative Binomial)

pred_dummy <- function( x,
            #                  original_class = "glmmTMB",
            category_name = "group",
            newdata,
            origdata = x$frame,
            se.fit = FALSE,
            exclude = NULL) {
    # Error checks
    # if( any(c("gam","glmmTMB") %in% class(x)) ){
    #   if( tolower(substr(family(x)$family,1,7)) != "tweedie" ) warning("`predict.mvtweedie` only implemented for a Tweedie distribution")
    #   if( family(x)$link != "log" ) stop("`predict.mvtweedie` only implemented for a log link")
    # }else if( "fit_model"%in%class(x) ){
    #   if( se.fit==TRUE ) error("se.fit not implemented for predict using VAST")
    # }else{
    #   stop("`predict.mvtweedie` only implemented for mgcv, glmmTMB and VAST")
    # }
    
    # Check and account for tibbles
    if( "package:tibble" %in% search() ){
      if( is_tibble(origdata) ){
        warning("Converting `origdata` from tibble to data.frame")
        origdata = as.data.frame(origdata)
      }
    }
    
    # Defaults
    if(missing(newdata) || is.null(newdata)) newdata = origdata
    
    # Predict each observation for each class
    se_pred_ic = pred_ic = array(
      NA, 
      dim = c(nrow(newdata), nlevels(origdata[,category_name]))
      )
    for(cI in 1:ncol(pred_ic)){
      
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
      
      # Modify class
      if ("mvtweedie" %in% class(x)) {
        class(x) = setdiff( class(x), "mvtweedie" )
      }
      # simulate to generate predictions from sdmTMB
      pred = predict(x,
                     newdata = data,
                     type="response",
                     se.fit = se.fit,
                     exclude = exclude)
      if( se.fit==TRUE ){
        pred_ic[,cI] = pred$fit
        se_pred_ic[,cI] = pred$se.fit
      }else{
        pred_ic[,cI] = pred
      }
    }
    
    
    # Normalize probability for each observation and class
    rowsum_pred_ic = outer( rowSums(pred_ic), rep(1,ncol(pred_ic)) )
    # rowsum_pred2_ic = outer( rowSums(pred2_ic), rep(1,ncol(pred2_ic)) )
    prob_ic = pred_ic / rowsum_pred_ic
    # prob2_ic = pred2_ic / rowsum_pred2_ic
    prob_i = prob_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ]
    
    # return prediction
    if( se.fit==TRUE ){
      # Normalize SE-squared for each observation and class
      rowsum_se2_ic = outer( rowSums(se_pred_ic^2), rep(1,ncol(pred_ic)) )
      se2_prob_ic = prob_ic^2 * ( se_pred_ic^2/pred_ic^2 - 2*se_pred_ic^2/(pred_ic*rowsum_pred_ic) + rowsum_se2_ic/rowsum_pred_ic^2 )
      se_i = sqrt(se2_prob_ic[ cbind(1:nrow(se2_prob_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ])
      out = list("fit"=prob_i, "se.fit"=se_i)
    }else{
      out = prob_i
    }

    return(out)
  }
