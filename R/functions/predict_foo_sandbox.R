fit2 <- gam(
  agg_prob ~ 0 + stock_group + s(sg_year, bs = "re"),
  data = agg_dat, family = "tw", method = "REML"
)
fit2b <- gam(
  agg_prob ~ 0 + stock_group +  s(year, bs = "re"),
  data = agg_dat, family = "tw", method = "REML"
)


nd <- newdata %>% 
  filter(
    # year %in% c("2022", "2023"),
    week_n == "30",
    strata == "Renfrew"
  )
category_name = "stock_group"
origdata = agg_dat
newdata = nd

se_pred_ic = pred_ic = array(NA, dim=c(nrow(newdata),nlevels(origdata[,category_name])))
for(cI in 1:ncol(pred_ic)){
  
  # Modify data
  data = newdata
  data[,category_name] = factor( levels(origdata[,category_name])[cI], levels=levels(origdata[,category_name]) )
  data$sg_year <- paste(
    factor( levels(origdata[,category_name])[cI], levels=levels(origdata[,category_name]) ),
    data$year_n,
    sep = "_"
  ) %>% 
    as.factor()
  
  # Modify class
  class(x) = setdiff( class(x), "mvtweedie" )
  #class(x) = original_class
  
  # Apply predict.original_class
  if( "fit_model" %in% class(x) ){
    # if using VAST
    pred_ic[,cI] = predict(x,
                           what="D_i",
                           Lat_i=x$data_frame[,'Lat_i'],
                           Lon_i=x$data_frame[,'Lon_i'],
                           t_i=x$data_frame[,'t_i'],
                           a_i=x$data_frame[,'a_i'],
                           c_iz=rep(cI-1,nrow(x$data_frame)),
                           v_i=x$data_frame[,'v_i'] )
  }else{
    pred = predict(x,
                   newdata = data,
                   type="response",
                   se.fit = se.fit )
    if( se.fit==TRUE ){
      pred_ic[,cI] = pred$fit
      se_pred_ic[,cI] = pred$se.fit
    }else{
      pred_ic[,cI] = pred
    }
  }
}

# Normalize probability for each observation and class
rowsum_pred_ic = outer( rowSums(pred_ic), rep(1,ncol(pred_ic)) )
prob_ic = pred_ic / rowsum_pred_ic
prob_i = prob_ic[ cbind(1:nrow(pred_ic), match(newdata[,category_name],levels(origdata[,category_name]))) ]
