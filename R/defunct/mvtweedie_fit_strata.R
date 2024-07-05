## Model Fitting
# Fit data to estimate seasonal changes in stock composition 
# Adaptation of rkw_strata_fit to fit a spatially explicit model using mvtweedie
# and mgcv or glmmtmb packages


library(tidyverse)
library(mvtweedie)
library(mgcv)
library(glmmTMB)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() 


dat <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    !legal == "sublegal",
    !strata == "other"
  ) %>% 
  mutate(
    sample_id = paste(strata, week_n, year, sep = "_"),
  ) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "km",
    utm_names = c("utm_x", "utm_y")
  ) %>% 
  droplevels() %>% 
  group_by(sample_id) %>% 
  mutate(
    sample_id_n = sum(prob),
    # calculate average location within a sampling event
    utm_y = mean(utm_y),
    utm_x = mean(utm_x),
    strata = factor(
        strata,
        levels = c("swiftsure", "swiftsure_nearshore", "renfrew", "vic",
                   "haro", "saanich"),
        labels = c("Swiftsure", "Nitinat", "Renfrew", "Sooke\n/Victoria",
                   "S. Gulf\nIslands", "Saanich")
      )
  ) %>% 
  ungroup()

sample_key <- dat %>% 
  select(sample_id, sample_id_n, strata, year, week_n, utm_y, utm_x) %>% 
  distinct()


# add zero observations
agg_dat <- expand.grid(
  sample_id = unique(dat$sample_id),
  stock_group = unique(dat$stock_group)
) %>% 
  left_join(., sample_key, by = "sample_id") %>% 
  left_join(
    ., 
    dat %>% 
      group_by(sample_id, stock_group) %>% 
      summarize(
        agg_prob = sum(prob)
      ) %>% 
      ungroup(),
    by = c("sample_id", "stock_group")
  ) %>% 
  mutate(
    agg_prob = ifelse(is.na(agg_prob), 0, agg_prob),
    agg_ppn = agg_prob / sample_id_n,
    year = as.factor(year),
    stock_group = as.factor(stock_group),
    utm_x_m = utm_x * 1000,
    utm_y_m = utm_y * 1000
  ) 
saveRDS(
  agg_dat, here::here("data", "rec", "cleaned_ppn_data_rec.rds")
)


## FIT MODEL -------------------------------------------------------------------

system.time(
  fit <- gam(
    agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 5, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds") +
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds") ,
    data = agg_dat, family = "tw",
    knots = list(week_n = c(0, 52))
  )
)
system.time(
  fit2 <- gam(
    agg_prob ~ 0 + stock_group + s(week_n, by = stock_group, k = 5, bs = "cc") +
      s(utm_y, utm_x, m = c(0.5, 1), bs = "ds") +
      s(utm_y, utm_x, by = stock_group, m = c(0.5, 1), bs = "ds") +
      s(year, bs = "re"),
    data = agg_dat, family = "tw",
    knots = list(week_n = c(0, 52))
  )
)

data( "Middleton_Island_TUPU", package="mvtweedie" )
Middleton_Island_TUPU$Year = as.numeric(as.character( Middleton_Island_TUPU$Year_factor ))

# Fit GAM
fit = mgcv::gam( Response ~ group + s(Year, by=group), data=Middleton_Island_TUPU, family="tw" )


## CHECKS ----------------------------------------------------------------------

## simulate based on:
# https://gist.github.com/dantonnoriega/ad2081c39b26d0f523ba3464f4a90282

fit_raw <- fit2

# use fitted GAM for mean estimates
phi.hat <- fit_raw$deviance/sum(fit_raw$prior.weights)
mu.hat <- fitted(fit_raw)
p.hat <- fit_raw$family$getTheta(TRUE)
prob_zero <- exp(-mu.hat^(2-p.hat) / phi.hat / (2-p.hat))

# single sim
y.tw <- mgcv::rTweedie(mu.hat, p = p.hat, phi = phi.hat)

# plot generated y vs simulated y from fitted values
y <- agg_dat$agg_prob
# y <- Middleton_Island_TUPU$Response
brks = seq(0, ceiling(max(max(y), max(y.tw))), by = 0.5)
hist(y, breaks = brks, col = scales::alpha('red', .9))
par(new = TRUE)
hist(y.tw, breaks = brks, col = scales::alpha('blue', .5), axes = FALSE, 
     xlab = NULL, ylab = NULL, main = NULL)


# full simulate
nsims <- 50
sim_mat <- matrix(NA, nrow = length(mu.hat), ncol = nsims)
for (i in 1:ncol(sim_mat)) {
  sim_mat[ , i] <- mgcv::rTweedie(mu.hat, p = p.hat, phi = phi.hat)
}

mu_pred <- predict(fit_raw, newdata = agg_dat) %>%
  fit_raw$family$linkinv(.)
# mu_pred <- predict(fit_raw, newdata = Middleton_Island_TUPU) %>% 
#   fit_raw$family$linkinv(.)

dharma_res <- DHARMa::createDHARMa(
  simulatedResponse = sim_mat,
  observedResponse = y,
  fittedPredictedResponse = mu_pred
)
plot(dharma_res)


# ppn zeros
sum(agg_dat$agg_prob == 0) / nrow(agg_dat)
sum(sim_mat < 0.001) / length(sim_mat)
sum(y == 0) / length(y)


agg_dat$sim_dat <- sim_mat[ , 2]
agg_dat$resid <- agg_dat$sim_dat - agg_dat$agg_prob

ggplot(agg_dat) +
  geom_point(aes(x = utm_x, y = utm_y, colour = resid)) +
  facet_wrap(
    ~ stock_group
  ) +
  scale_colour_gradient2() +
  ggsidekick::theme_sleek()
