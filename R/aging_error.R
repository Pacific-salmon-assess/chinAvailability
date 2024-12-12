# Calculate Age Classification Errror
# Uses data exported from clean_data.R
# Nov 22, 2024
# 1) Use marine fishery and freshwater escapement samples as classification
# accuracy data
# 2) Fit multinomial model then generate stock specific aging accuracy estimates
# 3) Export posterior to combine with size-at-age predictions to generate 
# dataset of join posterior that can be sampled to derive size distribution 
# probabilities and ultimately size class assignment vector

library(tidyverse)
library(brms)

rstan::rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

age_dat <- readRDS(here::here("data", "rec", "aging_data.rds"))

# make age_key 
age_key <- data.frame(
  age_gr = c("21", "31", "32", "33", "41", "42", "51", "52", "53", "62"),
  est_age = c(2, 3, 3, 3, 4, 4, 5, 5, 5, 6)
)


# import confusion matrices from Fraser and convert to long DF
fr_stocks <- c("FR_Sum_4.1", "FR_Fall", "FR_Sum_5.2", "FR_Spr_4.2")
stock_list <- vector(mode = "list", length = length(fr_stocks))
for (i in seq_along(fr_stocks)) {
  stock_list[[i]] <- readxl::read_xlsx(
    here::here("data", "rec", "aging_error.xlsx"),
    sheet = i
  ) %>%
    pivot_longer(
      cols = starts_with("cwt_"),
      names_to = "cwt_age",
      names_prefix = "cwt_",
      values_to = "count"
    ) %>% 
    mutate(
      stock_group = fr_stocks[i],
      age_source = "cwt",
      location = "freshwater",
      cwt_age = as.numeric(cwt_age)
    )
}
fr_stock_dat <- bind_rows(stock_list) %>% 
  rename(est_age = scale_age, age = cwt_age) %>% 
  #expand to individual observations
  uncount(weights = count) 


dum <- age_dat %>% 
  left_join(., age_key, by = "age_gr") %>%
  mutate(age = age,
         location = "marine") %>% 
  select(est_age, age, stock_group, age_source, location) %>%
  rbind(., fr_stock_dat) %>% 
  mutate(
    bias_n = age - est_age,
    bias = case_when(
      age - est_age == 0 ~ "zero",
      age - est_age > 0 ~ "under",
      age - est_age < 0 ~ "over"
    ) %>% 
      factor(., levels = c("zero", "under", "over")),
    est_age = as.factor(est_age),
    age = as.factor(age)
  ) %>% 
  mutate(
    stock_group = factor(
      stock_group, 
      labels = c(
        "other", "Col_Spr", "Col_Sum/Fall", "PSD", "WCVI", "ECVI_SOMN",
        "FR_Spr_4sub2", "FR_Spr_5sub2", "FR_Sum_5sub2", "FR_Sum_4sub1", 
        "FR_Fall"
      ))
  ) %>% 
  # remove samples with greater than one year difference 
  filter(
    !abs(bias_n) > 1
  )

samp_size <- dum %>% 
  group_by(stock_group, est_age) %>% 
  tally()

ppn_dat <- dum %>%  
  group_by(stock_group, est_age) %>% 
  mutate(age_n = n(), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(stock_group, bias, est_age, age_n) %>% 
  tally() %>% 
  mutate(prop = n / age_n)


obs_error <- ggplot(ppn_dat) +
  geom_bar(aes(x = est_age, y = prop, fill = bias),
           position="stack", stat="identity") +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  geom_text(
    data = samp_size, aes(x = est_age, y = -Inf, label = paste(n)),
    vjust = -1.1, size = rel(2.5)
  ) +
  labs(x = "Estimated Total Age", y = "Proportion of Samples") +
  theme(
    legend.position = "top"
  )


png(
  here::here("figs", "stock_size_age", "age_error.png"),
  height = 5.5, width = 6.5, units = "in", res = 250
)
obs_error
dev.off()


## fit multinomial model
prior_in <- c(
  prior(normal(-1, 1), class = "Intercept", dpar = "muover"),
  prior(normal(-1, 1), class = "Intercept", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "est_age3", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "est_age4", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "est_age5", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "est_age6", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "est_age3", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "est_age4", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "est_age5", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "est_age6", dpar = "muover"),
  prior(exponential(1), class = "sd", group = "stock_group", dpar = "muover"),
  prior(exponential(1), class = "sd", group = "stock_group", dpar = "muunder")
)

fit <- brm(
  formula = bf(bias ~ est_age + (1 | stock_group),
               family = categorical(link = "logit")),
  data = dum,
  prior = prior_in,
  chains = 4, cores = 4, iter = 2000,
  control = list(adapt_delta = 0.96)
)
saveRDS(fit, here::here("data", "model_fits", "age_error_est.rds"))
fit <- readRDS(here::here("data", "model_fits", "age_error_est.rds"))


new_data <- expand.grid(
  est_age = levels(dum$est_age),
  stock_group = levels(dum$stock_group)
) %>% 
  mutate(
    asg = paste(est_age, stock_group, sep = "-")
  )

# Generate posterior predicted probabilities for the new data
posterior_probs_new <- posterior_epred(fit, newdata = new_data)
colnames(posterior_probs_new) <- new_data$asg

posterior_probs_new %>%
  apply(c(2, 3), mean) %>%  # Take mean across posterior samples
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "bias", 
               values_to = "mean_prob") %>%
  mutate(asg = rep(new_data$asg, each = length(unique(dum$bias))),
         bias = factor(bias, levels = c("zero", "under", "over"))) %>% 
  left_join(., new_data, by = "asg") %>% 
  ggplot(.) +
  geom_bar(aes(x = est_age, y = mean_prob, fill = bias),
           position="stack", stat="identity") +
  facet_wrap(~stock_group) +
  ggsidekick::theme_sleek()


# export posterior samples
bias_names <- dimnames(posterior_probs_new)[[3]]
post_list <- vector(mode = "list", length = 3)
for (i in 1:3) {
  post_list[[i]] <- posterior_probs_new[ , , i] %>% 
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "asg", 
                 values_to = "prob") %>% 
    left_join(., new_data, by = "asg") %>%
    mutate(iter = rep(seq(1, 55, by = 1), times = 4000),
           bias = bias_names[[i]]) 
}

post_out <- bind_rows(post_list)
saveRDS(post_out, here::here("data", "rec", "age_bias_post_draws.rds"))


# same as above but ignore stock group in predictions so that average bias 
# is used when a given age/stock group was not observed
new_data2 <- expand.grid(
  est_age = levels(dum$est_age),
  stock_group = levels(dum$stock_group)[1]
)

# Generate posterior predicted probabilities for the new data
posterior_probs_new2 <- posterior_epred(fit, newdata = new_data2, re_formula = NA)
colnames(posterior_probs_new2) <- new_data2$est_age

posterior_probs_new2 %>%
  apply(c(2, 3), mean) %>%  # Take mean across posterior samples
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "bias", 
               values_to = "mean_prob") %>%
  mutate(est_age = rep(new_data2$est_age, each = length(unique(dum$bias))),
         bias = factor(bias, levels = c("zero", "under", "over"))) %>%
  # left_join(., new_data2, by = "asg") %>% 
  ggplot(.) +
  geom_bar(aes(x = est_age, y = mean_prob, fill = bias),
           position="stack", stat="identity") +
  # facet_wrap(~stock_group) +
  ggsidekick::theme_sleek()


# export posterior samples
bias_names <- dimnames(posterior_probs_new2)[[3]]
post_list2 <- vector(mode = "list", length = 3)
for (i in 1:3) {
  post_list2[[i]] <- posterior_probs_new2[ , , i] %>% 
    as.data.frame() %>%
    pivot_longer(cols = everything(), names_to = "age", 
                 values_to = "prob") %>% 
    mutate(iter = rep(seq(1, 4000, by = 1), each = 5),
           bias = bias_names[[i]]) 
}

post_out2 <- bind_rows(post_list2)
saveRDS(post_out2, 
        here::here("data", "rec", "age_bias_post_draws_no_stock.rds"))




## simulation check

# 1. Simulate Multinomial Response Data with Categorical Predictor
set.seed(123)  # For reproducibility

# Parameters
n <- 500  # Number of observations
k <- 3    # Number of categories (response)
levels_x <- 3  # Number of levels for the categorical predictor

# Generate categorical predictor (factor)
x <- factor(sample(1:levels_x, size = n, replace = TRUE), labels = c("Level1", "Level2", "Level3"))

# Define true coefficients for each level of `x` (reference level for x is "Level1")
beta <- matrix(c(-2, -2,      # Coefficients for Level1 (reference)
                 -1, -1,   # Coefficients for Level2
                 -1, -1), # Coefficients for Level3
               ncol = k - 1, byrow = TRUE)

# Design matrix for the categorical predictor
# Each row of beta corresponds to a level of x, excluding the reference level
design_matrix <- model.matrix(~ x - 1)

# Compute linear predictors for categories 2 and 3
eta <- design_matrix %*% beta

# Add baseline category (category 1, eta = 0)
eta_full <- cbind(0, eta)

# Convert to probabilities using the softmax function
probs <- exp(eta_full) / rowSums(exp(eta_full))

# Simulate categorical response
y <- apply(probs, 1, function(p) sample(1:k, size = 1, prob = p))
y <- factor(y, labels = c("Cat1", "Cat2", "Cat3"))

# Combine into a data frame
sim_data <- data.frame(y, x)

sim_data %>% 
  group_by(x, y) %>% 
  tally()

fit <- brm(
  formula = bf(y ~ x, #+ (1 | stock_group), 
               family = categorical(link = "logit")),
  data = sim_data,
  chains = 4, cores = 4, iter = 2000
)


# 3. Summarize the Model
summary(fit)

# 4. Compare True and Estimated Coefficients
# Extract estimated coefficients
estimated_beta <- summary(fit)$coefficients$cond

cat("\nTrue Coefficients:\n")
print(beta)
cat("\nEstimated Coefficients:\n")
print(estimated_beta)
