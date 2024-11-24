# Calculate Age Classification Errror
# Uses data exported from clean_data.R
# Nov 22, 2024
#1) make age_gr vs age key 
#2) assign scale_age based on age_gr and resolved_age
#3) define true age based on cwt, or pbt
#4) subset to samples that have both a true age and scale age estimate
#5) calculate classification accuracy

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

dum <- age_dat %>% 
  left_join(., age_key, by = "age_gr") %>% 
  mutate(
    bias = case_when(
      age - est_age == 0 ~ "zero",
      age - est_age > 0 ~ "under",
      age - est_age < 0 ~ "over"
    ) %>% 
      factor(., levels = c("zero", "under", "over")),
    age = as.factor(age)
  )

samp_size <- dum %>% 
  group_by(stock_group, age) %>% 
  summarize(
    n = length(unique(biokey))
  )

ppn_dat <- dum %>% 
  group_by(age) %>% 
  mutate(age_n = n(), .groups = "drop") %>% 
  ungroup() %>% 
  group_by(bias, age, age_n) %>% 
  tally() %>% 
  mutate(prop = n / age_n)

ggplot(ppn_dat) +
  geom_bar(aes(x = age, y = prop, fill = bias),
           position="stack", stat="identity") +
  # facet_wrap(~stock_group) +
  ggsidekick::theme_sleek() +
  geom_text(
    data = samp_size, aes(x = age, y = -Inf, label = paste(n)),
    vjust = -1.1
  )


## fit initial multinomial model
prior_in <- c(
  prior(normal(-1, 1), class = "Intercept", dpar = "muover"),
  prior(normal(-1, 1), class = "Intercept", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "age3", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "age4", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "age5", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "age6", dpar = "muunder"),
  prior(normal(0, 2), class = "b", coef = "age3", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "age4", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "age5", dpar = "muover"),
  prior(normal(0, 2), class = "b", coef = "age6", dpar = "muover"),
  prior(exponential(1), class = "sd")  # Prior for random effects
)

fit <- brm(
  formula = bf(bias ~ age + (1 | stock_group), 
               family = categorical(link = "logit")),
  data = dum,
  prior = prior_in,
  chains = 4, cores = 4, iter = 2000
)


new_data <- data.frame(
  age = levels(dum$age)
)

# Generate posterior predicted probabilities for the new data
posterior_probs_new <- posterior_epred(fit, newdata = new_data)

posterior_probs_new %>%
  apply(c(2, 3), mean) %>%  # Take mean across posterior samples
  as.data.frame() %>%
  pivot_longer(cols = everything(), names_to = "bias", 
               values_to = "mean_prob") %>%
  mutate(age = rep(new_data$age, each = length(unique(dum$bias)))) %>% 
  ggplot(.) +
  geom_bar(aes(x = age, y = mean_prob, fill = bias),
           position="stack", stat="identity") +
  ggsidekick::theme_sleek()



# Load necessary library
library(glmmTMB)

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


# 2. Fit the Multinomial Model using glmmTMB
fit <- glmmTMB(
  y ~ x,  # Formula for fixed effect
  family = multinomial(link = "logit"),
  data = sim_data
)
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
