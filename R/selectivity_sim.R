## Selectivity Simulation
# Use estimated parameters from mvtweedie_fit.R to test for evidence of 
# selectivity in RKW diet observations
# March 28, 2024


library(tidyverse)

rkw_dat <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
)