### MVN Random Intercepts Dirichlet
## Adapt stockseasonr model with MVN intercepts
## NOTE: attempted to incorporate random smooths, but unclear how to proceed 
## given parameters are a matrix, not vector
## Nov. 29, 2021
## Updated March 19 to explore composition at area levels

library(tidyverse)
library(TMB)


# tmb models
compile(here::here("src", "negbin_rsplines_dirichlet_mvn.cpp"))
dyn.load(dynlib(here::here("src", "negbin_rsplines_dirichlet_mvn")))


# utility functions for prepping smooths 
source(here::here("R", "functions", "utils.R"))
# data prep and model fitting functions
source(here::here("R", "functions", "fit.R"))



