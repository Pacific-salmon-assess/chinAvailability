## Prep Spatial Data
# Import subarea shapefiles and calculate centroid for each; then calculate 
# distance to arbitrary point outside of Juan de Fuca Strait as proxy for
# dispersal distance

library(sf)
library(ggplot2)
library(tidyverse)


# parallelize based on operating system (should speed up some of the spatial
# processing calculations)
library("parallel")
ncores <- detectCores() - 2
if (Sys.info()['sysname'] == "Windows") {
  library("doParallel")
  cl <- makeCluster(ncores)
  registerDoParallel(cl)
  shp_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/creel_areas/"
} else {
  doMC::registerDoMC(ncores)
  shp_path <- "/Users/cam/Google Drive/spatial/creel_areas/"
}


creel_utm <- st_read(here::here(shp_path, "creelareaspfma_2021.shp")) %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>% 
  janitor::clean_names()

# calculate centroid 
sf_cent <- st_centroid(creel_utm)

ggplot() +
  geom_sf(data = creel_utm)


# subset and check
creel_utm2 <- creel_utm %>% 
  filter(statarea %in% c(123, 121, 21, 20, 19, 18, 29, 28)) 
cent_inside <- sf_cent %>% 
  filter(statarea %in% c(121, 21, 20, 19, 18, 29, 28)) 
cent_start <- sf_cent %>% 
  filter(statarea %in% c(123), creelsub == "123-I") 

ggplot() +
  geom_sf(data = creel_utm2 %>% filter(statarea == "123"), 
          aes(fill = creelsub)) +
  geom_sf(data = cent_inside, col = "red")


# calculate distance and convert to km
swvi_dist <- st_distance(cent_inside, cent_start) %>% 
  as.numeric() 
cent_inside$dist_123i <- swvi_dist / 1000

saveRDS(cent_inside, here::here("data", "rec", "creel_subarea_spatial.rds"))
