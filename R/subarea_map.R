## Study Area Map
# Make map of coastal zone and study areas based on PFMA subarea shapefile

library(sf)
library(ggplot2)
library(tidyverse)
library(maptools)
library(rmapshaper)
library(mapdata)


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

# import stock key generated in dirichlet fit
area_key <- readRDS(here::here("data", "rec", "subarea_key.RDS"))

creel_sub <- st_read(here::here(shp_path, "creelareaspfma_2021.shp")) %>% 
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>%
  janitor::clean_names() %>% 
  # adjust so that subareas in 20D are pooled
  mutate(
    subareaid = ifelse(grepl("20D", subareaid), "20D", subareaid)
  )

dum <- left_join(area_key, 
                 creel_sub %>% 
                   dplyr::select(subarea_original = subareaid, creelsub) ,
          by = "subarea_original") %>% 
  st_as_sf()


# import coastline SF dataframe
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

subarea_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_sf(data = dum, aes(fill = core_area)) +
  scale_fill_discrete(name = "RKW Core\nArea") +
  coord_sf(xlim = c(-126, -122), ylim = c(48, 49.75)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

pdf(here::here("figs", "afs_subarea_preds", "subarea_map.pdf"),
    height = 5, width = 6.5)
subarea_map
dev.off()
