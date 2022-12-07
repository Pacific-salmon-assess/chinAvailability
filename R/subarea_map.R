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
  st_as_sf() %>% 
  mutate(
    focal_subarea = ifelse(core_area == "yes", subarea, NA)
  )


# import coastline SF dataframe
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84"))

subarea_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_sf(data = dum, aes(fill = focal_subarea)) +
  scale_fill_brewer(type = "qual", palette = "Set1", na.value="grey60",
                    name = "RKW Core\nArea") +
  # scale_fill_discrete(name = "RKW Core\nArea") +
  coord_sf(xlim = c(-126, -122), ylim = c(48, 49.75)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

png(here::here("figs", "afs_subarea_preds", "subarea_map.png"),
    height = 5, width = 6.5, units = "in", res = 250)
subarea_map
dev.off()


### NITINAT FOCUSED MAP --------------------------------------------------------

nitinat_areas <- creel_sub %>% 
  filter(
    # subareaid %in% c("23I", "23J", "123I") |
      statarea %in% c("123", "23")
  ) 

ggplot() +
  geom_sf(data = coast, color = "black", fill = "white") +
  geom_sf(data = nitinat_areas,
          aes(fill = as.factor(subareaid))) +
  scale_alpha_manual(values = alpha_pal) +
  coord_sf(xlim = c(-126, -122), ylim = c(48, 49.75)) +
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())
