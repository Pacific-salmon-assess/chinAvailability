## Prep Spatial Data
# 1) Import subarea shapefiles and calculate centroid for each; then calculate 
# distance to arbitrary point outside of Juan de Fuca Strait as proxy for
# dispersal distance
# 2) Import and clean RKW habitat use predictions (shape files generated using
# J. Watson's model and provided by S. Toews)

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
  shp_path <- "C:/Users/FRESHWATERC/Documents/drive docs/spatial/"
} else {
  doMC::registerDoMC(ncores)
  shp_path <- "/Users/cam/Google Drive/spatial/"
}


# RKW HABITAT USE --------------------------------------------------------------

srkw_hab <- st_read(
  here::here(shp_path, "srkw_foraging_areas", "creelareaspfma_2021.shp")
)

swift_hab <- raster::raster(
    here::here(shp_path, "srkw_foraging_areas", 
               "Swiftsure_post_forage_exc_0.25_NAD83_BCAlbers.tif"))
haro_hab <- raster::raster(
  here::here(shp_path, "srkw_foraging_areas", 
             "Haro_post_forage_exc_0.25_NAD83_BCAlbers.tif"))
hab_list <- list(swift_hab, haro_hab)

hab_full <- raster::merge(swift_hab, haro_hab, tolerance = 0.1)


template <- projectRaster(from = swift_hab, to= haro_hab, alignOnly=TRUE)
#template is an empty raster that has the projected extent of r2 but is aligned with r1 (i.e. same resolution, origin, and crs of r1)
r2_aligned<- projectRaster(from = swift_hab, to = template)
r_merged <- merge(haro_hab, r2_aligned) 
r_merged2 <- mosaic(r1,r2_aligned, fun=mean, na.rm=TRUE)


# switch to sf 
hab_sf <- purrr::map(hab_list, ~ {
  as(.x, 'SpatialGridDataFrame') %>%
    as.data.frame() %>%
    rename(lon = s1, lat = s2, for_prob = layer) %>%
    st_as_sf(., coords = c("lon", "lat"),
             # specify Albers based on readme
             crs = sp::CRS("+init=epsg:5070"))
}) %>% 
  bind_rows()

coast_albers <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -127.5, ymin = 46, xmax = -122, ymax = 49.5) %>% 
  sf::st_transform(., crs = sp::CRS("+init=epsg:5070"))

ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  geom_sf(data = hab_sf, 
          aes(fill = for_prob), shape = 21) +
  geom_sf(data = coast_utm) +
  scale_fill_viridis_c(name = "Mean Depth",
                       direction = -1) +
  theme(legend.position = "top")


# CREEL SUBAREA COV ------------------------------------------------------------

creel_utm <- st_read(
  here::here(shp_path, "creel_areas", "creelareaspfma_2021.shp")
  ) %>% 
  st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>% 
  janitor::clean_names()

# calculate centroid 
sf_cent <- st_centroid(creel_utm)

ggplot() +
  geom_sf(data = creel_utm %>% 
            filter(statarea %in% c("20", "18", "19", "29", "121",
                                   "21")),
          aes(fill = subareaid))


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
