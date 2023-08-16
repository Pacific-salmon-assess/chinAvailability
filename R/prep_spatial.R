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
  shp_path <- "/Users/cam/Google Drive/spatial"
}


# RKW HABITAT USE --------------------------------------------------------------

# import rasters (replace for now with polygon)
# swift_hab <- raster::raster(
#     here::here(shp_path, "srkw_foraging_areas", 
#                "Swiftsure_post_forage_exc_0.25_NAD83_BCAlbers.tif"))
# haro_hab <- raster::raster(
#   here::here(shp_path, "srkw_foraging_areas", 
#              "Haro_post_forage_exc_0.25_NAD83_BCAlbers.tif"))
# hab_list <- list(swift_hab, haro_hab)
# 
# # switch to sf 
# hab_sf <- purrr::map(hab_list, ~ {
#   as(.x, 'SpatialGridDataFrame') %>%
#     as.data.frame() %>%
#     rename(lon = s1, lat = s2, for_prob = layer) %>%
#     st_as_sf(., coords = c("lon", "lat"),
#              # specify Albers based on readme
#              crs = sp::CRS("+init=epsg:3005"))
# }) %>% 
#   bind_rows()


# shapefile boundaries (most conservative 70% boundary based on posterior draws)
swift_shp <- st_read(
  here::here(shp_path, "srkw_foraging_areas", 
             "swiftsure.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp"))
haro_shp <- st_read(
  here::here(shp_path, "srkw_foraging_areas", 
             "haro.forage.0.25exc.0.7prop.poly_NAD83_BCAlbers.shp"))
hab_sf <- rbind(swift_shp, haro_shp) %>%
  sf::st_transform(., crs = sp::CRS("+init=epsg:3005"))


# coastline cropped 
coast_albers <- rbind(rnaturalearth::ne_states( "United States of America", 
                                             returnclass = "sf"), 
                   rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
  sf::st_crop(., 
              xmin = -125.5, ymin = 48.15, xmax = -122.25, ymax = 49.25) %>% 
  sf::st_transform(., crs = sp::CRS("+init=epsg:3005"))


# recreational fishing data
rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region, region = cap_region)  %>% 
  mutate(
    lon = ifelse(lon > 0, -1 * lon, lon)
  ) %>% 
  filter(
    !is.na(lat)
  )

# convert to sf and add presence in cutoff
rec_sf <- rec_raw %>% 
  st_as_sf(
    ., 
    coords = c("lon", "lat"), 
    crs = sp::CRS("+proj=longlat +datum=WGS84")
  ) %>% 
  sf::st_transform(., crs = sp::CRS("+init=epsg:3005"))

rec_sf_sub <- rec_sf %>% 
  # select pts in habitat
  st_intersection(., hab_sf)

rec_out <- rec_sf %>% 
  mutate(
    rkw_habitat = ifelse(id %in% rec_sf_sub$id, "yes", "no")
  ) 

ggplot() + 
  ggsidekick::theme_sleek() +
  theme(axis.title = element_blank()) +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # geom_sf(data = hab_sf %>% filter(for_prob > 0.7),
  #         aes(colour = for_prob), shape = 1) +
  geom_sf(data = hab_sf, color = "red", fill = NA, size = 1.25) +
  geom_sf(data = coast_albers) +
  geom_jitter(data = rec_out,
              aes(x = x, y = y, fill = rkw_habitat),
              shape = 23, alpha = 0.05, width = 0.05) 


## convert back to WGS84 and export
saveRDS(
  hab_sf %>% 
    st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")),
  here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
)

saveRDS(
  rec_out %>% 
    st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
    sfheaders::sf_to_df(., fill = TRUE) %>% 
    rename(lon = x, lat = y),
  here::here("data", "rec_gsi_spatial.rds")
)

# # CREEL SUBAREA COV ------------------------------------------------------------
# 
# creel_utm <- st_read(
#   here::here(shp_path, "creel_areas", "creelareaspfma_2021.shp")
#   ) %>% 
#   st_transform(., crs = sp::CRS("+proj=utm +zone=9 +units=m")) %>% 
#   janitor::clean_names()
# 
# # calculate centroid 
# sf_cent <- st_centroid(creel_utm)
# 
# ggplot() +
#   geom_sf(data = creel_utm %>% 
#             filter(statarea %in% c("20", "18", "19", "29", "121",
#                                    "21")),
#           aes(fill = subareaid))
# 
# 
# # subset and check
# creel_utm2 <- creel_utm %>% 
#   filter(statarea %in% c(123, 121, 21, 20, 19, 18, 29, 28)) 
# cent_inside <- sf_cent %>% 
#   filter(statarea %in% c(121, 21, 20, 19, 18, 29, 28)) 
# cent_start <- sf_cent %>% 
#   filter(statarea %in% c(123), creelsub == "123-I") 
# 
# ggplot() +
#   geom_sf(data = creel_utm2 %>% filter(statarea == "123"), 
#           aes(fill = creelsub)) +
#   geom_sf(data = cent_inside, col = "red")
# 
# 
# # calculate distance and convert to km
# swvi_dist <- st_distance(cent_inside, cent_start) %>% 
#   as.numeric() 
# cent_inside$dist_123i <- swvi_dist / 1000
# 
# saveRDS(cent_inside, here::here("data", "rec", "creel_subarea_spatial.rds"))
