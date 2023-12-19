## Sample Coverage
# Explore alternative strata assignments


library(tidyverse)
library(maptools)
library(rmapshaper)
library(mapdata)


# map data
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  # sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.75, ymax = 48.8)

hab_sf <- readRDS(
  here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
) %>% 
  # sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.75, ymax = 48.8)


rec_raw <- readRDS(here::here("data", "rec", "rec_gsi.rds")) %>% 
  janitor::clean_names() %>% 
  rename(stock_region = region) %>% 
  filter(!is.na(lat), 
         !is.na(lon),
         #remove due to convergence issues and unique stock comp
         !strata == "saanich",
         !subarea == "19-8") %>% 
  mutate(
    #redefine single station at S border of n_haro as s_haro
    strata = ifelse(grepl("n_haro", strata) & lat < 48.55,
                    "s_haro",
                    strata),
    strata = ifelse(grepl("n_haro", strata), "n_haro", strata) %>% 
      factor(),
    strata2 = case_when(
      lon < -125 & strata != "wcvi_outside" ~ "barkley_corner",
      strata == "swiftsure" & lat > 48.55 ~ "nitinat_midshore",
      # strata == "nitinat_midshore" & lon > -124.82 ~ "renfrew_habitat",
      grepl("nitinat", strata) & lon > -124.83 ~ "renfrew_habitat",
      strata == "victoria" & lon < -123.48 ~ "sooke_nonhabitat",
      strata == "s_haro" ~ "victoria",
      TRUE ~ strata
    ),
    strata3 = case_when(
      grepl("renfrew", strata2) ~ "renfrew",
      grepl("nitinat", strata2) ~ "nitinat",
      grepl("sooke", strata2) ~ "sooke",
      TRUE ~ strata2
    )
  )

# strata = original designation based on PFMA and overlap w habitat
# strata2 = pooling based on eyeball + overlap w/ habitat
# strata3 = pooling based on eyeball + no difference w/ habitat

# all observed locations
obs_stations <- rec_raw %>%
  filter(lat < 48.8) %>% 
  select(id, lat, lon, rkw_habitat, strata, strata2, strata3) %>% 
  distinct() %>% 
  group_by(lat, lon, rkw_habitat, strata, strata2, strata3) %>% 
  tally()


# map including observed locations
base_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = NA) +
  geom_sf(data = hab_sf, color = "red") +
  # geom_sf(data = pfma_subareas, aes(colour = rkw_overlap), fill = NA) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0)) +
  # scale_shape_manual(values = shape_pal) +
  theme(
    legend.position = "top",
    axis.text = element_blank(),
    axis.title = element_blank()
  ) +
  scale_fill_brewer(type = "qual", palette = "Paired") +
  guides(
    fill = guide_legend(
      override.aes = list(shape = 21),
      nrow = 2, byrow = TRUE,
      title = "Habitat\nStrata"
    ),
    size = "none",
    colour = "none",
    shape = "none"
    # fill = "none",
    # size = guide_legend(nrow = 2, byrow = TRUE, title = "GSI\nSample\nSize")#,
    # colour = guide_legend(nrow = 2, byrow = TRUE, 
    #                       title = "PFMA\nIncludes\nPoly."),
    # shape = guide_legend(nrow = 2, byrow = TRUE, title = "Region")
  )


strata1 <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = rkw_habitat, fill = strata),
    alpha = 0.7,
    shape = 21
  ) 
strata2 <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = rkw_habitat, fill = strata2),
    alpha = 0.7,
    shape = 21
  ) 
strata3 <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = rkw_habitat, fill = strata3),
    alpha = 0.7,
    shape = 21
  ) 

png(here::here("figs", "strata_breakdown", "strata1.png"), 
    height = 3.5, width = 7, units = "in", res = 250)
strata1
dev.off()

png(here::here("figs", "strata_breakdown", "strata2.png"), 
    height = 3.5, width = 7, units = "in", res = 250)
strata2
dev.off()

png(here::here("figs", "strata_breakdown", "strata3.png"), 
    height = 3.5, width = 7, units = "in", res = 250)
strata3
dev.off()
