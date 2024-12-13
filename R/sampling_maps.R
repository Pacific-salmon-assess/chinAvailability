## Sampling Maps
# Maps for locations of SRKW diet and fisheries-dependent data
# Uses cleaned and aggregated data generated in model fitting scripts

library(tidyverse)
library(sf)


diet_dat_in <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) 
diet_dat <- diet_dat_in$stock %>% 
  mutate(
    dataset = "Diet",
    utm_x_m = utm_x * 1000,  
    utm_y_m = utm_y * 1000 
  ) %>% 
  group_by(utm_x_m, utm_y_m, strata, era, dataset) %>% 
  summarize(
    n_samples = round(sum(n_samples), digits = 0)
  ) %>% 
  ungroup() %>% 
  select(utm_x_m, utm_y_m, era, strata, n_samples, dataset) %>% 
  distinct()

rec_dat_in <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
) 
rec_dat <- rec_dat_in %>% 
  group_by(utm_x_m, utm_y_m, strata) %>% 
  summarize(
    era = "current",
    dataset = "Fishery",
    n_samples = round(sum(agg_prob), digits = 0)
  ) %>% 
  select(utm_x_m, utm_y_m, era, strata, n_samples, dataset) %>% 
  distinct() %>% 
  ungroup()


# combine and trim for mapping
dat <- rbind(diet_dat, rec_dat) 


# make a strata key that defines "representative" location within each
# strata so that predictions can be made for both datasets
strata_key <- data.frame(
  strata = levels(dat$strata),
  lon = c(-124.92151, -124.85158, -124.52353, -124.12741, -123.56583, 
          -123.06742, -123.35288, -123.50943),
  lat = c(48.55033, 48.65752, 48.53, 48.38, 48.29965,
          48.45791, 48.70844, 48.67728),
  dataset = "Fishery"
) %>% 
  sdmTMB::add_utm_columns(
    ., ll_names = c("lon", "lat"), ll_crs = 4326, units = "m",
    utm_names = c("utm_x_m", "utm_y_m")
  ) %>% 
  filter(
    !strata %in% c("San Juan\nIslands", "Central\nJDF")
  )
# saveRDS(
#   strata_key,
#   here::here("data", "spatial", "strata_key.rds")
# )


strata_colour_pal <- c(
  "#e41a1c",  "#377eb8",  "#4daf4a",  "#984ea3",  "#ff7f00",  "#ffff33", 
  "#a65628",  "#f781bf"
  )
names(strata_colour_pal) <- levels(dat$strata) 
era_shape_pal <- c(22, 21)
names(era_shape_pal) <- levels(dat$era) 

  
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")
  ) %>% 
  sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>%
  sf::st_transform(., crs = sf::st_crs("+proj=utm +zone=10 +units=m")) %>% 
  sf::st_crop(
    ., 
    xmin = min(dat$utm_x_m) - 1500, 
    ymin = min(dat$utm_y_m) - 2000,
    xmax = max(dat$utm_x_m) + 1500, 
    ymax = max(dat$utm_y_m) + 2000
  )

diet_samp_map <- ggplot() +
  geom_sf(data = coast,  fill = "white", colour = "white") +
  geom_point(
    data = dat,
    aes(x = utm_x_m, y = utm_y_m, fill = strata, shape = era, 
        size = n_samples), 
    alpha = 0.7
  ) +
  geom_point(
    data = strata_key, aes(x = utm_x_m, y = utm_y_m), 
    shape = 4, stroke = 1.5
  ) +
  coord_sf(expand = FALSE) +
  scale_fill_manual(values = strata_colour_pal,
                    name = "Spatial\nStrata") +
  scale_size_continuous(guide = "none") +
  scale_shape_manual(values = era_shape_pal, 
                     name = "Diet\nSampling\nEra") +
  ggsidekick::theme_sleek() +
  guides(
    fill = guide_legend(override.aes = list(shape = 21))
  ) +
  theme(
    strip.background = element_rect(colour="white", fill="white"),
    panel.background = element_rect(fill = "grey40"),
    legend.position = "top",
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
  ) +
  facet_wrap(
    ~dataset, ncol = 1
  )


png(
  here::here("figs", "sampling_map.png"),
  height = 7.5, width = 8, units = "in", res = 250
)
diet_samp_map
dev.off()
