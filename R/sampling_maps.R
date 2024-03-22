## Sampling Maps
# Maps for locations of SRKW diet and fisheries-dependent data
# Uses cleaned and aggregated data generated in model fitting scripts

library(tidyverse)
library(sf)


diet_dat_in <- readRDS(
  here::here("data", "rkw_diet", "cleaned_ppn_dat.rds")
) 
diet_dat <- diet_dat_in %>% 
  select(utm_x_m = utm_x, utm_y_m = utm_y, era, strata, n_samples) %>% 
  mutate(
    dataset = "diet"
  ) %>% 
  distinct()

rec_dat_in <- readRDS(
  here::here("data", "rec", "cleaned_ppn_data_rec_xy.rds")
) 
rec_dat <- rec_dat_in %>% 
  # filter(
  #   !(strata == "Nitinat" & utm_x_m < 35000)
  # ) %>% 
  mutate(
    era = "current",
    dataset = "fishery",
    n_samples = round(sample_id_n, digits = 0)
  ) %>% 
  select(utm_x_m, utm_y_m, era, strata, n_samples, dataset)


# combine and trim for mapping
dat <- rbind(diet_dat, rec_dat) 


strata_colour_pal <- RColorBrewer::brewer.pal(
  length(levels(dat$strata)),
  "Paired"
)
names(strata_colour_pal) <- levels(dat$strata) 
era_shape_pal <- c(15, 16)
names(era_shape_pal) <- levels(dat$era) 

  
coast <- rbind(rnaturalearth::ne_states( "United States of America", 
                                           returnclass = "sf"), 
                 rnaturalearth::ne_states( "Canada", returnclass = "sf")) %>% 
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
  geom_sf(data = coast, color = "black", fill = "grey") +
  geom_point(
    data = dat,
    aes(x = utm_x_m, y = utm_y_m, colour = strata, shape = era, 
        size = n_samples), 
    alpha = 0.4
  ) +
  coord_sf(expand = FALSE) +
  scale_colour_manual(values = strata_colour_pal, name = element_blank()) +
  scale_size_continuous(name = "Number\nof\nSamples") +
  scale_shape_manual(values = era_shape_pal, name = element_blank()) +
  ggsidekick::theme_sleek() +
  theme(
    panel.background = element_rect(fill = "white"),
    axis.title = element_blank(),
    legend.position = "top"
  ) +
  facet_wrap(
    ~dataset,
    ncol = 1
  )


png(
  here::here("figs", "ms_figs", "study_area.png"),
  height = 7.5, width = 7.5, units = "in", res = 250
)
diet_samp_map
dev.off()



