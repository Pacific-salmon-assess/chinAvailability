## Sample Coverage
# Explore alternative strata assignments


library(tidyverse)
library(maptools)
library(rmapshaper)
library(mapdata)
library(cluster)
library(factoextra)


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
    strata = ifelse(grepl("n_haro", strata), "n_haro", strata),
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
    ) %>%
      fct_reorder(., lon)
  )

# strata = original designation based on PFMA and overlap w habitat
# strata2 = pooling based on eyeball + overlap w/ habitat
# strata3 = pooling based on eyeball + no difference w/ habitat

## CLUSTER ANALYSIS ------------------------------------------------------------

## alternative means of identifying strata based on clustering
rec_trim <- rec_raw %>% 
  filter(
    fl > 600 | legal == "legal",
    month_n > 5 & month_n < 9
  ) %>% 
  mutate(
    # split out columbia again for greater resolution
    stock_group = case_when(
      grepl("CR-", pst_agg) ~ "Col",
      # pst_agg == "WACST" ~ "other",
      grepl("Spring", stock_group) | stock_group == "Fraser_Summer_5.2" ~ 
        "FR_yearling",
      TRUE ~ stock_group
    )
  ) %>% 
  group_by(strata3, fishing_site, lat, lon, stock_group) %>% 
  summarize(
    sum_prob = sum(prob)
  ) %>% 
  group_by(strata3, fishing_site, lat, lon) %>%
  mutate(
    n_samps = sum(sum_prob),
    ppn = sum_prob / n_samps
  ) %>% 
  ungroup() %>% 
  # remove sites with fewer than 10 samples
  filter(
    !n_samps < 10
  ) %>%
  # remove sum_prob column so pivot possible
  select(
    -sum_prob, -n_samps
  ) %>% 
  pivot_wider(
    names_from = stock_group, values_from = ppn
  ) 


rec_mat <- rec_trim %>% 
  select(-fishing_site, -strata3) %>% 
  as.matrix()

# replace nil observations with 0s
rec_mat[is.na(rec_mat)] <- 0

rec_mat_scaled <- scale(rec_mat)


## explore number of clusters using two different distance matrices
# 1) includes spatial location plus stock comp; centered and scaled to make
# comparable
# 2) includes only stock comp with Bray curtis dissimilarity
fviz_nbclust(rec_mat[ , -c(1:2)], kmeans, method = "silhouette")
fviz_nbclust(rec_mat_scaled, kmeans, method = "silhouette")

dist_matrix <- dist(rec_mat_scaled)
hclust_fit <- hclust(dist_matrix, method = "ward.D2")
hclust_fit$labels <- rec_trim$strata3
plot(hclust_fit, cex = 0.6, hang = -1)


dist_matrix_ppns <- vegan::vegdist(rec_mat[ , -c(1:2)], method = "bray")
hclust_fit2 <- hclust(dist_matrix_ppns, method = "ward.D2")
hclust_fit2$labels <- rec_trim$strata3
plot(hclust_fit2, cex = 0.6, hang = -1)


# add cluster IDs
rec_trim$cluster_1 <- cutree(hclust_fit, k = 5) %>% as.factor()
rec_trim$cluster_2 <- cutree(hclust_fit2, k = 5) %>% as.factor()

  

## MAPS ------------------------------------------------------------------------

# all observed locations
obs_stations <- rec_raw %>%
  filter(lat < 48.8) %>% 
  select(id, lat, lon, rkw_habitat, fishing_site, strata, strata2, strata3) %>% 
  distinct() %>% 
  group_by(lat, lon, rkw_habitat, fishing_site, strata, strata2, strata3) %>% 
  tally() %>% 
  left_join(
    ., 
    rec_trim %>% 
      select(fishing_site, cluster_1, cluster_2)
  ) 

pfma_subareas <- readRDS(
  here::here("data", "spatial", "pfma_subareas_sBC.rds")
) 

# map including observed locations
base_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = NA) +
  # geom_sf(data = hab_sf, color = "red") +
  geom_sf(data = pfma_subareas, aes(colour = rkw_overlap), fill = NA) +
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

cluster1 <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = rkw_habitat, fill = cluster_1),
    alpha = 0.7,
    shape = 21
  ) 
cluster2 <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = rkw_habitat, fill = cluster_2),
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


png(here::here("figs", "strata_breakdown", "cluster1.png"), 
    height = 3.5, width = 7, units = "in", res = 250)
cluster1
dev.off()

png(here::here("figs", "strata_breakdown", "cluster2.png"), 
    height = 3.5, width = 7, units = "in", res = 250)
cluster2
dev.off()



## export key used in model fitting
strata_out <- rec_raw %>% 
  select(
    fishing_site, lat, lon, strata, strata2, strata3
  ) %>% 
  distinct() 

saveRDS(
  strata_out,
  here::here(
    "data", "rec", "strata_key.rds"
  )
)
