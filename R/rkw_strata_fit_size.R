## Model Fitting
# Fit data to estimate seasonal changes in size composition 
# Extension of rkw_habitat_fit including more granular spatial strata
# 1) Fit 3 models that include different spatial strata 


library(tidyverse)
library(stockseasonr)

source(here::here("R", "utils.R"))

rec_raw <- readRDS(here::here("data", "rec", "rec_size.rds")) %>% 
  janitor::clean_names() 

rec_trim <- rec_raw %>% 
  filter(
    !is.na(lat),
    !is.na(lon),
    #exclude samples collected outside areas in relatively close proximity to 
    # SRKW foraging areas
    !rkw_habitat == "outside",
    !strata == "saanich",
    !is.na(strata),
    !size_bin2 == "sublegal"
  ) %>%
  mutate(
    sample_id = paste(strata, week_n, year, sep = "_"),
    #redefine single station at S border of n_haro as s_haro
    strata = case_when(
      grepl("haro", strata) & lat < 48.55 ~ "victoria",
      grepl("n_haro", strata) ~ "n_haro",
      grepl("nitinat", strata) ~ "nitinat",
      TRUE ~ strata
    ),
    strata_region = ifelse(
      grepl("vic", strata) | grepl("sooke", strata) |
        grepl("n_haro", strata) | grepl("s_haro", strata),
      "east",
      "west"
    ) %>% 
      factor(),
    size_bin2 = paste("sizebin", size_bin2, sep = "-")
  )

size_in <- rec_trim  %>% 
  group_by(sample_id) %>% 
  mutate(
    nn = length(unique(id)) %>% as.numeric
  ) %>% 
  ungroup() %>% 
  group_by(sample_id, 
           strata, strata_region,
           week_n, month_n, year, nn, 
           size_bin2) %>% 
  summarize(prob = length(unique(id)) %>% as.numeric,
            .groups = "drop") 


## MAP OF LOCATIONS ------------------------------------------------------------

# map data
coast <- readRDS(
  here::here("data", "spatial", "coast_major_river_sf_plotting.RDS")) %>% 
  # sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.9, ymax = 49)

hab_sf <- readRDS(
  here::here("data", "spatial", "rkw_critical_habitat_0.25exc_0.7prop.rds")
) %>% 
  # sf::st_transform(., crs = sp::CRS("+proj=longlat +datum=WGS84")) %>% 
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122.9, ymax = 49)

pfma_subareas <- readRDS(
  here::here("data", "spatial", "pfma_subareas_sBC.rds")
) %>% 
  sf::st_crop(xmin = -126, ymin = 48.2, xmax = -122, ymax = 49)


# all observed locations
obs_stations <- rec_trim %>%
  select(id, lat, lon, strata_region, strata) %>% 
  distinct() %>% 
  group_by(lat, lon, strata_region, strata) %>% 
  tally()

# map including observed locations
base_map <- ggplot() +
  geom_sf(data = coast, color = "black", fill = NA) +
  geom_sf(data = hab_sf, color = "red") +
  geom_sf(data = pfma_subareas, aes(colour = rkw_overlap), fill = NA) +
  ggsidekick::theme_sleek() +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0))

shape_pal <- c(21, 22)
names(shape_pal) <- unique(obs_stations$strata_region)
fill_pal <- pals::polychrome(n = length(unique(obs_stations$strata)))
names(fill_pal) <- unique(obs_stations$strata)

samp_locs <- base_map +
  geom_point(
    data = obs_stations, 
    aes(x = lon, y = lat, size = n, shape = strata_region, fill = strata),
    alpha = 0.6
  ) +
  scale_shape_manual(values = shape_pal) +
  scale_fill_manual(values = fill_pal) +
  theme(
    legend.position = "top"
  ) +
  guides(
    fill = "none",
    shape = "none",
    size = guide_legend(nrow = 2, byrow = TRUE, title = "GSI\nSample\nSize")
  )

samp_locs


## MODEL FIT  ------------------------------------------------------------------

size_in %>% group_by(strata_region, month_n) %>% tally()

# make distinct geographic datasets for model fitting
swift_dat <- size_in %>% 
  filter(strata_region == "west",
         #remove may due to low sample size
         !month_n < "6") %>% 
  droplevels()
east_dat <- size_in %>% 
  filter(strata_region == "east") %>% 
  droplevels()


east_dat %>% 
  group_by(strata, size_bin2, month_n) %>% 
  summarize(nn = sum(prob)) %>% 
  ggplot(.) +
  geom_bar(aes(x  = month_n, y = nn, fill = size_bin2), stat = "identity") +
  facet_wrap(~strata)


# make tibble
dat_tbl <- tibble(
  group = c("swiftsure", "east"),
  data = list(swift_dat, east_dat)
) %>% 
  mutate(
    pred_dat = purrr::map(
      data,
      ~ expand.grid(
        strata = unique(.x$strata),
        week_n = seq(min(.x$week_n), max(.x$week_n), by = 0.1),
        strata_region = unique(.x$strata_region)
      ) 
    )
  )


## stockseasonr fits -----------------------------------------------------------

# number of knots only has an impact when poorly sampled months (May) are included
swift_fit <- fit_stockseasonr(
  comp_formula = size_bin2 ~ 1 + strata + 
    s(week_n, bs = "tp", k = 4, m = 2) +
    (1 | year),
  comp_dat = dat_tbl$data[[1]],
  pred_dat = dat_tbl$pred_dat[[1]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)

east_fit <- fit_stockseasonr(
  comp_formula = size_bin2 ~ 1 + strata +
    s(week_n, bs = "cc", k = 4, m = 2) +
    (1 | year),
  comp_dat = dat_tbl$data[[2]],
  pred_dat = dat_tbl$pred_dat[[2]],
  model = "dirichlet",
  random_walk = FALSE,
  fit = TRUE,
  newton_loops = 1,
  silent = FALSE
)


# make predictions
pred_swift <- clean_pred_foo_size(fit = swift_fit, preds = dat_tbl$pred_dat[[1]])
pred_east <- clean_pred_foo_size(fit = east_fit, preds = dat_tbl$pred_dat[[2]])

preds <- rbind(
  pred_swift$preds %>% 
    mutate(region = "swift"),
  pred_east$preds %>%
    mutate(region = "east")
  ) %>% 
  mutate(
    size_bin2 = factor(
      size_bin2,
      levels = c("small", "medium", "large")
    )
  ) 
obs <- rbind(
  pred_swift$obs_dat %>% 
    mutate(region = "swift"),
  pred_east$obs_dat %>%
    mutate(region = "east")
) %>% 
  mutate(
    size_bin2 = factor(
      size_bin2,
      levels = c("small", "medium", "large")
    )
  ) 

colour_pal <- pals::polychrome(n = length(unique(obs_stations$strata)))
names(colour_pal) <- unique(obs_stations$strata)


p <- ggplot(
  data = preds, 
  aes(x = week_n, colour = strata)
) +
  labs(y = "Predicted Size Proportion", x = "Week") +
  facet_grid(region~size_bin2) +
  ggsidekick::theme_sleek() +
  geom_line(aes(y = pred_prob_est)) +
  scale_colour_manual(values = colour_pal)

p_ribbon <- p +
  geom_ribbon(data = preds,
              aes(ymin = pred_prob_low, ymax = pred_prob_up, fill = strata),
              alpha = 0.2) 

p_obs <- ggplot() +
  geom_jitter(data = obs ,
              aes(x = week_n, y = obs_ppn, colour = strata, 
                  size = samp_nn),
              alpha = 0.1) +
  geom_line(data = preds,
            aes(x = week_n, colour = strata, y = pred_prob_est),
            linewidth = 1.25) +
  labs(y = "Predicted Stock Proportion", x = "Month") +
  facet_grid(region~size_bin2) +
  ggsidekick::theme_sleek() +
  scale_size_continuous() +
  theme(legend.position = "top") +
  scale_colour_manual(values = colour_pal)


west_stacked <- ggplot(data = preds %>% 
                         filter(region == "swift"), 
                       aes(x = week_n)) +
  geom_area(aes(y = pred_prob_est, colour = size_bin2, fill = size_bin2), 
            stat = "identity") +
  scale_fill_brewer(name = "Stock", palette = "Spectral") +
  scale_colour_brewer(name = "Stock", palette = "Spectral") +
  labs(y = "Predicted Composition", x = "Week") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size=9),
        plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points")
  ) +
  # coord_cartesian(expand = FALSE, ylim = c(0, NA), xlim = c(1, 12)) +
  facet_wrap(~strata) 

east_stacked <- ggplot(data = preds %>% 
                         filter(region == "east"), 
                       aes(x = week_n)) +
  geom_area(aes(y = pred_prob_est, colour = size_bin2, fill = size_bin2), 
            stat = "identity") +
  scale_fill_brewer(name = "Stock", palette = "Spectral") +
  scale_colour_brewer(name = "Stock", palette = "Spectral") +
  labs(y = "Predicted Composition", x = "Week") +
  theme_classic() +
  theme(legend.position = "right",
        axis.text = element_text(size=9),
        plot.margin = unit(c(2.5, 11.5, 5.5, 5.5), "points")
  ) +
  # coord_cartesian(expand = FALSE, ylim = c(0, NA), xlim = c(1, 12)) +
  facet_wrap(~strata)


pdf(here::here("figs", "rkw_habitat", "size_comp_preds.pdf"),
    height = 7, width = 8.5)
p
p_ribbon
p_obs
west_stacked
east_stacked
dev.off()
