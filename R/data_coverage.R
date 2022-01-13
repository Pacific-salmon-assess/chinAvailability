## Sample Coverage
# Use GSI samples through 2019 to explore composition coverage and determine
# spatial scale of models

library(tidyverse)

rec_raw <- readRDS(here::here("data", "rec", "recIndProbsLong.rds"))


# sample coverage for GSI ------------------------------------------------------

subarea_list <- rec_raw %>%
  group_by(subarea, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>%
  ungroup() %>% 
  split(., .$region)
area_list <- rec_raw %>%
  group_by(area, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>%
  ungroup() %>% 
  split(., .$region)


# visualize coverage by subarea or area
subarea_coverage <- purrr::map2(
  subarea_list, names(subarea_list),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = legal)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_grid(subarea~year) +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)

area_coverage <- purrr::map2(
  area_list, names(area_list),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = legal)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_grid(area~year) +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)

# export
pdf(here::here("figs", "data_coverage", "gsi_samples_subarea.pdf"))
subarea_coverage
dev.off()

pdf(here::here("figs", "data_coverage", "gsi_samples_area.pdf"))
area_coverage
dev.off()



# sample coverage for GSI ------------------------------------------------------

area_list_fl <- rec_raw %>%
  filter(!is.na(fl)) %>% 
  glimpse()
  group_by(area, region, year, month_n, legal) %>%
  summarize(n = length(unique(id)),
            .groups = "drop") %>%
  ungroup() %>% 
  split(., .$region)

area_coverage_fl <- purrr::map2(
  area_list_fl, names(area_list_fl),
  function (x, y) {
    p <- ggplot(x, aes(x = as.factor(month_n), y = n, fill = legal)) +
      geom_bar(position="stack", stat="identity") +
      ggsidekick::theme_sleek() +
      facet_grid(area~year) +
      labs(x = "Month", y = "Samples", title = y)
    print(p)
  }
)

pdf(here::here("figs", "data_coverage", "fl_samples_area.pdf"))
area_coverage_fl
dev.off()