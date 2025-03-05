## Aggregate Abundance
# Figures of escapement/terminal return for relevant stock groups
# Oct 19, 2024


library(tidyverse)


# esc_raw_dat <- read.csv(
#   here::here(
#     "data", "run_size", "run_size_by_stock_group.csv"
#   ),
#   skip = 1
# ) %>% 
#   janitor::clean_names() %>% 
#   pivot_longer(cols = -year, names_to = "stock", values_to = "esc")

# import header key 
esc_raw_dat_headers <- readxl::read_xlsx(
  here::here(
    "data", "run_size", "run_size_by_stock_v2.xlsx"
  ),
  sheet = "header_key"
  ) %>% 
  select(
    stock_group = new_header1, stock = new_header2
  )

esc_dat <- readxl::read_xlsx(
  here::here(
    "data", "run_size", "run_size_by_stock_v2.xlsx"
  ),
  sheet = "summary_data", 
  skip = 1,
  na = "NA"
) %>%
  pivot_longer(cols = -year, names_to = "stock", values_to = "esc") %>% 
  left_join(., esc_raw_dat_headers, by = "stock") %>% 
  mutate(
    stock_group = ifelse(
      stock_group %in% c("washington_coast", "oregon_coast"),
      "wa_or_coast",
      stock_group
    ) %>% 
      factor(
      .,
      levels = c("wa_or_coast", "columbia_spring",
                 "columbia_summer_fall", "puget_sound", "wcvi", "ecvi_somn",
                 "fraser_spring_5sub2", "fraser_summer_5sub2",
                 "fraser_spring_4sub2", "fraser_summer_4sub1", "fraser_fall"),
      labels = c(
        "WA/OR\nCoastal", "Columbia\nSpring", "Columbia\nSummer/Fall",  
        "Puget Sound", "West Coast\nVan. Island", "ECVI\nand SOMN",
        "Fraser\nSpring 4sub2", "Fraser\nSpring 5sub2", "Fraser\nSummer 5sub2", 
        "Fraser\nSummer 4sub1", "Fraser \nFall"
      )
    ),
    region = ifelse(
      stock_group %in% c(
        "Puget Sound", "ECVI\nand SOMN",
        "Fraser\nSpring 4sub2", "Fraser\nSpring 5sub2", "Fraser\nSummer 5sub2", 
        "Fraser\nSummer 4sub1", "Fraser \nFall"
      ),
      "salish_sea",
      "coastal"
    ) %>% 
      factor(
        ., labels = c("Outer Coast", "Salish Sea")
      )#,
    # obs = ifelse(is.na(esc), 0, 1)
  ) %>%
  group_by(stock) %>% 
  filter(
    # remove stocks that are missing esc data after 1983 or 2019
    !any((is.na(esc) & year == "1983") | (is.na(esc) & year == "2019"))#,
    # remove years with gappy data
    # !(year < 1983 | year == "2023")
  ) %>% 
  ungroup() %>% 
  filter(year > 1981)


# # visualize data presence
# ggplot(
#   esc_dat,
#   aes(x = year, y = stock, colour = obs)
# ) +
#   geom_point()


sg_esc <- esc_dat %>% 
  group_by(year, stock_group) %>% 
  summarize(
    sum_esc = sum(esc)
  )  %>% 
  #remove incomplete data (southern US not available)
  filter(
    !is.na(sum_esc)
  ) 

re_esc <- esc_dat %>%
  group_by(year, region) %>% 
  summarize(
    sum_esc = sum(esc)
  ) 


## FIGURES ---------------------------------------------------------------------

smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(esc_dat$stock_group)

sg_plot <- ggplot(sg_esc) +
  geom_line(aes(x = year, y = sum_esc / 1000, colour = stock_group)) +
  scale_colour_manual(values = smu_colour_pal) +
  facet_wrap(~ stock_group, scales = "free_y") +
  ggsidekick::theme_sleek() +
  labs(y = "Terminal Abundance (thousands)") +
  theme(
    axis.title.x = element_blank(),
    legend.position = "none"
  )
  
png(
  here::here("figs", "run_size", "escapement_stock_group.png"),
  height = 5, width = 7.5, units = "in", res = 250
)
sg_plot
dev.off()


reg_plot <- ggplot(re_esc) +
  geom_bar(aes(x = year, y = sum_esc / 1000, fill = region), colour = "black",
           stat = "identity") +
  # scale_alpha_manual(values = average_pal) +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_d() +
  labs(y = "Terminal Abundance (thousands)") +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  labs(fill = NULL)

png(
  here::here("figs", "run_size", "escapement_region.png"),
  height = 5, width = 5.5, units = "in", res = 250
)
reg_plot
dev.off()


## TREND TEST ------------------------------------------------------------------

re_esc$year_z <- scale(re_esc$year) %>% as.numeric()

fit_region <- lm(sum_esc ~ year_z * region, data = re_esc)
summary(fit_region)


## CALCULATE SYNCHONY INDEX ----------------------------------------------------

library(synchrony)

sg_esc_wide <- sg_esc %>% 
  ungroup() %>% 
  select(
    year, sum_esc, stock_group
  ) %>% 
  pivot_wider(names_from = stock_group, values_from = sum_esc) 

esc_mat <- sg_esc_wide %>% 
  select(-year) %>% 
  as.matrix()

rolling_synchrony <- function(data, window_size) {
  n_years <- nrow(data)
  results <- rep(NA, n_years - window_size + 1)  # Store results
  
  for (i in 1:(n_years - window_size + 1)) {
    window_data <- data[i:(i + window_size - 1), ]  # Extract window
    results[i] <- community.sync(window_data, method = "pearson")$obs  # Compute synchrony
  }
  
  return(results)
}


# Calculate rolling 8-year synchrony
sync_values <- rolling_synchrony(esc_mat, window_size = 8)
total_esc <- re_esc %>% 
  summarize(total_esc = sum(sum_esc) / 1000000)
roll_esc <- zoo::rollmean(total_esc$total_esc, k = 8, fill = NA, 
                          align = "right")

synch_dat <- data.frame(
  year = unique(sg_esc$year)[8:length(unique(sg_esc$year))],
  synch = sync_values,
  total_esc = roll_esc[8:length(roll_esc)]
)
max_esc <- max(synch_dat$total_esc)
min_esc <- min(synch_dat$total_esc)
max_synch <- max(synch_dat$synch)
min_synch <- min(synch_dat$synch)

synch_dat$scaled_esc <- (synch_dat$total_esc - min_esc) / (max_esc - min_esc) * 
      (max_synch - min_synch) + min_synch


ggplot(synch_dat, aes(x = year)) +
  geom_line(aes(y = synch), colour = "red") +
  geom_line(aes(y = scaled_esc), colour = "black") +
  scale_y_continuous(
    name = "Synchrony Index (0-1)",  # Left y-axis label
    sec.axis = sec_axis(~ (. - min_synch) / (max_synch - min_synch) * 
                          (max_esc - min_esc) + min_esc,
                        name = "Total Abundance (millions)")  # Right y-axis label and transformation
  ) +
  theme_minimal() +
  theme(
    axis.title.y.right = element_text(color = "black"),  # Right axis color
    axis.title.y.left = element_text(color = "red")  # Left axis color
  ) +
  ggsidekick::theme_sleek()


## CORRELATE WITH WCVI INDEX ---------------------------------------------------

aabm_dat <- read.csv(
  here::here("data", "run_size", "old", "aabm_fishery_abundance_indices.csv"),
  stringsAsFactors = FALSE
) %>% 
  group_by(aabm_fishery) %>% 
  mutate(
    rolling_mean = zoo::rollmean(abundance_index, k = 5, fill = NA, align = "right")
  )


wcvi_aabm <- aabm_dat %>% 
  filter(aabm_fishery == "wcvi") 
nbc_aabm <- aabm_dat %>% 
  filter(aabm_fishery == "nbc") 


dd <- re_esc %>% 
  group_by(year) %>% 
  summarize(total_esc = sum(sum_esc)) %>% 
  left_join(., wcvi_aabm, by = "year") %>% 
  filter(!is.na(aabm_fishery))

png(
  here::here("figs", "run_size", "escapement_aabm_corr_aggregate.png"),
  height = 5, width = 5.5, units = "in", res = 250
)
plot(total_esc ~ abundance_index, data = dd)
dev.off()


cor(dd$total_esc, dd$abundance_index)


sg_esc2 <- sg_esc %>% 
  left_join(., wcvi_aabm, by = "year")  %>% 
  filter(!is.na(aabm_fishery))

cor_df <- sg_esc2 %>%
  group_by(stock_group) %>%
  summarise(correlation = cor(abundance_index, sum_esc)) %>%
  mutate(label = paste0("r = ", round(correlation, 2)))

png(
  here::here("figs", "run_size", "escapement_aabm_corr.png"),
  height = 7, width = 7.5, units = "in", res = 250
)
ggplot(data = sg_esc2, aes(x = abundance_index, y = sum_esc / 1000)) +
  geom_point() +
  facet_wrap(~stock_group, scales = "free_y") +
  geom_text(data = cor_df, aes(x = -Inf, y = Inf, label = label),
            hjust = -0.1, vjust = 1.5, inherit.aes = FALSE) + 
  labs(x = "WCVI AABM Index", y = "Terminal Run Size") +
  ggsidekick::theme_sleek()
dev.off()

