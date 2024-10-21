## Aggregate Abundance
# Figures of escapement/terminal return for relevant stock groups
# Oct 19, 2024


library(tidyverse)


esc_raw_dat <- read.csv(
  here::here(
    "data", "run_size", "run_size_by_stock_group.csv"
  ),
  skip = 1
) %>% 
  janitor::clean_names() %>% 
  pivot_longer(cols = -year, names_to = "stock", values_to = "esc")
  
  
  
# only mean values available for early years; expand accordingly
esc_dat <- expand.grid(
  year_n = seq(1981, 2023, by = 1),
  stock = unique(esc_raw_dat$stock)
) %>% 
  mutate(
    year = case_when(
      year_n < 1986 ~ "1981-1985",
      year_n < 1991 ~ "1986-1990",
      year_n < 1996 ~ "1991-1995",
      year_n < 2001 ~ "1996-2000",
      year_n < 2006 ~ "2001-2005",
      year_n < 2011 ~ "2006-2010",
      TRUE ~ as.character(year_n)
    )
  ) %>% 
  left_join(
    ., esc_raw_dat, by = c("year", "stock")
  ) %>% 
  mutate(
    stock_group = case_when(
      grepl("spring_5_2", stock) ~ "Fraser_Spring_5.2",
      grepl("spring_4_2", stock) ~ "Fraser_Spring_4.2",
      grepl("summer_5_2", stock) ~ "Fraser_Summer_5.2",
      grepl("summer_4_1", stock) ~ "Fraser_Summer_4.1",
      grepl("fr_fall", stock) ~ "Fraser_Fall",
      stock %in% c("lower_gst_a17_19_28_29", "upper_gst_a13_16", "jst_a12") ~ 
        "ECVI_SOMN",
      stock %in% c("nwvi_fall_a25_27", "swvi_fall_a20_24") ~ "WCVI",
      stock == "columbia_river_spring" ~ "Col_Spring",
      grepl("cr_", stock) ~ "Col_Summer_Fall",
      stock %in% c("willapa", "grays_h_spring", "grays_h_fall", "queets", 
                   "hoh_s_s", "hoh_fall", "quillayute_s_s", "quillayute_fall",
                   "hoko_s_f", "juan_de_fuca") ~ "WA_Coastal",
      stock %in% c("nooksack_samish", "skagit", "hood_canal", 
                   "stillaquamish_snohomish", "south_puget_sound") ~ "PSD"
    ) %>% 
      factor(
        .,
        levels = c("WA_Coastal", "Col_Spring", "Col_Summer_Fall", "PSD",  
                   "WCVI", "ECVI_SOMN", "Fraser_Spring_4.2",
                   "Fraser_Spring_5.2", "Fraser_Summer_5.2", "Fraser_Summer_4.1",
                   "Fraser_Fall"),
        labels = c(
          "Washington\nCoastal", "Columbia\nSpring", "Columbia\nSummer/Fall",  
          "Puget Sound", "Westcoast\nVan. Island", "Eastcoast\nVan. Island",
          "Fraser\nSpring 4.2", "Fraser\nSpring 5.2", "Fraser\nSummer 5.2", 
          "Fraser\nSummer 4.1", "Fraser \nFall"
        )
      ),
    region = ifelse(
      stock_group %in% c(
        "Puget Sound",
        "Fraser\nSpring 4.2", "Fraser\nSpring 5.2", "Fraser\nSummer 5.2", 
        "Fraser\nSummer 4.1", "Fraser \nFall"
        ) |
        stock %in% c("lower_gst_a17_19_28_29", "upper_gst_a13_16", 
                     "juan_de_fuca"),
      "salish_sea",
      "coastal"
    ) %>% 
      factor(
        ., labels = c("Coastal", "Salish Sea")
      ),
    average = ifelse(
      year_n < 2011, TRUE, FALSE
    )
  )

sg_esc <- esc_dat %>% 
  group_by(year_n, average, stock_group) %>% 
  summarize(
    sum_esc = sum(esc, na.rm = T)
  ) # %>% 
  # remove incomplete data (southern US not available)
  # filter(
  #   !is.na(sum_esc)
  # )

re_esc <- esc_dat %>%
  group_by(year_n, average, region) %>% 
  summarize(
    sum_esc = sum(esc, na.rm = TRUE)
  ) %>% 
  # remove incomplete data (southern US not available)
  filter(
    !year_n == "2023"
  )


## FIGURES ---------------------------------------------------------------------

smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(esc_dat$stock_group)
average_pal <- c(0.5, 1.0)
names(average_pal) <- c(TRUE, FALSE)

sg_plot <- ggplot(sg_esc) +
  geom_line(aes(x = year_n, y = sum_esc / 1000, colour = stock_group, 
                lty = average)) +
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
  geom_bar(aes(x = year_n, y = sum_esc / 1000, fill = region, alpha = average),
           stat = "identity") +
  scale_alpha_manual(values = average_pal) +
  ggsidekick::theme_sleek() +
  scale_fill_viridis_d() +
  labs(y = "Terminal Abundance (thousands)") +
  theme(
    legend.position = "top",
    axis.title.x = element_blank()
  ) +
  labs(fill = NULL) +
  guides(
    alpha = "none"
  )

png(
  here::here("figs", "run_size", "escapement_region.png"),
  height = 5, width = 5.5, units = "in", res = 250
)
reg_plot
dev.off()
