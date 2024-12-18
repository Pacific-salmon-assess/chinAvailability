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
        "Puget Sound", "East Coast\nVan. Island",
        "Fraser\nSpring 4.2", "Fraser\nSpring 5.2", "Fraser\nSummer 5.2", 
        "Fraser\nSummer 4.1", "Fraser \nFall"
      ),
      "salish_sea",
      "coastal"
    ) %>% 
      factor(
        ., labels = c("Outer Coast", "Salish Sea")
      )
    # obs = ifelse(is.na(esc), 0, 1)
  ) %>%
  group_by(stock) %>% 
  filter(
    # remove stocks that are missing esc data after 1983 or 2019
    !any((is.na(esc) & year == "1983") | (is.na(esc) & year == "2019"))#,
    # remove years with gappy data
    # !(year < 1983 | year == "2023")
  ) %>% 
  ungroup()


# # visualize data presence
# ggplot(
#   esc_raw_dat,
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


reg_plot <- ggplot(re_esc %>% filter(year > 1982 & year < 2023)) +
  geom_bar(aes(x = year, y = sum_esc / 1000, fill = region),
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
