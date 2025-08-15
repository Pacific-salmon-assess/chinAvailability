## Salmonid Enhancement Program Figures
# Estimates of PBT rates and pHOS using data provided by Eric Rondeau or Brock
# Ramshaw
# Oct 9, 2024

# Set French language option
FRENCH <- FALSE

# Create appropriate figure directories
if (FRENCH) {
  dir.create("figs-french", showWarnings = FALSE)
  dir.create("figs-french/pbt_phos", showWarnings = FALSE)
  fig_dir <- "figs-french"
} else {
  dir.create("figs/pbt_phos", showWarnings = FALSE)
  fig_dir <- "figs"
}

# Translation helper function
tr <- function(english, french) {
  if (FRENCH) french else english
}

# Helper function for figure paths
fig_path <- function(filename) {
  file.path(fig_dir, filename)
}

library(tidyverse)


# pbt_rate <- readRDS(here::here("data", "sep", "cleaned_pbt.rds"))
long_rec <- readRDS(here::here("data", "rec", "rec_gsi.rds"))


# import SEP release estimates to identify major stocks (defined as facilities
# w/ greater than 140k average releases per year, following Brock's advice)
# sep_rel <- readxl::read_xlsx(
#   here::here("data", "sep", "2018-2022BY avg release.xlsx")
# ) %>% 
#   janitor::clean_names() %>% 
#   filter(
#     avg_rel > 140000 | stock_name %in% c("Nechako R", "Chilko R")
#   ) %>% 
#   mutate(
#     #rename to match pbt_rate df
#     stock_name = case_when(
#       stock_name == "Big Qualicum R" ~ "Qualicum River",
#       stock_name == "L Qualicum R" ~ "Little Qualicum River",
#       stock_name == "Shuswap R Low" ~ "Shuswap River Lower",
#       stock_name == "Shuswap R Middle" ~ "Shuswap River Middle",
#       stock_name == "Chehalis R" ~ "Chehalis River Summer",
#       stock_name == "Kitsumkalum R" ~ "Kitsumkalum River Lower",
#       stock_name == "Chilliwack R" ~ "Chilliwack River Fall",
#       stock_name == "Nanaimo R" ~ "Nanaimo River Fall",
#       grepl(" R", stock_name) ~ paste(stock_name, "iver", sep = ""),
#       grepl(" Cr", stock_name) ~ paste(stock_name, "eek", sep = ""),
#       TRUE ~ stock_name
#     ) %>% 
#       toupper(),
#     pbt_stock = gsub(" ", "_", stock_name)
#   )


# pbt_rate_plot <- pbt_rate %>% 
#   mutate(
#     brood_year = as.factor(brood_year),
#     stock_group = case_when(
#       pbt_stock %in% 
#         c("CAPILANO_RIVER", "COWICHAN_RIVER", "LITTLE_QUALICUM_RIVER",
#           "NANAIMO_RIVER_FALL", "QUALICUM_RIVER", "PUNTLEDGE_RIVER") ~ "ECVI_SOMN",
#       pbt_stock %in% c("CHEHALIS_RIVER_SUMMER", "CHILLIWACK_RIVER_FALL",
#                        "HARRISON_RIVER") ~ "FR_Fall",
#       pbt_stock %in% c("SHUSWAP_RIVER_LOWER", "SHUSWAP_RIVER_MIDDLE") ~ 
#         "FR_Sum_4sub1",
#       pbt_stock %in% c("CHILKO_RIVER", "NECHAKO_RIVER") ~ "FR_Sum_5sub2",
#       pbt_stock == "NICOLA_RIVER" ~ "FR_Spr_4sub2",
#       TRUE ~ "WCVI"
#     ) %>% 
#       factor(., levels = levels(long_rec$stock_group)),
#     pbt_stock = fct_reorder(pbt_stock, as.numeric(stock_group))
#   ) %>% 
#   # focus on large facilities
#   filter(
#     pbt_stock %in% sep_rel$pbt_stock,
#     # remove northern stocks 
#     !pbt_stock %in% c("ATNARKO_RIVER", "BURMAN_RIVER", 
#                       "GOLD_RIVER", "KITIMAT_RIVER", "KITSUMKALUM_RIVER_LOWER",
#                       "LANG_CREEK", "LEINER_RIVER", "MARBLE_RIVER", 
#                       "NIMPKISH_RIVER", "QUINSAM_RIVER", "TAHSIS_RIVER", 
#                       "WANNOCK_RIVER", "WOSS_RIVER", "YAKOUN_RIVER")
#   )
# saveRDS(
#   pbt_rate_plot,
#   here::here(
#     "data", "sep_plot_data.rds"
#   )
# )
pbt_rate_plot <- readRDS(
  here::here(
    "data", "sep_plot_data.rds"
  )
)

smu_colour_pal <- c("grey30", "#3182bd", "#bdd7e7", "#bae4bc", "#6A51A3",
                    "#CBC9E2", "#67000D", "#A50F15", "#EF3B2C", "#FC9272", 
                    "#FCBBA1")
names(smu_colour_pal) <- levels(long_rec$stock_group)


pbt_coverage <- ggplot(pbt_rate_plot, aes(x = brood_year, y = tag_rate)) +
  geom_point(aes(fill = stock_group), shape = 21) +
  geom_hline(aes(yintercept = 0.8), colour = "red") +
  facet_wrap(~pbt_stock, ncol = 3) +
  labs(y = tr("Proportion of Brood Stock with PBT", "Proportion du stock reproducteur avec PBT"), 
       x = tr("Brood Year", "Année de fraie"), 
       fill = tr("Stock", "Stock")) +
  scale_x_discrete(
    breaks = seq(2013, 2021, by = 2)
  ) +
  scale_fill_manual(values = smu_colour_pal, name = tr("Stock", "Stock")) +
  ggsidekick::theme_sleek() +
  theme(legend.position = "top")


png(
  fig_path("pbt_phos/pbt_rate.png"),
  height = 6.5, width = 5.75, units = "in", res = 250
)
pbt_coverage
dev.off()



# use SEP data to fit simple LM to estimate mean hatchery contribution for WCVI 
# and ECVI systems
# phos_plot_dat <- readxl::read_xlsx(
#   here::here("data", "sep", "2024-08-20 annual pHOS summary.xlsx")
# ) %>% 
#   filter(
#     estimate_type %in% c("Direct", "Partial Direct"),
#     cu_acronym %in% c("SWVI", "Chil_transp_F", "QP-fall", "MFR-summer",
#                       "CWCH-KOK", "LFR-fall", "EVI-fall", "EVIGStr-sum",
#                       "LTh", "NEVI", "STh-1.3", "STh-SHUR", "NEVI")
#   ) %>% 
#   mutate(
#     stock_group = case_when(
#       cu_acronym == "SWVI" ~ "WCVI",
#       cu_acronym %in% c("Chil_transp_F", "LFR-fall") ~ "FR_Fall",
#       grepl("EVI", cu_acronym) | cu_acronym %in% c("CWCH-KOK", "QP-fall") ~ 
#         "ECVI_SOMN",
#       cu_acronym == "MFR-summer" ~ "FR_Sum_5sub2",
#       cu_acronym == "LTh" ~ "FR_Spr_4sub2",
#       cu_acronym == "STh-1.3" ~ "FR_Spr_5sub2",
#       cu_acronym == "STh-SHUR" ~ "FR_Sum_4sub1"
#     ),
#     stock_group = factor(stock_group, levels = levels(long_rec$stock_group)),
#     population = gsub(" River", "", population),
#     population = gsub(" Creek", "", population),
#     population = fct_reorder(population, as.numeric(stock_group)),
#     pHOS = ifelse(pHOS == "1", 0.999, pHOS),
#     pHOS = ifelse(pHOS == "0", 0.001, pHOS)
#   ) 
# saveRDS(
#   phos_plot_dat,
#   here::here(
#     "data", "phos_plot_data.rds"
#   )
# )
phos_plot_dat <- readRDS(
  here::here(
    "data", "phos_plot_data.rds"
  )
)

phos_box <- ggplot(phos_plot_dat) +
  geom_boxplot(aes(x = population, y = pHOS, fill = stock_group)) +
  scale_fill_manual(values = smu_colour_pal, name = tr("Stock", "Stock")) +
  labs(y = tr("Proportion Hatchery Origin Spawners", "Proportion de géniteurs d'origine d'écloserie"), 
       fill = tr("Stock", "Stock")) +
  ggsidekick::theme_sleek() +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "top")


library(glmmTMB)
# 
# # only one Fraser Fall population
fit <- glmmTMB(
  pHOS ~ 0 + stock_group + (1 | population), data = phos,
  family = beta_family()
)

# # intercept is mean after accounting for interannual and among population 
# # variation
dd <- summary(fit)
ci <-  confint(fit, parm = dd$coefficients$cond %>% rownames(), level = 0.95)
boot::inv.logit(ci)

png(
  fig_path("pbt_phos/phos_box.png"),
  height = 5, width = 8.25, units = "in", res = 250
)
phos_box
dev.off()
