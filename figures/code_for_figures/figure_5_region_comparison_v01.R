# 16 March 2023
# Figure 5
# Visualize differences in number-of-groups and group-size between regions

library(here)
library(tidyverse)
library(MetBrewer)
library(MCMCvis)

setwd(here::here("results"))
load("herbivore_case_study_results_v01.RData")

sp_key <- tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's", "hartebeest",
                "impala",  "Thomson's", "topi",
                "warthog", "waterbuck" ) )

MCMCvis::MCMCsummary(out, params = c("alpha0_contrast", "beta0_contrast")) %>% 
  as_tibble(rownames = "param") %>% 
  separate( param, into = c("param", "sp"), sep = "_") %>% 
  dplyr::mutate(sp = readr::parse_number(sp)) %>% 
  mutate( type = ifelse(grepl("alpha", param), "Number of groups", "Group size")) %>% 
  mutate(type = factor(type, levels = c("Number of groups", "Group size"))) %>% 
  dplyr::full_join( sp_key ) %>% 
  dplyr::select( sp_name, type, mean, l95 = `2.5%`, u95 = `97.5%`) %>% 
  dplyr::group_by( type ) %>% 
  dplyr::arrange( mean ) %>% 
  dplyr::mutate( sp_name = factor(sp_name, levels = unique(sp_name ) ) ) %>% 
  
  ggplot( aes( x = mean, y = sp_name ) ) +
  facet_wrap( ~type ) + 
  geom_vline(xintercept = 0,
             color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
             linetype = "dashed") + 
  geom_errorbar(aes(xmin = l95, xmax = u95),
                width = 0, 
                size = 1,
                color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  geom_point(size = 3,
             color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  theme_minimal() +
  labs(x = "Intercept contrast ( Talek - Mara )") + 
  theme(axis.title.y = element_blank(),
        axis.text = element_text(color = "black",
                                 size = 8),
        axis.title.x = element_text(color = "black", 
                                    size = 9),
        strip.text = element_text(color = "black", 
                                  size = 9),
        axis.line = element_line(color = "black",
                                 size = 0.1),
        plot.background = element_rect(color = NA, 
                                       fill = "white"), 
        panel.background = element_rect(color = NA, 
                                        fill = "white"))
  
setwd(here::here("figures"))
ggsave("figure_05_contrasts_v01.png",
       width = 5, 
       height = 2, 
       units = "in", 
       dpi = 300)