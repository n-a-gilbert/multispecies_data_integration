# 16 March 2023
# Figure 5
# Visualize differences in number-of-groups and group-size between regions

library(here)
library(tidyverse)
library(MetBrewer)
library(MCMCvis)
library(ggdist)
library(patchwork)

setwd(here::here("results"))
load("herbivore_case_study_results_v02.RData")

sp_key <- tibble::tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's", "hartebeest",
                "impala",  "Thomson's", "topi",
                "warthog", "waterbuck" ) )

community <- MCMCvis::MCMCsummary( out, params = c( "mu_alpha0_contrast", 
                                                    "mu_beta0_contrast" ) ) %>% 
  tibble::as_tibble( rownames = "param" ) %>% 
  dplyr::mutate( type = base::ifelse( base::grepl("alpha", param),
                                      "Mean: number of groups", "Mean: group size")) %>% 
  dplyr::mutate( type = base::factor( type, levels = c("Mean: number of groups", "Mean: group size"))) %>% 
  dplyr::select( type, mean, l95 = `2.5%`, u95 = `97.5%` )

species <- MCMCvis::MCMCsummary(out, params = c("alpha0_contrast", "beta0_contrast")) %>% 
  tibble::as_tibble( rownames = "param" ) %>% 
  tidyr::separate( param, into = c("param", "sp"), sep = "_") %>% 
  dplyr::mutate(sp = readr::parse_number(sp)) %>% 
  dplyr::mutate( type = base::ifelse(grepl("alpha", param), "Mean: number of groups", "Mean: group size")) %>% 
  dplyr::mutate(type = base::factor(type, levels = c("Mean: number of groups", "Mean: group size"))) %>% 
  dplyr::full_join( sp_key ) %>% 
  dplyr::select( sp_name, type, mean, l95 = `2.5%`, u95 = `97.5%`) %>% 
  dplyr::group_by( type ) %>% 
  dplyr::arrange( mean ) %>% 
  dplyr::mutate( sp_name = factor(sp_name, levels = unique(sp_name ) ) ) 

# plot differences in average number-of-groups and group size between region, both species- and community-level
( average_plot <-
    ggplot() +
    geom_rect( data = community, 
               aes( xmin = l95, 
                    xmax = u95), 
               ymin = -Inf, 
               ymax = Inf,
               fill = MetBrewer::MetPalettes$Hiroshige[[1]][7],
               color = NA, 
               alpha = 0.3 ) +
    geom_vline( data = community, 
                aes( xintercept = mean ), 
                color =  MetBrewer::MetPalettes$Hiroshige[[1]][7],
                size = 1,
                alpha = 0.3 ) +
    geom_vline( xintercept = 0,
                color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                linetype = "dashed" ) +
    facet_wrap( ~type ) +
    geom_errorbar( data = species,
                   aes( y = sp_name, 
                        xmin = l95, 
                        xmax = u95 ), 
                   width = 0, 
                   size = 1, 
                   color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    geom_point( data = species, 
                aes( y = sp_name, 
                     x = mean ), 
                size = 3, 
                color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    theme_minimal() +
    scale_x_continuous( limits = c(-3.01, 3.53)) +
    labs(title = "(a)") +
    theme(axis.title = element_blank(),
          axis.text.y = element_text(color = "black", size = 10),
          axis.text.x = element_blank(),
          plot.title = element_text(color = "black", size = 11),
          strip.text = element_text(color = "black", size = 11),
          axis.line.y = element_line(color = "black", size = 0.1),
          axis.line.x = element_blank(),
          plot.background = element_rect(color = NA, fill = "white"), 
          panel.background = element_rect(color = NA, fill = "white"),
          panel.grid = element_blank() ) )

# create a ggplot showing region differences in among-species variation in number-of-groups and group size
# variation_plot <-
#   MCMCvis::MCMCpstr( out, params = c("sd_alpha0_contrast"), type = "chains" )[[1]] %>% 
#   tibble::as_tibble() %>% 
#   tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "value") %>% 
#   tibble::add_column( type = "Interspecific variation: number of groups" ) %>% 
#   dplyr::full_join(
#     MCMCvis::MCMCpstr( out, params = c("sd_beta0_contrast"), type = "chains" )[[1]] %>% 
#       tibble::as_tibble() %>% 
#       tidyr::pivot_longer( 1:3000, names_to = "iter", values_to = "value" ) %>% 
#       tibble::add_column( type = "Interspecific variation: group size" )
#   ) %>%
#   dplyr::mutate(type = factor(type, levels = c("Interspecific variation: number of groups",
#                                                "Interspecific variation: group size" ))) %>%
#   ggplot( aes( x = value ) ) +
#   facet_wrap(~type) +
#   ggdist::stat_halfeye(
#     point_interval = "mean_qi",
#     adjust = 2,
#     .width = c(0.95),
#     fill = MetBrewer::MetPalettes$Hiroshige[[1]][7],
#     slab_alpha = 0.3, 
#     color = MetBrewer::MetPalettes$Hiroshige[[1]][7] ) +
#   geom_vline( xintercept = 0,
#               color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
#               linetype = "dashed" ) +
#   theme_minimal() +
#   scale_x_continuous( limits = c(-3.01, 3.53), 
#                       breaks = c( -3, -1.5, 0, 1.5, 3)) +
#   labs( x = "Contrast ( Talek - Mara )",
#         title = "(b)" ) + 
#   theme( panel.grid = element_blank(),
#          axis.title.y = element_blank(),
#          axis.line.y = element_blank(),
#          axis.text.y = element_blank(), 
#          plot.title = element_text(color = "black", size = 9),
#          axis.text.x = element_text(color = "black", size = 8),
#          axis.title.x = element_text(color = "black", size = 9),
#          strip.text = element_text(color = "black", size = 9),
#          axis.line = element_line(color = "black", size = 0.1),
#          plot.background = element_rect(color = NA, fill = "white"), 
#          panel.background = element_rect(color = NA, fill = "white")) +
#   ylim( c(-0.5, 1))


( variation_plot <- MCMCvis::MCMCpstr( out, params = c("sd_alpha0_contrast"), type = "chains" )[[1]] %>% 
    tibble::as_tibble() %>% 
    tidyr::pivot_longer(1:3000, names_to = "iter", values_to = "value") %>% 
    tibble::add_column( type = "SD: number of groups" ) %>% 
    dplyr::full_join(
      MCMCvis::MCMCpstr( out, params = c("sd_beta0_contrast"), type = "chains" )[[1]] %>% 
        tibble::as_tibble() %>% 
        tidyr::pivot_longer( 1:3000, names_to = "iter", values_to = "value" ) %>% 
        tibble::add_column( type = "SD: group size" )
    ) %>%
    dplyr::mutate(type = factor(type, levels = c("SD: number of groups",
                                                 "SD: group size" ))) %>% 
    
    group_by(type) %>% 
    summarise( mean = mean(value), 
               l95 = quantile(value, c(0.025)), 
               u95 = quantile(value, c(0.975))) %>% 
    
    ggplot( aes( x = mean ) ) +
    facet_wrap(~type) +
    geom_rect( aes( xmin = l95, 
                    xmax = u95), 
               ymin = -Inf, 
               ymax = Inf,
               fill = MetBrewer::MetPalettes$Hiroshige[[1]][7],
               color = NA, 
               alpha = 0.3 ) +
    geom_vline( aes( xintercept = mean ), 
                color =  MetBrewer::MetPalettes$Hiroshige[[1]][7],
                size = 1,
                alpha = 0.3 ) +
    geom_vline( xintercept = 0,
                color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                linetype = "dashed" ) +
    theme_minimal() +
    scale_x_continuous( limits = c(-3.01, 3.53), 
                        breaks = c( -3, -1.5, 0, 1.5, 3)) +
    labs( x = "Contrast ( Talek - Mara )",
          title = "(b)" ) + 
    theme( panel.grid = element_blank(),
           axis.title.y = element_blank(),
           axis.line.y = element_blank(),
           axis.text.y = element_blank(), 
           plot.title = element_text(color = "black", size = 11),
           axis.text.x = element_text(color = "black", size = 10),
           axis.title.x = element_text(color = "black", size = 11),
           strip.text = element_text(color = "black", size = 11),
           axis.line = element_line(color = "black", size = 0.1),
           plot.background = element_rect(color = NA, fill = "white"), 
           panel.background = element_rect(color = NA, fill = "white")) +
    ylim( c(-0.5, 1)) )

# combine with patchwork 
average_plot + variation_plot + plot_layout( nrow = 2, heights = c(10, 1))

setwd(here::here("figures"))
ggsave("figure_05_contrasts_v01.png",
       width = 4.75, 
       height = 3.5, 
       units = "in", 
       dpi = 300)

