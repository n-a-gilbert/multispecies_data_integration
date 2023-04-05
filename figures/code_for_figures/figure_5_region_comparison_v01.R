# 16 March 2023
# Figure 5
# Visualize differences in number-of-groups and group-size between regions

library(here)
library(tidyverse)
library(MetBrewer)
library(MCMCvis)
library(ggdist)
library(patchwork)
library(reshape)

setwd(here::here("results"))
load("herbivore_case_study_results_v03.RData")

sp_key <- tibble::tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's", "hartebeest",
                "impala",  "Thomson's", "topi",
                "warthog", "waterbuck" ) )

sp_n <- MCMCvis::MCMCpstr( out, params = "N_region", type = "chains")[[1]] %>% 
  reshape::melt(., varnames = c("sp", "region", "iter")) %>% 
  tibble::as_tibble() %>% 
  dplyr::mutate( value = log( value ) ) %>% 
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>% 
  dplyr::full_join(sp_key) %>% 
  dplyr::group_by( sp_name, region ) %>% 
  dplyr::summarise( mean = mean(value), 
                    l95 = quantile(value, c(0.025)), 
                    u95 = quantile(value, c(0.975))) %>% 
  dplyr::mutate(sp_name = factor(sp_name, 
                                 levels = c("Grant's", "hartebeest", "eland", "giraffe", "elephant", "waterbuck", 
                                            "warthog", "buffalo", "topi", "Thomson's", "impala")))

sp_sig <- MCMCvis::MCMCpstr( out, params = "N_region", type = "chains")[[1]] %>% 
  reshape::melt(., varnames = c("sp", "region", "iter")) %>% 
  tibble::as_tibble() %>%  
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  pivot_wider(names_from = region, values_from = value) %>% 
  mutate( diff = Talek - Mara) %>% 
  group_by(sp) %>% 
  summarise( mean = mean(diff), 
             l95 = quantile(diff, c(0.025)), 
             u95 = quantile(diff, c(0.975))) %>% 
  full_join(sp_key) %>% 
  mutate(region = ifelse( mean > 0, "Talek", "Mara")) %>% 
  mutate(sig = ifelse( l95 < 0 & u95 > 0, NA, "*")) %>% 
  dplyr::select(sp_name, sig, region) %>% 
  add_column( mean = -3.9)

com_n <- MCMCvis::MCMCpstr( out, 
                            params = c("mu_N"),
                            type = "chains")[[1]] %>%
  reshape::melt(., varnames = c("region", "iter")) %>% 
  as_tibble() %>% 
  mutate(region = parse_number(as.character(region))) %>% 
  mutate(region = ifelse(region == 1, "Mara", "Talek")) %>% 
  mutate(value = log(value)) %>% 
  group_by(region) %>% 
  summarise( mean = mean(value), 
             l95 = quantile(value, c(0.025)), 
             u95 = quantile(value, c(0.975)))

( nplot <- ggplot() +
    geom_rect(
      data = com_n, 
      aes( xmin = l95, 
           xmax = u95,
           fill = region), 
      ymin = -Inf, 
      ymax = Inf,
      alpha = 0.1) +
    geom_vline(
      data = com_n, 
      aes(xintercept = mean, 
          color = region), 
      size = 1, 
      alpha = 0.2) +
    geom_errorbar(
      data = sp_n, 
      aes(y = sp_name, 
          xmin = l95, 
          xmax = u95, 
          color = region),
      width = 0,
      size = 1,
      position = position_dodge(width = 0.5)) +
    geom_point(
      data = sp_n, 
      aes(y = sp_name,
          x = mean, 
          color = region),
      size = 3,
      position = position_dodge(width = 0.5)) +
    geom_text( data = sp_sig, 
               aes(x = mean, 
                   y = sp_name, 
                   label = sig, 
                   color = region),
               size = 6,
               vjust = 0.75,
               fontface = "bold",
               show.legend = FALSE) +
    scale_fill_manual("Region",
                      values = MetBrewer::MetPalettes$Greek[[1]][c(2, 5)]) +
    scale_color_manual( "Region",
                        values = MetBrewer::MetPalettes$Greek[[1]][c(2, 5)]) +
    scale_x_continuous( limits = c(-4, 5.75), 
                        expand = c(0.01, 0)) +
    theme_minimal() +
    labs( x = " Log ( Abundance ) ",
          title = "(a)") + 
    theme( panel.grid = element_blank(),
           axis.title.y = element_blank(),
           legend.text = element_text(color = "black", size = 10), 
           legend.title = element_text(color = "black", size = 11),
           legend.position = c(0.8, 0.1),
           plot.title = element_text(color = "black", size = 10),
           axis.text = element_text(color = "black", size = 10),
           axis.title.x = element_text(color = "black", size = 11),
           axis.line = element_line(color = "black", size = 0.1),
           plot.background = element_rect(color = NA, fill = "white"), 
           panel.background = element_rect(color = NA, fill = "white")) +
    guides(color = guide_legend(reverse = TRUE), 
           fill = guide_legend(reverse = TRUE)) )

community <- MCMCvis::MCMCsummary( out, params = c( "mu_alpha0_contrast",
                                                    "mu_beta0_contrast" ) ) %>%
  tibble::as_tibble( rownames = "param" ) %>%
  dplyr::mutate( type = base::ifelse( base::grepl("alpha", param),
                                      "Number of groups", "Group size")) %>%
  dplyr::mutate( type = base::factor( type, levels = c("Number of groups", "Group size"))) %>%
  dplyr::select( type, mean, l95 = `2.5%`, u95 = `97.5%` )

species <- MCMCvis::MCMCsummary(out, params = c("alpha0_contrast", "beta0_contrast")) %>%
  tibble::as_tibble( rownames = "param" ) %>%
  tidyr::separate( param, into = c("param", "sp"), sep = "_") %>%
  dplyr::mutate(sp = readr::parse_number(sp)) %>%
  dplyr::mutate( type = base::ifelse(grepl("alpha", param), "Number of groups", "Group size")) %>%
  dplyr::mutate(type = base::factor(type, levels = c("Number of groups", "Group size"))) %>%
  dplyr::full_join( sp_key ) %>%
  dplyr::select( sp_name, type, mean, l95 = `2.5%`, u95 = `97.5%`) %>%
  dplyr::group_by( type ) %>%
  dplyr::arrange( mean ) %>%
  mutate(sp_name = factor(sp_name,
                          levels = c("Grant's", "hartebeest", "eland", "giraffe", "elephant", "waterbuck",
                                     "warthog", "buffalo", "topi", "Thomson's", "impala")))

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
    facet_wrap( ~type, ncol = 1 ) +
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
    scale_x_continuous( limits = c(-3.03, 3.54),
                        breaks = c( -3, -1.5, 0, 1.5, 3)) +
    labs(x = "Contrast ( Talek - Mara )",
         title = "(b)") +
    theme(axis.title.y = element_blank(),
          axis.title.x = element_text(size = 11, color = "black"),
          axis.text = element_text(color = "black", size = 10),
          plot.title = element_text(color = "black", size = 10),
          strip.text = element_text(color = "black", size = 11),
          axis.line = element_line(color = "black", size = 0.1),
          plot.background = element_rect(color = NA, fill = "white"),
          panel.background = element_rect(color = NA, fill = "white"),
          panel.grid = element_blank() ) )

nplot + average_plot + plot_layout(widths = c(1.5, 1))

setwd(here::here("figures"))
ggsave(
  "region_comparsion_v01.png", 
  width = 7,
  height = 5, 
  units = "in", 
  dpi = 300)