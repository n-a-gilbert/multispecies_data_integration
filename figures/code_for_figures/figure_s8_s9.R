# 16 March 2023
# code to visualize estimates of abundance underlying distance sampling and count data

library(here)
library(tidyverse)
library(MCMCvis)
library(MetBrewer)
library(patchwork)

setwd(here::here("data"))

load("distance_sampling_data_v01.RData")
load("count_data_v01.RData")

ng_data <- final2 |>
  dplyr::select( sp, site, rep, ng, gs, area, region ) |> 
  dplyr::filter( !base::is.na( ng ) ) |> 
  dplyr::group_by( sp, site, rep, area, region ) |> 
  dplyr::summarise( ng = base::unique( ng ),
                    yN_DS = base::sum(gs, na.rm = TRUE ) ) |> 
  dplyr::arrange( sp, site, rep )

setwd(here::here("results"))
load("herbivore_case_study_results_v01.RData")

sp_key <- tibble::tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's", "hartebeest",
                "impala",  "Thomson's", "topi",
                "warthog", "waterbuck" ) )

N_DS <- MCMCvis::MCMCsummary(out, params = "N_DS") |> 
  tibble::as_tibble(rownames = "param") |> 
  tibble::add_column( sp = constants$SP_NG , 
                      site = ng_data$site , 
                      rep = ng_data$rep , 
                      region = constants$REGION_NG ,
                      area = ng_data$area ) |> 
  dplyr::full_join( sp_key ) |> 
  dplyr::select( sp, sp_name, site, area, rep, region, mean, sd ) |> 
  dplyr::mutate( density_mean = mean / area,  # convert to density to aid comparison with count data
                 density_sd = sd / area )

( panela <- ggplot2::ggplot( N_DS, 
                           aes( y = site, x = rep, fill = log1p(density_mean))) +
  ggplot2::facet_wrap( ~sp_name ) +
  ggplot2::geom_tile(color = NA) + 
  ggplot2::scale_fill_gradientn(colors = met.brewer(name="Hiroshige", n=100,
                                                    type="continuous",
                                                    direction = -1),
                                limits = c(0, 6.83)) +
  ggplot2::geom_hline(yintercept = 13.5, color = "black", size = 0.7) + 
  ggplot2::theme_minimal() +
  ggplot2::scale_y_continuous(breaks = c(1, 6, 12, 18))+
  ggplot2::labs(x = "Temporal replicate", 
                y = "Distance sampling site",
                fill = "Log1P ( Animals / km2 )",
                title = "Distance sampling")+
  ggplot2::theme(panel.grid = element_blank(),
                 axis.line = element_line(color = "black", 
                                          size = 0.1),
                 axis.title = element_text(color = "black", 
                                           size = 9),
                 axis.text = element_text(color = "black", 
                                          size = 8),
                 plot.title = element_text(hjust = 0.5, 
                                           color = "black", 
                                           size = 10),
                 strip.text = element_text(color = "black", 
                                           size = 8),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 8, color = "black")) +
  ggplot2::guides(fill = guide_colorbar(title.hjust = 0.75,
                                        ticks = FALSE,
                                        barheight = 0.5)) )

N_C <- MCMCvis::MCMCsummary(out, params = "N_TC") |> 
  tibble::as_tibble(rownames = "param") |> 
  tibble::add_column( sp = constants$SP_TC , 
                      site = transect_data$site , 
                      rep = transect_data$rep , 
                      region = constants$REGION_TC ,
                      area = transect_data$area ) |> 
  dplyr::full_join( sp_key ) |> 
  dplyr::select( sp, sp_name, site, area, rep, region, mean, sd ) |> 
  dplyr::mutate( density_mean = mean / area, 
                 density_sd = sd / area )

( panelb <- ggplot2::ggplot( N_C, 
                           aes( y = site, x = rep, fill = log1p(density_mean))) +
  ggplot2::facet_wrap( ~sp_name ) +
  ggplot2::geom_tile(color = NA) + 
  ggplot2::scale_fill_gradientn(colors = met.brewer(name="Hiroshige", n=100,
                                                    type="continuous",
                                                    direction = -1),
                                limits = c(0, 6.83)) +
  ggplot2::geom_hline(yintercept = 6.5, color = "black", size = 0.7) + 
  ggplot2::theme_minimal() +
  ggplot2::scale_y_continuous(breaks = c(1, 4, 8))+
  ggplot2::labs(x = "Temporal replicate", 
                y = "Count site",
                fill = "Log1P ( Animals / km2 )",
                title = "Counts")+
  ggplot2::theme(panel.grid = element_blank(),
                 axis.line = element_line(color = "black", 
                                          size = 0.1),
                 axis.title = element_text(color = "black", 
                                           size = 9),
                 axis.text = element_text(color = "black", 
                                          size = 8),
                 plot.title = element_text(hjust = 0.5, 
                                           color = "black", 
                                           size = 10),
                 strip.text = element_text(color = "black", 
                                           size = 8),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 8, color = "black")) +
  ggplot2::guides(fill = guide_colorbar(title.hjust = 0.75,
                                        ticks = FALSE,
                                        barheight = 0.5)) )

# posterior means
mean_density <- panela + panelb + plot_layout(guides = "collect") &
  theme(legend.position='bottom',
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)) 

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s8.png", 
  mean_density, 
  width = 7, 
  height = 4, 
  units = "in", 
  dpi = 300
)

panela_sd <- ggplot2::ggplot( N_DS, 
                              aes( y = site, x = rep, fill = log1p(density_sd))) +
  ggplot2::facet_wrap( ~sp_name ) +
  ggplot2::geom_tile(color = NA) + 
  ggplot2::scale_fill_gradientn(colors = met.brewer(name="Hiroshige", n=100,
                                                    type="continuous",
                                                    direction = -1),
                                limits = c(0, 4.56)) +
  ggplot2::geom_hline(yintercept = 13.5, color = "black", size = 0.7) + 
  ggplot2::theme_minimal() +
  ggplot2::scale_y_continuous(breaks = c(1, 6, 12, 18))+
  ggplot2::labs(x = "Temporal replicate", 
                y = "Distance sampling site",
                fill = "Log1P ( Animals / km2 )",
                title = "Distance sampling")+
  ggplot2::theme(panel.grid = element_blank(),
                 axis.line = element_line(color = "black", 
                                          size = 0.1),
                 axis.title = element_text(color = "black", 
                                           size = 9),
                 axis.text = element_text(color = "black", 
                                          size = 8),
                 plot.title = element_text(hjust = 0.5, 
                                           color = "black", 
                                           size = 10),
                 strip.text = element_text(color = "black", 
                                           size = 8),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 8, color = "black")) +
  ggplot2::guides(fill = guide_colorbar(title.hjust = 0.75,
                                        ticks = FALSE,
                                        barheight = 0.5))

panelb_sd <- ggplot2::ggplot( N_C, 
                              aes( y = site, x = rep, fill = log1p(density_sd))) +
  ggplot2::facet_wrap( ~sp_name ) +
  ggplot2::geom_tile(color = NA) + 
  ggplot2::scale_fill_gradientn(colors = met.brewer(name="Hiroshige", n=100,
                                                    type="continuous",
                                                    direction = -1),
                                limits = c(0, 4.56)) +
  ggplot2::geom_hline(yintercept = 6.5, color = "black", size = 0.7) + 
  ggplot2::theme_minimal() +
  ggplot2::scale_y_continuous(breaks = c(1, 4, 8))+
  ggplot2::labs(x = "Temporal replicate", 
                y = "Count site",
                fill = "Log1P ( Animals / km2 )",
                title = "Counts")+
  ggplot2::theme(panel.grid = element_blank(),
                 axis.line = element_line(color = "black", 
                                          size = 0.1),
                 axis.title = element_text(color = "black", 
                                           size = 9),
                 axis.text = element_text(color = "black", 
                                          size = 8),
                 plot.title = element_text(hjust = 0.5, 
                                           color = "black", 
                                           size = 10),
                 strip.text = element_text(color = "black", 
                                           size = 8),
                 legend.title = element_text(size = 9, color = "black"),
                 legend.text = element_text(size = 8, color = "black")) +
  ggplot2::guides(fill = guide_colorbar(title.hjust = 0.75,
                                        ticks = FALSE,
                                        barheight = 0.5))

#posterior SDs
sd_density <- panela_sd + panelb_sd + plot_layout(guides = "collect") &
  theme(legend.position='bottom',
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)) 

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s9.png", 
  sd_density, 
  width = 7, 
  height = 4, 
  units = "in", 
  dpi = 300)
