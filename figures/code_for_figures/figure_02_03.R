library(here)
library(tidyverse)
library(MetBrewer)
library(patchwork)

setwd(here::here("results"))
load("main_simulation_results_v01.RData")

prop_labs <- icm_n_ds_rb |> 
  tibble::add_column( source = "Distance sampling") |> 
  dplyr::full_join(icm_info) |> 
  dplyr::full_join( add_column( full_join(icm_n_tc_rb, icm_info), source = "Counts" ) ) |> 
  dplyr::mutate(tot = ifelse(source == "Counts", totTC, totDS)) |>
  dplyr::select( simrep, sp, mean, sd, l95, u95, contain0, tot, nsites_tc_fact, p_bias, source) |> 
  dplyr::mutate(p_bias = ifelse(p_bias == 0, "Same    \ndetectability", "Count    \ndetectability\n20% lower")) |> 
  dplyr::mutate(p_bias = factor(p_bias, levels = c(
    "Count    \ndetectability\n20% lower",
    "Same    \ndetectability"))) |> 
  dplyr::mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) |> 
  dplyr::mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
    "4x more count data",
    "Same amount of both data"))) |> 
  dplyr::mutate( contain0 = ifelse(contain0 == 0, "False", "True")) |> 
  dplyr::ungroup() |> 
  dplyr::group_by( source, nsites_tc_fact, p_bias ) |> 
  dplyr::mutate(n = n()) |> 
  dplyr::summarise( prop = round (sum( contain0 == "False", na.rm = TRUE) / n, 2 )) |> 
  dplyr::mutate(prop = paste0( (100 - (prop * 100)), "% unbiased")) |> 
  dplyr::distinct() |> 
  tibble::add_column( mean = c( 1.65, 1.65, 
                                1.65, 1.65, 
                                1.65, 1.65, 
                                1.65, 1.65))

( panela <- icm_n_ds_rb |> 
    tibble::add_column( source = "Distance sampling") |> 
    dplyr::full_join(icm_info) |> 
    dplyr::full_join( add_column( full_join(icm_n_tc_rb, icm_info), source = "Counts" ) ) |> 
    dplyr::mutate(tot = ifelse(source == "Counts", totTC, totDS)) |>
    dplyr::select( simrep, sp, mean, sd, l95, u95, contain0, tot, nsites_tc_fact, p_bias, source) |> 
    dplyr::mutate(p_bias = ifelse(p_bias == 0, "Same    \ndetectability", "Count    \ndetectability\n20% lower")) |> 
    dplyr::mutate(p_bias = factor(p_bias, levels = c(
      "Count    \ndetectability\n20% lower",
      "Same    \ndetectability"))) |> 
    dplyr::mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) |> 
    dplyr::mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
      "4x more count data",
      "Same amount of both data"))) |> 
    dplyr::mutate( contain0 = ifelse(contain0 == 0, "False", "True")) |> 
    
    ggplot2::ggplot( aes( x = mean, y = p_bias, fill = nsites_tc_fact )) +
    ggplot2::facet_wrap(~source, scales = "free_x") +
    ggplot2::geom_vline(xintercept = 0,
                        color = "gray60",
                        linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.alpha = 0.2,
                          outlier.size = 0.75, 
                          size = 0.25) +
    ggplot2::geom_label(data = prop_labs,
                        aes(label = prop,
                            color = nsites_tc_fact),
                        position = position_dodge(width = 1.2),
                        size = 3,
                        fill = "white",
                        show.legend = FALSE,
                        label.size = NA) +
    ggplot2::scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    ggplot2::scale_color_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    ggplot2::theme_classic() +
    ggplot2::xlim(c(-1, 3.55)) +
    ggplot2::labs(x = "Relative bias (%)",
                  title = "(a)                    Abundance") +
    ggplot2::theme(legend.position = "bottom",
                   legend.title = element_blank(), 
                   axis.title.y =element_blank(),
                   axis.text = element_text(size = 10, color = "black"), 
                   axis.title = element_text(size = 10, color = "black"), 
                   strip.text = element_text(size = 11, color = "black"),
                   legend.text = element_text(size = 10, color = "black"),
                   plot.title = element_text(size = 11, color = "black"),
                   panel.background = element_rect(fill = "white", color = NA), 
                   plot.background = element_rect(fill = "white", color = NA),
                   legend.margin = margin(0, 0, 0, 0), 
                   legend.box.margin = margin(-5, 0, 0, 0),
                   strip.background = element_rect(color = NA),
                   axis.line = element_line(size = 0.1, color = "black"),
                   axis.ticks = element_line(size = 0.1, color = "black")) +
    ggplot2::guides( fill = guide_legend(nrow = 1,
                                         reverse = TRUE)))

( panelb <- icm_a1 |> 
    dplyr::mutate(p_bias = ifelse(p_bias == 0, 
                                  "Same     \ndetectability", 
                                  "Count   \ndetectability\n20% lower")) |> 
    dplyr::mutate(p_bias = factor(p_bias, 
                                  levels = c(
                                    "Count   \ndetectability\n20% lower",
                                    "Same     \ndetectability"))) |> 
    dplyr::mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) |> 
    dplyr::mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
      "4x more count data",
      "Same amount of both data"))) |> 
    tibble::add_column( type = "Covariate effect") |> 
    
    ggplot2::ggplot( aes( x = mean, y = p_bias, fill = nsites_tc_fact)) +
    ggplot2::facet_wrap(~type) +
    ggplot2::geom_vline(xintercept = 0,
                        color = "gray60",
                        linetype = "dashed") +
    ggplot2::geom_boxplot(outlier.alpha = 0.2,
                          outlier.size = 0.75, 
                          size = 0.25) +
    ggplot2::labs( x = "Truth - estimate") +
    ggplot2::scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    ggplot2::ggtitle("(b)    Covariate effect") +
    ggplot2::theme_classic() +
    ggplot2::theme(legend.position = "bottom",
                   legend.title = element_blank(), 
                   axis.title.y =element_blank(),
                   axis.text.x = element_text(size = 10, color = "black"),
                   axis.text.y = element_blank(), 
                   axis.title.x = element_text(size = 10, color = "black"), 
                   strip.text = element_text(size = 11, color = "white"),
                   legend.text = element_text(size = 10, color = "black"),
                   panel.background = element_rect(fill = "white", color = NA), 
                   plot.background = element_rect(fill = "white", color = NA),
                   legend.margin = margin(0, 18, 0, 0), 
                   plot.title = element_text(size = 11, color = "black"),
                   legend.box.margin = margin(-5, 18, 0, 0),
                   strip.background = element_rect(color = NA),
                   axis.line = element_line(size = 0.1, color = "black"),
                   axis.ticks.x = element_line(size = 0.1, color = "black"),
                   axis.ticks.y = element_blank()) +
    ggplot2::guides(fill = guide_legend(nrow = 1, 
                                        reverse = TRUE)))

panela + panelb + plot_layout( nrow = 1, widths = c(2, 1), guides = "collect") &
  theme(legend.position='bottom',
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)) 

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_02.png",
  width = 6, 
  height = 4, 
  units = "in", 
  dpi = 300
)

# Table S2 - accuracy and precision of ICM's covariate estimates
icm_a1 |> 
  dplyr::mutate(p_bias = ifelse(p_bias == 0, 
                                "Same detectability", 
                                "Count detectability 20% lower")) |> 
  dplyr::mutate(p_bias = factor(p_bias, 
                                levels = c(
                                  "Same detectability",
                                  "Count detectability 20% lower"
                                ))) |> 
  dplyr::mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, 
                                        "Same amount of both data", 
                                        "4x more count data")) |> 
  dplyr::mutate(nsites_tc_fact = factor(nsites_tc_fact, 
                                        levels = c(
                                          "Same amount of both data",
                                          "4x more count data"
                                        ))) |> 
  dplyr::group_by( nsites_tc_fact, p_bias ) |>
  dplyr::mutate(n = n()) |> 
  dplyr::summarise( unbiased = round( 100 * ( sum(contain0) / n), 1),
                    biased = round( 100 - unbiased, 1),
                    mean_sd = round( mean(sd), 2)) |> 
  dplyr::distinct() |> 
  dplyr::arrange(nsites_tc_fact, p_bias) |> 
  dplyr::mutate(unbiased = paste0(unbiased, "%"), 
                biased = paste0(biased, "%")) |> 
  dplyr::rename( `Data amount` = nsites_tc_fact, 
                 `Detectability` = p_bias, 
                 Unbiased = unbiased, 
                 Biased = biased, 
                 SD = mean_sd )


setwd(here::here("results"))
load("simulation_alternative_model_results_v01.RData")

icm_for_comparison <- icm_info |>
  dplyr::group_by( simrep ) |> 
  dplyr::arrange(simrep, totDS) |> 
  dplyr::filter(! totDS == max(totDS) ) |> 
  dplyr::slice(1, n()) |> 
  dplyr::mutate(species = ifelse( totDS == min(totDS), "rare", "common")) |> 
  dplyr::full_join( distinct( dplyr::select(icm_n_ds_rb, simrep, nsites_tc_fact, p_bias))) |> 
  dplyr::filter(nsites_tc_fact == 1 & p_bias == 0) |> 
  tibble::add_column(model = "ICM") |> 
  dplyr::select( model, species, sp, simrep ) |> 
  dplyr::left_join(icm_a1) |> 
  dplyr::ungroup() |> 
  dplyr::mutate(contain0 = ifelse(contain0 == 1, "contains 0", "does not contain 0")) |> 
  dplyr::select(model, species, mean, sd, l95, u95, contain0 )

key <- tibble::tribble(
  ~model, ~name,
  "ICM", "Community, integrated",
  "CC", "Community, counts",
  "CDS", "Community, distance sampling",
  "ISS", "Single species, integrated", 
  "SSC", "Single species, counts", 
  "SSDS", "Single species, distance sampling")

alpha1_truth_minus_estimate |> 
  dplyr::full_join(icm_for_comparison) |> 
  dplyr::full_join(key) |> 
  dplyr::mutate(species = ifelse(species == "common", "Common species", "Rare species")) |> 
  dplyr::mutate(species = factor(species, levels = c("Rare species", "Common species"))) |> 
  dplyr::mutate( name = factor(name, 
                               levels = c(
                                 "Community, integrated",
                                 "Community, distance sampling", 
                                 "Community, counts",
                                 "Single species, integrated", 
                                 "Single species, distance sampling",
                                 "Single species, counts"))) |> 
  ggplot2::ggplot(
    aes( x = mean, y = rev(name), fill = species)) +
  ggplot2::scale_y_discrete(expand = c(0, 0)) +
  ggplot2::geom_vline(xintercept = 0,
                      color = MetPalettes$Hiroshige[[1]][c(1)]) +
  ggplot2::geom_boxplot( outlier.alpha = 0.2,
                         outlier.size = 0.75, 
                         size = 0.25) +
  ggplot2::scale_color_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  ggplot2::scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Covariate effect: truth - estimate") +
  ggplot2::theme(
    legend.position = "bottom",
    legend.title = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text( size = c(10, 10, 10, 10, 10, 14), color = "black",
                                face = c("plain", "plain", "plain", "plain", "plain", "bold")),
    axis.title = element_text(size = 11, color = "black"), 
    legend.text = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(0, 12.5, 0, 0), 
    legend.box.margin = margin(-10, 12.5, 0, 0),
    axis.line = element_line(size = 0.1, color = "black"), 
    axis.ticks = element_line(size = 0.1, color = "black")) + 
  ggplot2::guides(color = guide_legend(reverse = TRUE)) +
  ggplot2::scale_y_discrete(limits = rev)

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_03.png",
  width = 4.75, 
  height = 4, 
  units = "in", 
  dpi = 300
)

# relative precision of models
# Table S3
alpha1_truth_minus_estimate |> 
  dplyr::full_join(icm_for_comparison) |> 
  dplyr::full_join(key) |> 
  dplyr::mutate(species = ifelse(species == "common", "Common species", "Rare species")) |> 
  dplyr::mutate(species = factor(species, levels = c("Rare species", "Common species"))) |> 
  dplyr::mutate( name = factor(name, 
                               levels = c(
                                 "Community, integrated",
                                 "Community, distance sampling", 
                                 "Community, counts",
                                 "Single species, integrated", 
                                 "Single species, distance sampling",
                                 "Single species, counts"))) |> 
  dplyr::group_by(name, species) |> 
  dplyr::summarise(mean_sd = mean(sd)) |> 
  dplyr::ungroup() |> 
  dplyr::group_by(species) |> 
  dplyr::mutate( rel_sd = mean_sd / min(mean_sd) ) |> 
  dplyr::select(name, species, rel_sd) |> 
  tidyr::pivot_wider(names_from = species, values_from = rel_sd) |> 
  dplyr::rename(Model = name)
