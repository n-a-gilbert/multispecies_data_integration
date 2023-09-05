# 16 March 2023
# Figure 5
# Visualize differences in number-of-groups and group-size between regions

library(here)
library(tidyverse)
library(MetBrewer)
library(MCMCvis)
library(patchwork)
library(reshape)

setwd(here::here("results"))
load("herbivore_case_study_results_v04.RData")

sp_key <- tibble::tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's gazelle", "hartebeest",
                "impala",  "Thomson's gazelle", "topi",
                "warthog", "waterbuck" ),
  mass = c(592666, 562592.7, 3824540, 964654.7, 55464.46,
           160937.9,  52591.69, 22907.43, 136000.3,
           82499.99, 204393.5)) |> 
  dplyr::arrange(mass) |> 
  dplyr::mutate( sp_name = factor(sp_name, levels = unique(sp_name)))

sp_n <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("sp", "region", "iter")) |> 
  dplyr::bind_rows(.id = "param") |> 
  dplyr::mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  dplyr::mutate( value = exp(value)) %>%
  tibble::as_tibble() |> 
  tidyr::pivot_wider(names_from = param, values_from = value) |> 
  dplyr::mutate( N = `Group size` * `Number of groups`) %>%
  dplyr::select(sp, region, iter, value = N) |> 
  dplyr::mutate(value = log(value)) |> 
  dplyr::full_join(sp_key) |> 
  dplyr::group_by(sp_name, region) |> 
  dplyr::summarise( mean = mean(value), 
                    l95 = quantile(value, c(0.025)), 
                    u95 = quantile(value, c(0.975)))

sp_sig <-
  lapply(
    X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
    FUN = reshape::melt, 
    varnames = c("sp", "region", "iter")) |> 
  dplyr::bind_rows(.id = "param") |> 
  dplyr::mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  dplyr::mutate( value = exp(value)) %>%
  tibble::as_tibble() |> 
  tidyr::pivot_wider(names_from = param, values_from = value) |> 
  dplyr::mutate( N = `Group size` * `Number of groups`) %>%
  dplyr::select(sp, region, iter, N) %>%
  tidyr::pivot_wider(names_from = region, values_from = N) %>%
  dplyr::mutate(diff = Talek - Mara) %>%
  dplyr::group_by(sp) %>%
  dplyr::summarise( mean = mean(diff),
                    l95 = quantile(diff, c(0.025)),
                    u95 = quantile(diff, c(0.975))) |> 
  dplyr::full_join(sp_key) |> 
  dplyr::mutate(region = ifelse( mean > 0, "Talek", "Mara")) |> 
  dplyr::mutate(sig = ifelse( l95 < 0 & u95 > 0, NA, "*")) |> 
  dplyr::select(sp_name, sig, region) |> 
  tibble::add_column( mean = -3.45)

com_n <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("mu_beta0", "mu_alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("region", "iter")) |> 
  dplyr::bind_rows(.id = "param") |> 
  tibble::as_tibble() |> 
  tidyr::separate( region, into = c("junk", "region"), sep = "\\[") |> 
  dplyr::mutate(region = parse_number(region)) |> 
  dplyr::mutate( param = ifelse(param == "mu_alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate( value = exp(value)) |> 
  dplyr::select( param, region, iter, value ) %>%
  tidyr::pivot_wider(names_from = param, values_from = value) %>%
  dplyr::mutate( value = `Group size` * `Number of groups`) %>%
  dplyr::mutate(value = log(value)) |> 
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  
  dplyr::group_by( region) |> 
  dplyr::summarise( mean = mean(value), 
                    l95 = quantile(value, c(0.025)), 
                    u95 = quantile(value, c(0.975))) 

( nplot <- ggplot2::ggplot() +
    ggplot2::geom_rect(
      data = com_n,
      aes( xmin = l95,
           xmax = u95,
           fill = region),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.2) +
    ggplot2::geom_vline(
      data = com_n,
      aes(xintercept = mean,
          color = region),
      size = 1,
      alpha = 0.2) +
    ggplot2::geom_errorbar(
      data = sp_n, 
      aes(y = sp_name, 
          xmin = l95, 
          xmax = u95, 
          color = region),
      width = 0,
      size = 1,
      position = position_dodge(width = 0.5)) +
    ggplot2::geom_point(
      data = sp_n, 
      aes(y = sp_name,
          x = mean, 
          color = region),
      size = 3,
      position = position_dodge(width = 0.5)) +
    ggplot2::geom_text( data = sp_sig,
                        aes(x = mean,
                            y = sp_name,
                            label = sig,
                            color = region),
                        size = 6,
                        vjust = 0.75,
                        fontface = "bold",
                        show.legend = FALSE) +
    ggplot2::scale_fill_manual("Region",
                               values = MetBrewer::MetPalettes$Tam[[1]][c(3, 7)]) +
    ggplot2::scale_color_manual( "Region",
                                 values = MetBrewer::MetPalettes$Tam[[1]][c(3, 7)]) +
    ggplot2::scale_x_continuous( limits = c(-3.5, 5.5),
                                 expand = c(0.02, 0)) +
    ggplot2::theme_minimal() +
    ggplot2::labs( x = " Log ( Abundance ) ",
                   title = "(a)") + 
    ggplot2::theme( panel.grid = element_blank(),
                    axis.title.y = element_blank(),
                    legend.text = element_text(color = "black", size = 9), 
                    legend.title = element_text(color = "black", size = 10),
                    legend.position = c(0.9, 0.9),
                    plot.title = element_text(color = "black", size = 9),
                    axis.text = element_text(color = "black", size = 9),
                    axis.title.x = element_text(color = "black", size = 10),
                    axis.line = element_line(color = "black", size = 0.1),
                    axis.ticks.x = element_line(color = "black", size = 0.1),
                    plot.background = element_rect(color = NA, fill = "white"), 
                    panel.background = element_rect(color = NA, fill = "white")) +
    ggplot2::guides(color = guide_legend(reverse = TRUE), 
                    fill = guide_legend(reverse = TRUE)) )

species_ng_gs <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("sp", "region", "iter")) |> 
  dplyr::bind_rows(.id = "param") |> 
  dplyr::mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  dplyr::mutate( value = exp(value)) |> 
  tidyr::pivot_wider(names_from = region, values_from = value) |> 
  dplyr::mutate( diff = Talek - Mara) |> 
  dplyr::group_by(sp, param) |> 
  dplyr::summarise( mean = mean(diff), 
                    l95 = quantile(diff, c(0.025)), 
                    u95 = quantile(diff, c(0.975))) |> 
  dplyr::full_join(sp_key) |> 
  dplyr::mutate( param = factor(param, 
                                levels = c("Number of groups", "Group size")))

com_ng_gs <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("mu_beta0", "mu_alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("region", "iter")) |> 
  dplyr::bind_rows() %>%
  tibble::as_tibble() |> 
  dplyr::mutate( param = ifelse( grepl("alpha0", region), "Number of groups", "Group size")) |> 
  tidyr::separate(region, into = c("junk", "region"), sep = "\\[") |> 
  dplyr::mutate(region = parse_number(region)) |> 
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) |> 
  dplyr::mutate( value = exp(value)) |> 
  tidyr::pivot_wider(names_from = region, values_from = value) |> 
  dplyr::mutate( diff = Talek - Mara)|> 
  dplyr::group_by(param) |> 
  dplyr::summarise( mean = mean(diff), 
                    l95 = quantile(diff, c(0.025)), 
                    u95 = quantile(diff, c(0.975))) |> 
  dplyr::mutate( param = factor(param, 
                                levels = c("Number of groups", "Group size")))

labs <- tibble::tibble(
  param = c("Number of groups", "Group size"),
  sp_name = c("hartebeest", "hartebeest"),
  mean = c(10, 10),
  label = c("More groups in Talek", "Larger groups in Talek")) |> 
  dplyr::mutate(sp_name = factor(sp_name,
                                 levels = c("Grant's gazelle", "hartebeest", "eland", "giraffe", "elephant", "waterbuck",
                                            "warthog", "buffalo", "topi", "Thomson's gazelle", "impala"))) |> 
  dplyr::mutate( param = factor(param, 
                                levels = c("Number of groups", "Group size")))

( diffplot <- ggplot2::ggplot(  ) +
    ggplot2::geom_rect( data = com_ng_gs, 
                        aes( xmin = l95, 
                             xmax = u95), 
                        ymin = -Inf, 
                        ymax = Inf,
                        alpha = 0.2,
                        fill = MetBrewer::MetPalettes$Hiroshige[[1]][7]) +
    ggplot2::geom_vline(data = com_ng_gs, 
                        aes( xintercept = mean),
                        color = MetBrewer::MetPalettes$Hiroshige[[1]][7],
                        alpha = 0.2, 
                        size = 1) +
    ggplot2::geom_vline( xintercept = 0,
                         color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                         linetype = "dashed" ) +
    ggplot2::geom_errorbar( data = species_ng_gs, 
                            aes( y = sp_name, xmin = l95, xmax = u95),
                            width = 0, size = 1,
                            color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    ggplot2::geom_point( data = species_ng_gs, 
                         aes(y = sp_name, x = mean), 
                         size  = 3, 
                         color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    ggplot2::geom_text(
      data = labs, 
      aes( y = sp_name, 
           x = mean, 
           label = label),
      size = 2.5,
      color = MetBrewer::MetPalettes$Hiroshige[[1]][1]) +
    ggplot2::facet_wrap(~param, ncol = 1)  +
    ggplot2::theme_minimal() +
    ggplot2::scale_x_continuous(breaks = c(-5, 0, 5, 10, 15)) +
    ggplot2::labs( x = " Talek - Mara ",
                   title = "(b)") + 
    ggplot2::theme( panel.grid = element_blank(),
                    axis.title.y = element_blank(),
                    plot.title = element_text(color = "black", size = 9),
                    axis.text = element_text(color = "black", size = 9),
                    axis.title.x = element_text(color = "black", size = 10),
                    axis.line = element_line(color = "black", size = 0.1),
                    axis.ticks.x = element_line(color = "black", size = 0.1),
                    plot.background = element_rect(color = NA, fill = "white"), 
                    panel.background = element_rect(color = NA, fill = "white")) )


nplot + diffplot + plot_layout(widths = c(1, 1))

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_04.png", 
  width = 7.15,
  height = 5, 
  units = "in", 
  dpi = 300)