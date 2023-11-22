library(here)
library(tidyverse)
library(MetBrewer)

setwd(here::here("results"))

load("ic.RData")

ic_clean <- ic |> 
  dplyr::mutate(diff = mean - truth)

rb_label <- ic_clean |> 
  dplyr::filter(param %in% c(
    "alpha0",
    "alpha1",
    "gamma0_ds",
    "gamma0_c",
    "N_DS", "N_TC")) |> 
  
  dplyr::mutate(param = ifelse(param == "alpha0", "Abundance intercept", 
                               ifelse(param == "alpha1", "Abundance covariate coefficient",
                                      ifelse(param == "gamma0_ds", "Detection intercept (distance sampling)", 
                                             ifelse(param == "gamma0_c", "Detection intercept (counts)",
                                                    ifelse(param == "N_DS", "Abundance (distance sampling)", 
                                                           "Abundance (counts)")))))) |> 
  dplyr::mutate(param = factor(param, 
                               levels = c(
                                 "Abundance intercept", 
                                 "Abundance covariate coefficient",
                                 "Detection intercept (distance sampling)",
                                 "Detection intercept (counts)",
                                 "Abundance (distance sampling)",
                                 "Abundance (counts)"))) |> 
  dplyr::group_by(param) |> 
  dplyr::summarise( rb = paste0( "Relative bias\n", sprintf("%.2f", round(100 * sum(diff) / sum(abs(truth)), 2)), "%")) |> 
  dplyr::arrange(param) |> 
  tibble::add_column( y = c(3379.6, 1940.4, 2088.1, 1539.3, 214349.1, 411098.8),
                      diff = c(1.62511451524311, 0.482292718683334, 0.571734517221783, 0.745914880139566, 
                               11.997, 11.9998)) |> 
  dplyr::mutate(param = factor(param, 
                               levels = c(
                                 "Detection intercept (distance sampling)",
                                 "Detection intercept (counts)",
                                 "Abundance intercept", 
                                 "Abundance covariate coefficient",
                                 "Abundance (distance sampling)",
                                 "Abundance (counts)")))
ic_clean |> 
  dplyr::filter(param %in% c(
    "alpha0",
    "alpha1",
    "gamma0_ds",
    "gamma0_c",
    "N_DS", "N_TC")) |> 
  
  dplyr::mutate(param = ifelse(param == "alpha0", "Abundance intercept", 
                               ifelse(param == "alpha1", "Abundance covariate coefficient",
                                      ifelse(param == "gamma0_ds", "Detection intercept (distance sampling)", 
                                             ifelse(param == "gamma0_c", "Detection intercept (counts)",
                                                    ifelse(param == "N_DS", "Abundance (distance sampling)", 
                                                           "Abundance (counts)")))))) |> 
  dplyr::mutate(param = factor(param, 
                               levels = c(
                                 "Detection intercept (distance sampling)",
                                 "Detection intercept (counts)",
                                 "Abundance intercept", 
                                 "Abundance covariate coefficient",
                                 "Abundance (distance sampling)",
                                 "Abundance (counts)"))) |> 
  # filter(abs(diff) < 20) |> 
  # group_by(param) |>
  # mutate(bin = cut(diff, breaks = 30)) |>
  # count(bin) |>
  # summarise(max_n = max(n),
  #           max_n70 = max_n*0.70) |>
  # pull(max_n70) |>
  # dput()
  
  dplyr::filter(abs(diff) < 20) |> 
  ggplot2::ggplot(aes(x = diff)) +
  ggplot2::facet_wrap(~param, scales = "free", ncol = 2) + 
  ggplot2::geom_histogram(fill = MetBrewer::MetPalettes$Hiroshige[[1]][8],
                          color = MetBrewer::MetPalettes$Hiroshige[[1]][9]) +
  ggplot2::xlab("Posterior mean - true value") +
  
  ggplot2::geom_label(data = rb_label,
                      aes(label = rb,
                          y = y),
                      size = 3,
                      fill = "white",
                      color = MetBrewer::MetPalettes$Hiroshige[[1]][9],
                      show.legend = FALSE,
                      label.size = NA) + 
  
  ggplot2::geom_vline(xintercept = 0,
                      color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                      linetype = "dashed",
                      linewidth = 1) +
  ggplot2::theme_classic() +
  ggplot2::theme(strip.background = element_blank(),
                 axis.text.y = element_blank(), 
                 axis.ticks.y = element_blank(),
                 axis.line.y = element_blank(),
                 axis.title.y = element_blank(),
                 axis.text.x = element_text(size = 10, color = "black"),
                 axis.title.x = element_text(size = 11, color = "black"),
                 axis.line.x = element_line(color = "black", size = 0.2),
                 axis.ticks.x = element_line(color = "black", size = 0.2),
                 strip.text = element_text(color = "black", size = 9))

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_02.png",
  width = 5,
  height = 6,
  units = "in",
  dpi = 300
)