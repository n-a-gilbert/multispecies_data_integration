# 12 October 2023
# This script creates Figures S1 & S2
# Visualizes the detection function and abundance-covariate relationship for the simulated community

library(tidyverse)

nsp <- 1000
mu_gamma <- 5.5
sd_gamma <- 0.25
gamma <- rnorm(nsp, mu_gamma, sd_gamma)
omega <- exp(gamma)

mu_gamma_c <- 5.0
sd_gamma_c <- 0.25
gamma_c <- rnorm(nsp, mu_gamma_c, sd_gamma_c)
omega_c <- exp(gamma_c)

distances <- seq(from = 0, to = 1000, by = 10)

sp_df <- tibble::tibble(
  sp = 1:nsp, 
  omega = omega,
  omega_c = omega_c)

com <- expand.grid(
  distances = seq(from = 0, to = 1000, by = 10)) |> 
  dplyr::mutate(p = exp( - distances * distances / (2 * exp(mu_gamma) * exp(mu_gamma))),
                pc = exp( - distances * distances / (2 * exp(mu_gamma_c) * exp(mu_gamma_c)))) |> 
  tidyr::pivot_longer(p:pc, names_to = "data", values_to = "value") |> 
  dplyr::mutate( data = ifelse(grepl("c", data), "Count", "Distance sampling")) |> 
  dplyr::mutate(data = factor(data, levels = c("Distance sampling", "Count")))

spp <- expand.grid(sp = 1:nsp, 
                   distances = seq(from = 0, to = 1000, by = 10)) |> 
  dplyr::full_join(sp_df) |> 
  dplyr::select(sp, omega, omega_c, distances) |> 
  tibble::as_tibble() |>
  tidyr::pivot_longer(omega:omega_c, names_to = "data", values_to = "value") |> 
  dplyr::mutate(p = exp ( - distances * distances / (2 * value * value))) |> 
  dplyr::mutate( data = ifelse(grepl("c", data), "Count", "Distance sampling")) |> 
  dplyr::mutate(data = factor(data, levels = c("Distance sampling", "Count")))

( det_fun <- ggplot2::ggplot() + 
    ggplot2::geom_line(data = spp, aes(x = distances, y = p, group = factor(sp)),
                       color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                       alpha = 0.1) + 
    ggplot2::facet_wrap(~data) +
    
    ggplot2::geom_line(data = com, aes(x = distances, y = value), color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                       size = 1.5) +
    ggplot2::labs(x = "Distance (m) from observer", 
                  y = "Detection probability") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = element_text(size = 8, color = "black"),
                   axis.title = element_text(size = 9, color = "black")) )

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s2.png",
  det_fun,
  width = 7,
  height = 3.25,
  units = "in",
  dpi = 300)

mu_alpha0 <- 0.87
sigma_alpha0 <- 1.95
mu_alpha1 <- 0.05    
sigma_alpha1 <- 0.25 

alpha0 <- rnorm(nsp, mean = mu_alpha0, sigma_alpha0)
alpha1 <- rnorm(nsp, mu_alpha1, sigma_alpha1)

sp_df <- tibble::tibble(
  sp = 1:nsp, 
  alpha0 = alpha0, 
  alpha1 = alpha1)

spp_cov <- expand.grid(sp = 1:nsp,
                       x = seq(from = -2, to = 2, by = 0.1)) |> 
  dplyr::full_join(sp_df) |> 
  dplyr::mutate(lp = alpha0 + alpha1 * x)   

com_cov <- expand.grid(x = seq(from = -2, to = 2, by = 0.1)) |> 
  tibble::add_column(alpha0 = mu_alpha0,
                     alpha1 = mu_alpha1) |> 
  dplyr::mutate(lp = alpha0 + alpha1 * x)

( cov_plot <- ggplot2::ggplot() + 
    ggplot2::geom_line(data = spp_cov, aes(x = x, y = lp, group = sp),
                       color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                       alpha = 0.1) + 
    ggplot2::geom_line(data = com_cov, aes(x = x, y = lp), 
                       color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                       size = 1.5) +
    ggplot2::labs(x = "Covariate", 
                  y = "Log( Abundance )") +
    ggplot2::theme_minimal() +
    ggplot2::theme(axis.text = element_text(size = 8, color = "black"),
                   axis.title = element_text(size = 9, color = "black")) )

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s1.png", 
  cov_plot, 
  width = 4.5, 
  height = 3.25, 
  units = "in", 
  dpi = 300)
