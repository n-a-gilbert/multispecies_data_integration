# 15 March 2023
# This script creates Figures S1 & S2
# Visualizes the detection function and number-of-groups/covariate relationship for the simulated community

library(tidyverse)

nsp <- 1000
mu_gamma <- 5.5
sd_gamma <- 0.25
gamma <- rnorm(nsp, mu_gamma, sd_gamma)
omega <- exp(gamma)
distances <- seq(from = 0, to = 1000, by = 10)

sp_df <- tibble::tibble(
  sp = 1:nsp, 
  omega = omega)

com <- expand.grid(
  distances = seq(from = 0, to = 1000, by = 10)) |> 
  dplyr::mutate(p = exp( - distances * distances / (2 * exp(mu_gamma) * exp(mu_gamma))))

spp <- expand.grid(sp = 1:nsp, 
                   distances = seq(from = 0, to = 1000, by = 10)) |> 
  dplyr::full_join(sp_df) |> 
  dplyr::select(sp, omega, distances) |> 
  tibble::as_tibble() |> 
  dplyr::mutate(p = exp ( - distances * distances / (2 * omega * omega)))

det_fun <- ggplot2::ggplot() + 
  ggplot2::geom_line(data = spp, aes(x = distances, y = p, group = factor(sp)),
                     color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                     alpha = 0.1) + 
  ggplot2::geom_line(data = com, aes(x = distances, y = p), color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                     size = 1.5) +
  ggplot2::labs(x = "Distance (m) from observer", 
                y = "Detection probability") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text = element_text(size = 8, color = "black"),
                 axis.title = element_text(size = 9, color = "black")) 

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s2.png", 
  det_fun, 
  width = 4.5, 
  height = 3.25, 
  units = "in", 
  dpi = 300)

mu_alpha0 <- -0.65
sigma_alpha0 <- 1.7
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

cov_plot <- ggplot2::ggplot() + 
  ggplot2::geom_line(data = spp_cov, aes(x = x, y = lp, group = sp),
                     color = MetBrewer::MetPalettes$Hiroshige[[1]][10],
                     alpha = 0.1) + 
  ggplot2::geom_line(data = com_cov, aes(x = x, y = lp), 
                     color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                     size = 1.5) +
  ggplot2::labs(x = "Covariate", 
                y = "Log(expected number of groups)") +
  ggplot2::theme_minimal() +
  ggplot2::theme(axis.text = element_text(size = 8, color = "black"),
                 axis.title = element_text(size = 9, color = "black")) 

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s1.png", 
  cov_plot, 
  width = 4.5, 
  height = 3.25, 
  units = "in", 
  dpi = 300)
