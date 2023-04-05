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

sp_df <- tibble(
  sp = 1:nsp, 
  omega = omega)

com <- expand.grid(
  distances = seq(from = 0, to = 1000, by = 10)) %>% 
  mutate(p = exp( - distances * distances / (2 * exp(mu_gamma) * exp(mu_gamma))))

spp <- expand.grid(sp = 1:nsp, 
                   distances = seq(from = 0, to = 1000, by = 10)) %>% 
  full_join(sp_df) %>% 
  dplyr::select(sp, omega, distances) %>% 
  as_tibble() %>% 
  mutate(p = exp ( - distances * distances / (2 * omega * omega)))

det_fun <- ggplot() + 
  geom_line(data = spp, aes(x = distances, y = p, group = factor(sp)), alpha = 0.1) + 
  geom_line(data = com, aes(x = distances, y = p), color = "red", size = 1.5) +
  labs(x = "Distance (m) from transect", 
       y = "Detection probability") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black")) 

setwd(here::here("figures"))
ggsave(
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

sp_df <- tibble(
  sp = 1:nsp, 
  alpha0 = alpha0, 
  alpha1 = alpha1)

spp_cov <- expand.grid(sp = 1:nsp,
                       x = seq(from = -2, to = 2, by = 0.1)) %>% 
  full_join(sp_df) %>% 
  mutate(lp = alpha0 + alpha1 * x)   

com_cov <- expand.grid(x = seq(from = -2, to = 2, by = 0.1)) %>% 
  add_column(alpha0 = mu_alpha0,
             alpha1 = mu_alpha1) %>% 
  mutate(lp = alpha0 + alpha1 * x)

cov_plot <- ggplot() + 
  geom_line(data = spp_cov, aes(x = x, y = lp, group = sp), alpha = 0.1) + 
  geom_line(data = com_cov, aes(x = x, y = lp), color = "red", size = 1.5) +
  labs(x = "Covariate", 
       y = "Log(expected number of groups)") +
  theme_minimal() +
  theme(axis.text = element_text(size = 8, color = "black"),
        axis.title = element_text(size = 9, color = "black")) 

setwd(here::here("figures"))
ggsave(
  "figure_s1.png", 
  cov_plot, 
  width = 4.5, 
  height = 3.25, 
  units = "in", 
  dpi = 300)