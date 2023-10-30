# 30 October 2023
# Code to demonstrate prior predictive check for prior for scale parameter intercept

library(tidyverse)
library(here)

mids <- seq(from = 12.5, to = 987.5, by = 25)
gamma0 <- seq(from =  0, to = 10, by = 1)

d <- list(list())
for(i in 1:length(gamma0)){
  p <- exp(-mids*mids/(2 * exp(gamma0[i])*exp(gamma0[i])))
  d[[i]] <- tibble::tibble(d = mids,
                           p = p) |> 
    tibble::add_column(gamma0 = gamma0[i])
}

dplyr::bind_rows(d) |> 
  dplyr::mutate(gamma0 = paste0("Gamma0: ", gamma0)) |> 
  dplyr::mutate(gamma0 = factor(gamma0, 
                         levels = paste0("Gamma0: ", 0:10))) |> 
  
  ggplot2::ggplot(aes(x = d, y = p)) +
  ggplot2::facet_wrap(~gamma0) +
  ggplot2::geom_line(color = "gray60") +
  ggplot2::geom_point(size = 1.5) +
  ggplot2::labs(x = "Distance", 
       y = "Detection probability") +
  ggplot2::theme_minimal()

setwd(here::here("figures"))
ggsave(
  "figure_s3.png", 
  width = 6, 
  height = 4, 
  units = "in", 
  dpi = 300
)
