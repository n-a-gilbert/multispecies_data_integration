library(here)
library(tidyverse)
library(sf)
library(maptools)
library(cowplot)
library(MetBrewer)

setwd(here::here("data/Shapefiles/reserve/"))

res <- sf::st_read("Reserve_polygon.shp") |> 
  dplyr::select(Id, geometry) |> 
  dplyr::mutate(region = ifelse(Id == 0, "Talek", 
                                ifelse(Id == 2, "Mara Triangle", "River")))

setwd(here::here("data/Shapefiles/DS"))

ds <-sf::st_read("DS_10kmpersite.shp") |> 
  tibble::add_column(data = "Distance sampling") |> 
  dplyr::select(data, geometry)

setwd(here::here("data/Shapefiles/Transects"))

tc <- sf::st_read("Transects.shp") |> 
  tibble::add_column(data = "Counts") |> 
  dplyr::select(data, geometry)

surveys <- rbind(ds, tc)

fillpal <- MetBrewer::MetPalettes$Hiroshige[[1]][c(4, 6, 6)]
linepal <- MetBrewer::MetPalettes$Hiroshige[[1]][c(1, 10)]

talek <- tibble::tibble(
  x = 35.214983,
  y = -1.441084) |> 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326) |> 
  sf::st_transform(st_crs(res)) |> 
  tibble::add_column(name = "Talek town") |> 
  ( function(z) dplyr::mutate(z, 
                              x = sf::st_coordinates(z)[,1] + 6000,
                              y = sf::st_coordinates(z)[,2] + 1000))()

mt_label <- tibble::tibble(
  x = c(34.909227),
  y = c(-1.446516), 
  name = c("Mara Triangle")) |> 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326) |> 
  sf::st_transform( sf::st_crs(res)) |> 
  ( function(z) dplyr::mutate(z, 
                              x = sf::st_coordinates(z)[,1],
                              y = sf::st_coordinates(z)[,2]))()

tr_label <- tibble::tibble(
  x = c(35.161652),
  y = c(-1.546157), 
  name = c("Talek region")) |> 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326) |> 
  sf::st_transform(sf::st_crs(res)) |> 
  ( function(z) dplyr::mutate(z, 
                              x = sf::st_coordinates(z)[,1],
                              y = sf::st_coordinates(z)[,2]))()

tt_label <- tibble::tibble(
  x = c(35.27),
  y = c(-1.425038), 
  name = c("Talek town")) |> 
  sf::st_as_sf(coords = c("x", "y"), 
               crs = 4326) |> 
  sf::st_transform(sf::st_crs(res)) |> 
  ( function(z) dplyr::mutate(z, 
                              x = sf::st_coordinates(z)[,1],
                              y = sf::st_coordinates(z)[,2]))()

res_plot <-
  ggplot2::ggplot(  ) + 
  ggplot2::geom_sf(data = res, 
                   aes( geometry = geometry, fill = region),
                   color = NA ) +
  
  ggplot2::geom_sf( data = filter(surveys, data != "Distance sampling"), 
                    aes( geometry = geometry, 
                         color = data),
                    lwd = 3) +  
  ggplot2::geom_sf( data = filter(surveys, data == "Distance sampling"), 
                    aes( geometry = geometry, 
                         color = data),
                    lwd = 1) +
  ggplot2::geom_sf(data = talek, 
                   aes(geometry = geometry), 
                   size = 4) +
  ggplot2::geom_text(data = mt_label, 
                     aes(x = x, 
                         y = y, 
                         label = name),
                     angle = 330,
                     color = MetBrewer::MetPalettes$Hiroshige[[1]][2]) +
  ggplot2::geom_text(data = tr_label, 
                     aes(x = x, 
                         y = y, 
                         label = name), 
                     color = MetBrewer::MetPalettes$Hiroshige[[1]][8]) +
  ggplot2::geom_text(data = tt_label,
                     aes(x = x, 
                         y = y, 
                         label = name), 
                     color = "black",
                     size = 3.1) +
  # scale_size_manual(values = c(10, 5)) +
  ggplot2::scale_color_manual(values = linepal) + 
  ggplot2::theme_minimal() +
  ggplot2::scale_fill_manual(
    values = fillpal,
    limits = c("Mara Triangle", "Talek"))  +  
  ggplot2::theme(
    legend.position = "none", 
    axis.text = element_text(color = "black", 
                             size = 9),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.5, "cm"),
    axis.title = element_blank(),
    plot.background = element_rect(color = NA, 
                                   fill = "white"),
    panel.background = element_rect(color = NA, 
                                    fill = "white"))

data("wrld_simpl")

wrld_sf <- sf::st_as_sf(wrld_simpl)

ken <- wrld_sf |>
  dplyr::filter(NAME == "Kenya") |> 
  sf::st_transform(sf::st_crs(res))

res_bbox <- sf::st_bbox(res) |> 
  ( function(x) sf::st_as_sfc(x))() |> 
  ( function(x) sf::st_buffer(x, 30000))()

ken_label <- sf::st_centroid(ken) |> 
  ( function(z) dplyr::mutate(z,
                              x = sf::st_coordinates(z)[,1],
                              y = sf::st_coordinates(z)[,2]))() |> 
  sf::st_drop_geometry() |> 
  tibble::add_column(name = "Kenya") |> 
  dplyr::ungroup() |> 
  dplyr::select(x, y, name)

inset <-
  ggplot2::ggplot() + 
  ggplot2::geom_sf(data = ken, aes(geometry = geometry),
                   size = 0.1) + 
  ggplot2::geom_sf(data = res, aes(geometry = geometry), 
                   color = "black", 
                   fill = "black") +
  ggplot2::geom_sf(data = res_bbox, aes(geometry = geometry), 
                   color = "red", 
                   fill = NA,
                   size = 0.75) +
  ggplot2::geom_text(data = ken_label, aes(x = x, y = y, 
                                           label = name),
                     color = "black") +
  ggplot2::theme_void() +
  ggplot2::theme(plot.background = element_rect(fill = "white",
                                                color = "gray20"))

leg_data <- tibble::tibble( x = c(1, 2), 
                            xend = c(2, 3), 
                            y = c(1, 2), 
                            yend = c( 2, 3),
                            data = factor(c("Counts", "Distance sampling"),
                                          levels = c("Distance sampling", "Counts")))


for_legend <-
  ggplot2::ggplot() +
  ggplot2::geom_segment(data = filter(leg_data, 
                                      data == "Counts"), 
                        aes(x = x, 
                            xend = xend, 
                            y = y, 
                            yend = yend, 
                            color = data, 
                            lwd = data),
                        lwd = 3) +
  
  ggplot2::geom_segment(data = filter(leg_data, 
                                      data == "Distance sampling"), 
                        aes(x = x, 
                            xend = xend, 
                            y = y, 
                            yend = yend, 
                            color = data,
                            lwd = data),
                        lwd = 1) +
  ggplot2::scale_color_manual(values = linepal) +
  ggplot2::theme_minimal() +
  ggplot2::theme(legend.title = element_blank(),
                 legend.background = element_rect(fill = "white", color = NA))

leg <- cowplot::get_legend(for_legend)

cowplot::ggdraw() + 
  cowplot::draw_plot(res_plot) + 
  cowplot::draw_plot(inset,
                     x = 0.65,
                     y = 0.67,
                     width = 0.3,
                     height = 0.32) + 
  cowplot::draw_plot(leg,
                     x = -0.15, 
                     y = -0.25)

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s4.png",
  width = 7.5, 
  height = 5, 
  units = "in", 
  dpi = 300)
