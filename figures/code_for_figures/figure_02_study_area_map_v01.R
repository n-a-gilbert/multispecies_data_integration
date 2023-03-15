library(here)
library(tidyverse)
library(sf)
library(maptools)
library(cowplot)
library(MetBrewer)

setwd(here::here("data/Shapefiles/reserve/"))

res <- sf::st_read("Reserve_polygon.shp") %>% 
  dplyr::select(Id, geometry) %>% 
  mutate(region = ifelse(Id == 0, "Talek", 
                         ifelse(Id == 2, "Mara Triangle", "River")))

setwd(here::here("data/Shapefiles/DS"))
ds <-sf::st_read("DS_10kmpersite.shp") %>% 
  add_column(data = "Distance sampling") %>% 
  dplyr::select(data, geometry)

setwd(here::here("data/Shapefiles/Transects"))
tc <- sf::st_read("Transects.shp") %>% 
  add_column(data = "Counts") %>% 
  dplyr::select(data, geometry)

surveys <- rbind(ds, tc)

fillpal <- MetBrewer::MetPalettes$Hiroshige[[1]][c(4, 6, 6)]
linepal <- MetBrewer::MetPalettes$Hiroshige[[1]][c(1, 10)]

talek <- tibble(
  x = 35.214983,
  y = -1.441084) %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = 4326) %>% 
  st_transform(st_crs(res)) %>% 
  add_column(name = "Talek town") %>% 
  mutate(x = st_coordinates(.)[,1] + 6000,
         y = st_coordinates(.)[,2] + 1000)

mt_label <- tibble(
  x = c(34.909227),
  y = c(-1.446516), 
  name = c("Mara Triangle")) %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = 4326) %>% 
  st_transform(st_crs(res)) %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])

tr_label <- tibble(
  x = c(35.161652),
  y = c(-1.546157), 
  name = c("Talek region")) %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = 4326) %>% 
  st_transform(st_crs(res)) %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])

tt_label <- tibble(
  x = c(35.291453),
  y = c(-1.425038), 
  name = c("Talek town")
) %>% 
  st_as_sf(coords = c("x", "y"), 
           crs = 4326) %>% 
  st_transform(st_crs(res)) %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2])

res_plot <-
  ggplot(  ) + 
  geom_sf(data = res, 
          aes( geometry = geometry, fill = region),
          color = NA ) +
  geom_sf(data = talek, 
          aes(geometry = geometry), 
          size = 4) +
  geom_sf( data = surveys, 
           aes( geometry = geometry, 
                color = data,
                lwd = data) ) +
  geom_text(data = mt_label, 
            aes(x = x, 
                y = y, 
                label = name),
            angle = 330,
            color = MetBrewer::MetPalettes$Hiroshige[[1]][2]) +
  geom_text(data = tr_label, 
            aes(x = x, 
                y = y, 
                label = name), 
            color = MetBrewer::MetPalettes$Hiroshige[[1]][8]) +
  geom_text(data = tt_label,
            aes(x = x, 
                y = y, 
                label = name), 
            color = "black",
            size = 3.1) +
  scale_size_manual(values = c(1.5, 0.5)) +
  scale_color_manual(values = linepal) + 
  theme_minimal() +
  scale_fill_manual(
    values = fillpal,
    limits = c("Mara Triangle", "Talek"))  +  
  theme(
    legend.position = "none", 
    axis.text = element_text(color = "black", 
                             size = 9),
    legend.key = element_rect(color = NA, fill = NA),
    legend.key.size = unit(0.5, "cm"),
    axis.title = element_blank() )

data("wrld_simpl")

wrld_sf <- st_as_sf(wrld_simpl)

ken <- wrld_sf %>%
  filter(NAME == "Kenya") %>% 
  st_transform(st_crs(res))

res_bbox <- st_bbox(res) %>% 
  st_as_sfc(.) %>% 
  st_buffer(., 30000)

ken_label <- st_centroid(ken) %>% 
  mutate(x = st_coordinates(.)[,1],
         y = st_coordinates(.)[,2]) %>% 
  st_drop_geometry() %>% 
  add_column(name = "Kenya") %>% 
  ungroup() %>% 
  dplyr::select(x, y, name)

inset <-
  ggplot() + 
  geom_sf(data = ken, aes(geometry = geometry),
          size = 0.1) + 
  geom_sf(data = res, aes(geometry = geometry), 
          color = "black", 
          fill = "black") +
  geom_sf(data = res_bbox, aes(geometry = geometry), 
          color = "red", 
          fill = NA,
          size = 0.75) +
  geom_text(data = ken_label, aes(x = x, y = y, 
                                  label = name),
            color = "black") +
  theme_void() +
  theme(plot.background = element_rect(fill = "white",
                                       color = "gray20"))

for_legend <- ggplot( data = 
                        tibble( x = c(1, 2), 
                                xend = c(2, 3), 
                                y = c(1, 2), 
                                yend = c( 2, 3),
                                data = c("Counts", "Distance sampling")
                        )) + 
  geom_segment(aes(x = x, 
                   xend = xend, 
                   y = y, 
                   yend = yend,
                   color = data, 
                   lwd = data)) +
  scale_size_manual(values = c(2, 1)) +
  scale_color_manual(values = linepal) +
  theme_minimal() +
  theme(legend.title = element_blank(),
        legend.background = element_rect(fill = "white", color = NA))

leg <- cowplot::get_legend(for_legend)

ggdraw() + 
  draw_plot(res_plot) + 
  draw_plot(inset,
            x = 0.7,
            y = 0.67,
            width = 0.3,
            height = 0.32) + 
  draw_plot(leg,
            x = -0.15, 
            y = -0.25)

setwd(here::here("figures"))
ggsave(
  "study_area_map_v02.png",
  width = 4.5, 
  height = 3, 
  units = "in", 
  dpi = 300)