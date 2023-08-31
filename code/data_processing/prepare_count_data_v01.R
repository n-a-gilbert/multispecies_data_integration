library(here)
library(tidyverse)
library(sf)
library(lubridate)

setwd(here::here("data"))

tc <- readr::read_csv("tblPreyCensus_2012to2014.csv")

tc_shape <- sf::st_read(dsn = "./Shapefiles/Transects", layer = "Transects")

min_coords <- sf::st_coordinates(sf::st_transform(tc_shape, 4326)) |>  
  tibble::as_tibble() |>  
  dplyr::group_by(L1) |>  
  dplyr::summarise(x = min(X), 
                   y = min(Y))

transect_key <- tibble::tribble(
  ~tc_name, ~tc_shape_name, 
  "W3",  "Prey 3",
  "RSP", "RSP",
  "WLOW", "Low Road",
  "WHIGH", "High Road",
  "S1", "South I",
  "S2", "South II",
  "SST", "SST",
  "North", "North",
  "HZT", "HZT")

offset <- tc_shape |>  
  tibble::add_column(minX = min_coords$x) |>  
  dplyr::mutate(area = ( as.numeric(st_length(geometry)) * 100 * 2 ) / 1E6,
                region = ifelse(minX > 35.1, 1, 0)) |>  
  sf::st_drop_geometry() |>  
  dplyr::select(tc_shape_name  = Transect, area, region) |>  
  dplyr::full_join(transect_key) |>  
  dplyr::select(transect = tc_name, area, region)

transect_data <- tc |>  
  dplyr::select(transect, date, baboon:waterbuck) |>  
  tidyr::pivot_longer(baboon:waterbuck, names_to = "animal", values_to = "count") |>  
  dplyr::filter(animal %in% c(
    "buffalo",
    "eland",
    "elephant",
    "giraffe",
    "grants",
    "hartebeest",
    "impala",
    "thomsons",
    "topi",
    "warthog",
    "waterbuck")) |>  
  dplyr::filter( transect %in% c(
    "WLOW",
    "WHIGH",
    "W3",
    "North",
    "S1",
    "S2",
    "RSP",
    "SST",
    "HZT")) |>  
  tidyr::separate(col = date, into = c("month", "day", "year"), sep = "/") |>  
  dplyr::mutate(across(month:year, as.numeric)) |>  
  dplyr::group_by(animal, transect, month, year) |>  
  dplyr::mutate(min_day = min(day), 
                max_day = max(day), 
                n_day = n_distinct(day)) |>  
  dplyr::filter(!n_day == 3 & max_day == day) |>  
  dplyr::mutate(k = ifelse(n_day == 1 & day < 16, 1,
                           ifelse(n_day == 2 & day == min_day, 1,
                                  ifelse(n_day == 3 & day == min_day, 1,0)))) |> 
  dplyr::ungroup() |>  
  dplyr::mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) |>  
  dplyr::group_by(animal, transect) |>  
  dplyr::arrange(transect, animal, date) |>  
  dplyr::mutate(rep = row_number()) |>  
  dplyr::arrange(animal, transect, rep) |>  
  dplyr::mutate(animal = factor(animal),
                transect = factor(transect)) |>  
  dplyr::mutate(spec = as.numeric(animal),
                site = as.numeric(transect)) |>  
  dplyr::ungroup() |>  
  dplyr::select(transect, sp_name = animal, date, sp = spec, site, rep, count) |>  
  dplyr::full_join(offset)

setwd(here::here("data"))
save(
  transect_data, 
  file = "count_data_v01.RData")