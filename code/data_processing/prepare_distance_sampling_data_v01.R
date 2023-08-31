library(here)
library(janitor)
library(tidyverse)
library(sf)
library(lubridate)

setwd(here::here("data"))

ds_shape <- sf::st_read(dsn = "./Shapefiles/DS", layer = "DS_10kmpersite") 

min_coords <- sf::st_coordinates(sf::st_transform(ds_shape, 4326)) |> 
  tibble::as_tibble() |> 
  dplyr::group_by(L2) |> 
  dplyr::summarise(x = min(X), 
                   y = min(Y))

ds_offset <- ds_shape |>
  tibble::add_column(minX = min_coords$x) |> 
  dplyr::mutate(site = 1:length(unique(X10kmsiteID)) ,
                area = ( as.numeric(sf::st_length(geometry)) * 1000 * 2 ) / 1E6,
                region = ifelse(minX > 35.1, 1, 0)) |> 
  sf::st_drop_geometry() |> 
  dplyr::select(site, area, region)

ds <- readr::read_csv("Herbivore Utilization Complete.csv") |> 
  janitor::clean_names() |> 
  dplyr::filter(animal %in% c(
    "Buffalo",
    "Eland",
    "Elephant",
    "Giraffe",
    "Grants",
    "Hartebeest",
    "Impala",
    "Thomsons",
    "Topi",
    "Warthog",
    "Waterbuck")) |> 
  dplyr::filter(! is.na(count)) |> 
  dplyr::filter(! count == 0) |>
  # new base pipe doesn't have placeholder: https://tinyurl.com/2p8zarxr
  (\(x) sf::st_as_sf(x, 
                    coords = c("adj_easting", "adj_northing"),
                    crs = sf::st_crs(ds_shape)))()
 
ds_matrix <- sf::st_distance(ds, ds_shape)

# assign transect to each observation 
ds$site <- apply(ds_matrix, 1, which.min)

# add corrected distances
ds$dst <- apply(ds_matrix, 1, min)

ds <- ds |> 
  dplyr::filter(dst <= 1000)

#-----------------------#
# Assign distance classes #
#-----------------------#

#Number of observations for distance sampling
nobs <- NULL

nobs[1] <- dim(ds)[1]

#ID for distance class
di <- seq(0, 1000, 25)

#Distance class
dclass <- rep(NA, nobs[1])

#Number of distance classes
nG <- length(di) - 1

#Minimum distance to assigned transect
dst <- ds$dst

for(i in 1:nobs[1]){
  for(k in 1:nG){
    if(di[k] < dst[i] && dst[i] <= di[k+1])
      dclass[i] <- k
    
  }
}

ds$dclass <- dclass

#Replicate
ds$reps[ds$territory == "North"] <- dplyr::filter(ds, territory == "North") |> 
  dplyr::group_by(year, month )|> 
  dplyr::group_indices()

ds$reps[ds$territory == "South"] <- dplyr::filter(ds, territory =="South") |> 
  dplyr::group_by(year, month) |> 
  dplyr::group_indices()

ds$reps[ds$territory == "West"] <- dplyr::filter(ds, territory == "West") |> 
  dplyr::group_by(year, month) |> 
  dplyr::group_indices()

final <- ds |> 
  sf::st_drop_geometry() |> 
  dplyr::mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) |> 
  dplyr::select(date, animal, site, rep = reps, count, dclass) |>
  dplyr::arrange(animal, site, rep) |> 
  dplyr::mutate(animal = factor(animal)) |> 
  dplyr::mutate(spec = as.numeric(animal)) |> 
  dplyr::select(date,
                sp_name = animal, 
                sp = spec, 
                site, 
                rep, 
                gs = count, 
                dclass) |> 
  dplyr::group_by(sp, site, rep) |> 
  dplyr::mutate(ng = n())

site_key <- final |> 
  dplyr::ungroup() |> 
  dplyr::select(site, rep, date) |> 
  dplyr::distinct() |> 
  dplyr::group_by(site, rep) |> 
  dplyr::summarise(date = min(date))

sp_key <- final |> 
  dplyr::ungroup() |> 
  dplyr::select(sp, sp_name) |> 
  dplyr::distinct() |> 
  dplyr::mutate(sp_name = tolower(sp_name))

#Width of distance classes
v <- 25 # meters

# Transect half-width
b <- 1000 # meters

#Distance class midpoint ID
mdpt <- seq( v/2, b, v)

#area of transects (m^2)
area <- as.numeric( sf::st_length(ds_shape)*1000*2)

#set baseline unit as 1 km^2
offset <- area / 1E6

final2 <- expand.grid(
  sp = sort(unique(final$sp)),
  site = sort(unique(final$site)),
  rep = sort(unique(final$rep))) |> 
  tibble::as_tibble() |> 
  dplyr::full_join(
    dplyr::select(
      dplyr::ungroup(final),
      sp:ng)
    ) |> 
  dplyr::full_join(ds_offset) |> 
  dplyr::full_join(site_key) |> 
  dplyr::full_join(sp_key) |> 
  dplyr::mutate(ng = ifelse(is.na(ng) & (!is.na(date)), 0, ng))

save(
  final2, 
  v, 
  b, 
  mdpt, 
  file = "distance_sampling_data_v01.RData")