library(here)
library(janitor)
library(tidyverse)
library(sf)

setwd(here::here("data"))

ds_shape <- st_read(dsn = "./Shapefiles/DS", layer = "DS_10kmpersite") #CB

min_coords <- st_coordinates(st_transform(ds_shape, 4326)) %>% 
  as_tibble() %>% 
  group_by(L2) %>% 
  summarise(x = min(X), 
            y = min(Y))

ds_offset <- ds_shape %>%
  add_column(minX = min_coords$x) %>% 
  mutate(site = 1:nrow(.) ,
         area = ( as.numeric(st_length(.)) * 1000 * 2 ) / 1E6,
         region = ifelse(minX > 35.1, 1, 0)) %>% 
  st_drop_geometry() %>% 
  dplyr::select(site, area, region)

ds <- read_csv("Herbivore Utilization Complete.csv") %>% 
  janitor::clean_names() %>% 
  filter(animal %in% c(
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
    "Waterbuck")) %>% 
  filter(! is.na(count)) %>% 
  filter(! count == 0) %>%
  st_as_sf(.,
           coords = c("adj_easting", "adj_northing"),
           crs = st_crs(ds_shape))

ds_matrix <- st_distance(ds, ds_shape)

# assign transect to each observation 
ds$site <- apply(ds_matrix, 1, which.min)

# add corrected distances
ds$dst <- apply(ds_matrix, 1, min)

ds <- ds %>% filter(dst <= 1000)

#-----------------------#
#Assign distance classes#
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
    if(di[k] < dst[i] && dst[i] <= di[k+1]) #why the k+1 argument? 
      dclass[i] <- k
    
  }
}

ds$dclass <- dclass

#Replicate
ds$reps[ds$territory == "North"] <- filter(ds, territory == "North") %>% 
  group_by(year, month )%>% 
  group_indices()

ds$reps[ds$territory == "South"] <- filter(ds, territory =="South") %>% 
  group_by(year, month) %>% 
  group_indices()

ds$reps[ds$territory == "West"] <- filter(ds, territory == "West") %>% 
  group_by(year, month) %>% 
  group_indices()

final <- ds %>% 
  st_drop_geometry() %>% 
  mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) %>% 
  dplyr::select(date, animal, site, rep = reps, count, dclass) %>%
  arrange(animal, site, rep) %>% 
  mutate(animal = factor(animal)) %>% 
  mutate(spec = as.numeric(animal)) %>% 
  dplyr::select(date,
                sp_name = animal, 
                sp = spec, 
                site, 
                rep, 
                gs = count, 
                dclass) %>% 
  group_by(sp, site, rep) %>% 
  mutate(ng = n())

site_key <- final %>% 
  ungroup() %>% 
  dplyr::select(site, rep, date) %>% 
  distinct(.) %>% 
  group_by(site, rep) %>% 
  summarise(date = min(date))

sp_key <- final %>% 
  ungroup() %>% 
  dplyr::select(sp, sp_name) %>% 
  distinct() %>% 
  mutate(sp_name = tolower(sp_name))

#Width of distance classes
v <- 25 # meters

# Transect half-width
b <- 1000 # meters

#Distance class midpoint ID
mdpt <- seq( v/2, b, v)

#area of transects (m^2)
area <- as.numeric( st_length(ds_shape)*1000*2)

#set baseline unit as 1 km^2
offset <- area / 1E6

final2 <- expand.grid(
  sp = sort(unique(final$sp)),
  site = sort(unique(final$site)),
  rep = sort(unique(final$rep))) %>% 
  as_tibble() %>% 
  full_join(dplyr::select(ungroup(final), sp:ng)) %>% 
  full_join(ds_offset) %>% 
  full_join(site_key) %>% 
  full_join(sp_key) %>% 
  mutate(ng = ifelse(is.na(ng) & (!is.na(date)), 0, ng))

save(
  final2, 
  v, 
  b, 
  mdpt, 
  file = "ds_data_ng_v06.RData")