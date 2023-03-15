library(here)
library(tidyverse)
library(sf)

setwd(here::here("data"))

tc <- read_csv("tblPreyCensus_2012to2014.csv")

tc_shape <- st_read(dsn = "./Shapefiles/Transects", layer = "Transects")

min_coords <- st_coordinates(st_transform(tc_shape, 4326)) %>% 
  as_tibble() %>% 
  group_by(L1) %>% 
  summarise(x = min(X), 
            y = min(Y))

transect_key <- tribble(
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

offset <- tc_shape %>% 
  add_column(minX = min_coords$x) %>% 
  mutate(area = ( as.numeric(st_length(.)) * 100 * 2 ) / 1E6,
         region = ifelse(minX > 35.1, 1, 0)) %>% 
  st_drop_geometry() %>% 
  dplyr::select(tc_shape_name  = Transect, area, region) %>% 
  full_join(transect_key) %>% 
  dplyr::select(transect = tc_name, area, region)

transect_data <- tc %>% 
  dplyr::select(transect, date, baboon:waterbuck) %>% 
  pivot_longer(baboon:waterbuck, names_to = "animal", values_to = "count") %>% 
  filter(animal %in% c(
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
    "waterbuck")) %>% 
  filter( transect %in% c(
    "WLOW",
    "WHIGH",
    "W3",
    "North",
    "S1",
    "S2",
    "RSP",
    "SST",
    "HZT")) %>% 
  separate(col = date, into = c("month", "day", "year"), sep = "/") %>% 
  mutate(across(month:year, as.numeric)) %>% 
  group_by(animal, transect, month, year) %>% 
  mutate(min_day = min(day), 
         max_day = max(day), 
         n_day = n_distinct(day)) %>% 
  filter(!n_day == 3 & max_day == day) %>% 
  mutate(k = ifelse(n_day == 1 & day < 16, 1,
                    ifelse(n_day == 2 & day == min_day, 1,
                           ifelse(n_day == 3 & day == min_day, 1,0)))) %>%
  ungroup() %>% 
  mutate(date = lubridate::ymd(paste(year, month, day, sep = "-"))) %>% 
  group_by(animal, transect) %>% 
  arrange(transect, animal, date) %>% 
  mutate(rep = row_number()) %>% 
  arrange(animal, transect, rep) %>% 
  mutate(animal = factor(animal),
         transect = factor(transect)) %>% 
  mutate(spec = as.numeric(animal),
         site = as.numeric(transect)) %>% 
  ungroup() %>% 
  dplyr::select(transect, sp_name = animal, date, sp = spec, site, rep, count) %>% 
  full_join(offset)

setwd(here::here("data"))
save(
  transect_data, 
  file = "count_data_v01.RData")
