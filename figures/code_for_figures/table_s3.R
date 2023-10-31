library(here)
library(tidyverse)
library(MetBrewer)
library(officer)
library(flextable)
library(magrittr)

setwd(here::here("results"))

load("ic.RData")
ic <- ic |> 
  dplyr::filter( (nobs > 1 & ndist > 1) | is.na(nobs)) |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(nobs, na.rm = TRUE),
                max_num_obs = max(nobs, na.rm = TRUE)) |> 
  dplyr::filter(nobs == min_num_obs | nobs == max_num_obs | is.na(nobs)) |> 
  dplyr::mutate( type = ifelse(nobs == min_num_obs, "rare", 
                               ifelse(nobs == max_num_obs, "common", NA)),
                 nsp = length(unique(sp))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter( sp == first_sp | is.na(first_sp)) |> 
  tibble::add_column(model = "ic") |> 
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`) 

load("dc.RData")
dc <- dc |>   
  dplyr::filter( (nobs > 1 & ndist > 1) | is.na(nobs)) |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(nobs, na.rm = TRUE),
                max_num_obs = max(nobs, na.rm = TRUE)) |> 
  dplyr::filter(nobs == min_num_obs | nobs == max_num_obs | is.na(nobs)) |> 
  dplyr::mutate( type = ifelse(nobs == min_num_obs, "rare", 
                               ifelse(nobs == max_num_obs, "common", NA))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter( sp == first_sp | is.na(first_sp)) |> 
  dplyr::filter(!is.na(mean)) |> 
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`) 

load("cc.RData")
cc <- cc |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(nobs, na.rm = TRUE), 
                max_num_obs = max(nobs, na.rm = TRUE)) |> 
  dplyr::filter(nobs == min_num_obs | nobs == max_num_obs | is.na(nobs)) |> 
  dplyr::mutate(type = ifelse(nobs == min_num_obs, "rare", 
                              ifelse(nobs == max_num_obs, "common", NA))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter(sp == first_sp | is.na(first_sp)) |> 
  dplyr::filter(!is.na(mean)) |> 
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`)

load("is.RData")
is <- is |> 
  dplyr::mutate(type = ifelse(model == "isr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean))

load("ds.RData")
ds <- ds |> 
  dplyr::mutate(type = ifelse(model == "dsr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

load("cs.RData")
cs <- cs  |> 
  dplyr::mutate(type = ifelse(model == "csr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(is) |> 
  dplyr::full_join(ds) |> 
  dplyr::full_join(cs) 

com_params <- all |> 
  dplyr::mutate(diff = mean - truth) |> 
  dplyr::mutate(model = toupper(model)) |> 
  dplyr::group_by(model, type, param) |> 
  dplyr::summarise(ab = mean(diff)) |> 
  dplyr::filter(is.na(type)) |> 
  dplyr::ungroup() |> 
  dplyr::select(-type) |> 
  dplyr::full_join( tibble::tibble(
    model = c("IC", "IC", "DC", "DC", "CC", "CC"), 
    type = c("rare", "common", "rare", "common", "rare", "common"))) |> 
  dplyr::select(model, type, param, ab)

ab_table <- all |> 
  dplyr::mutate(diff = mean - truth) |> 
  dplyr::mutate(model = toupper(model)) |> 
  dplyr::group_by(model, type, param) |> 
  dplyr::summarise(ab = mean(diff)) |> 
  dplyr::filter(!is.na(type)) |> 
  dplyr::full_join(com_params) |> 
  dplyr::mutate(type = toupper(str_sub(type, 1, 1))) |> 
  dplyr::mutate(model = ifelse(model == "ISC", "IS", 
                               ifelse(model == "ISR", "IS", 
                                      ifelse(model == "DSR", "DS", 
                                             ifelse(model == "DSC", "DS", 
                                                    ifelse(model == "CSR", "CS", 
                                                           ifelse(model == "CSC", "CS", model))))))) |> 
  dplyr::mutate(model = paste0(model, type)) |> 
  dplyr::mutate(model = factor(model, levels = c(
    "ICR", "ICC", 
    "DCR", "DCC", 
    "CCR", "CCC", 
    "ISR", "ISC", 
    "DSR", "DSC", 
    "CSR", "CSC"
  ))) |> 
  dplyr::arrange(model) |>
  dplyr::ungroup() |> 
  dplyr::select(-type) |> 
  dplyr::mutate( ab = paste0( sprintf("%.2f", round( ab, 2)))) |> 
  dplyr::mutate(param = factor(param, levels = c("N_DS", "N_TC",
                                                 "mu_gamma0", "sd_gamma0", "gamma0_ds",
                                                 "mu_gamma0_c", "sd_gamma0_c", "gamma0_c",
                                                 "mu_alpha0", "sd_alpha0",  "alpha0", 
                                                 "mu_alpha1", "sd_alpha1", "alpha1"))) |> 
  dplyr::arrange(model, param) |> 
  tidyr::pivot_wider(names_from = model, values_from = ab) |> 
  dplyr::rename(Parameter = param)

flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable( data = ab_table, cwidth = 0.7)  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)
