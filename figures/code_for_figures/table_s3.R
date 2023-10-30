library(here)
library(tidyverse)
library(MetBrewer)
library(officer)
library(flextable)
library(magrittr)

setwd(here::here("results"))

ic <- readr::read_csv("ic_simulation.csv") |> 
  dplyr::filter( (totDS_obs > 1 & ndistances > 1) | is.na(totDS_obs)) |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(totDS_obs, na.rm = TRUE),
               max_num_obs = max(totDS_obs, na.rm = TRUE)) |> 
  dplyr::filter(totDS_obs == min_num_obs | totDS_obs == max_num_obs | is.na(totDS_obs)) |> 
  dplyr::mutate( type = ifelse(totDS_obs == min_num_obs, "rare", 
                               ifelse(totDS_obs == max_num_obs, "common", NA)),
                 nsp = length(unique(sp))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter( sp == first_sp | is.na(first_sp)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  tibble::add_column(model = "ic")

dc <- readr::read_csv("dc_simulation.csv") |>  
  dplyr::filter( (totDS_obs > 1 & ndistances > 1) | is.na(totDS_obs)) |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(totDS_obs, na.rm = TRUE),
                max_num_obs = max(totDS_obs, na.rm = TRUE)) |> 
  dplyr::filter(totDS_obs == min_num_obs | totDS_obs == max_num_obs | is.na(totDS_obs)) |> 
  dplyr::mutate( type = ifelse(totDS_obs == min_num_obs, "rare", 
                               ifelse(totDS_obs == max_num_obs, "common", NA))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter( sp == first_sp | is.na(first_sp)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "dc")

cc <- readr::read_csv("cc_simulation.csv") |> 
  dplyr::group_by(simrep) |> 
  dplyr::mutate(min_num_obs = min(num_obs, na.rm = TRUE), 
                max_num_obs = max(num_obs, na.rm = TRUE)) |> 
  dplyr::filter(num_obs == min_num_obs | num_obs == max_num_obs | is.na(num_obs)) |> 
  dplyr::mutate(type = ifelse(num_obs == min_num_obs, "rare", 
                              ifelse(num_obs == max_num_obs, "common", NA))) |> 
  dplyr::group_by(simrep, type) |> 
  dplyr::mutate(first_sp = first(sp)) |> 
  dplyr::filter(sp == first_sp | is.na(first_sp)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "cc")

isr <- readr::read_csv("isr_simulation.csv") |> 
  tibble::add_column(type = "rare") |>
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "is")

isc <- readr::read_csv("isc_simulation.csv") |> 
  tibble::add_column(type = "common") |>
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "is")

dsr <- readr::read_csv("dsr_simulation.csv") |> 
  tibble::add_column(type = "rare") |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "ds")

dsc <- readr::read_csv("dsc_simulation.csv") |> 
  tibble::add_column(type = "common") |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "ds")

csr <- readr::read_csv("csr_simulation.csv") |> 
  tibble::add_column(type = "rare") |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "cs")

csc <- readr::read_csv("csc_simulation.csv") |> 
  tibble::add_column(type = "common") |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(param, type, truth, mean, sd, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) |> 
  tibble::add_column(model = "cs")

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(isr) |> 
  dplyr::full_join(isc) |> 
  dplyr::full_join(dsr) |> 
  dplyr::full_join(dsc) |> 
  dplyr::full_join(csr) |> 
  dplyr::full_join(csc)

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

rb_table <- all |> 
  dplyr::mutate(diff = mean - truth) |> 
  dplyr::mutate(model = toupper(model)) |> 
  dplyr::group_by(model, type, param) |> 
  dplyr::summarise(ab = mean(diff)) |> 
  dplyr::filter(!is.na(type)) |> 
  dplyr::full_join(com_params) |> 
  dplyr::mutate(type = toupper(str_sub(type, 1, 1))) |> 
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
  dplyr::mutate( ab = sprintf("%.2f", round( ab, 2))) |> 
  dplyr::mutate(param = factor(param, levels = c("N_DS", "N_TC",
                                                 "mu_gamma0", "sd_gamma0", "gamma0_ds",
                                                 "mu_gamma0_c", "sd_gamma0_c", "gamma0_c",
                                                 "mu_alpha0", "sd_alpha0",  "alpha0", 
                                                 "mu_alpha1", "sd_alpha1", "alpha1"))) |> 
  dplyr::arrange(model, param) |> 
  tidyr::pivot_wider(names_from = model, values_from = ab) |> 
  dplyr::rename(Parameter = param)

flextable::set_flextable_defaults(font.size = 10)
ft <- flextable::flextable( data = rb_table, cwidth = 0.7)  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)
