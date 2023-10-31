# 24 October 2023
# Create convergence summary table for ICM and alternative models
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
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat) 

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
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat) 

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
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat)

load("is.RData")
is <- is |> 
  dplyr::mutate(type = ifelse(model == "isr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat) |> 
  dplyr::filter(!is.na(mean))

load("ds.RData")
ds <- ds |> 
  dplyr::mutate(type = ifelse(model == "dsr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat) |> 
  dplyr::filter(!is.na(mean)) 

load("cs.RData")
cs <- cs  |> 
  dplyr::mutate(type = ifelse(model == "csr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, sd, `2.5%`, `97.5%`, Rhat) |> 
  dplyr::filter(!is.na(mean)) 

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(is) |> 
  dplyr::full_join(ds) |> 
  dplyr::full_join(cs) 

convergence_table <- all |> 
  dplyr::mutate(Rhat = ifelse(is.na(Rhat), 10, Rhat)) |> # give Rhats with NA (unconverged) an arbitrarily large number
  dplyr::group_by(model, simrep) |> 
  dplyr::summarise(max_rhat = max(Rhat)) |> 
  dplyr::summarise( n = dplyr::n(),
                    `Proportion unconverged replicates` = sum(max_rhat > 1.1) / n) |> 
  dplyr::select(-n) |> 
  dplyr::full_join(
    all |> 
      dplyr::mutate(Rhat = ifelse(is.na(Rhat), 10, Rhat)) |> # give Rhats with NA (unconverged) an arbitrarily large number
      dplyr::group_by( model ) |> 
      dplyr::summarise(
        n = dplyr::n(),
        `Proportion unconverged parameters` = sum(Rhat > 1.1) / n) |> 
      dplyr::select(-n)
  ) |> 
  dplyr::mutate(model = toupper(model)) |> 
  dplyr::mutate(model = factor(model, levels = c("IC", "DC", "CC", "ISR", "ISC", "DSR", "DSC", "CSR", "CSC"))) |> 
  dplyr::arrange(model) |> 
  tidyr::pivot_longer(`Proportion unconverged replicates`:`Proportion unconverged parameters`, names_to = "Metric", values_to = "value") |> 
  tidyr::pivot_wider(names_from = model, values_from = value) |>
  dplyr::mutate(across(IC:CSC, function(x) ifelse( x < 0.0005 & x > 0, 
                                                   format(x, scientific = TRUE, digits = 2), 
                                                   sprintf("%.3f", round(x, digits = 3)))))



ft <- flextable::flextable( data = convergence_table )  

tmp <- tempfile(fileext = ".docx")

officer::read_docx() |> 
  flextable::body_add_flextable(ft) |> 
  print(target = tmp)

utils::browseURL(tmp)
