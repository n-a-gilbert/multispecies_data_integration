# 24 October 2023
# Create convergence summary table for ICM and alternative models

library(here)
library(tidyverse)
library(officer)
library(flextable)
library(magrittr)

setwd(here::here("results"))

ic <- readr::read_csv("ic_simulation.csv") |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "ic")

dc <- readr::read_csv("dc_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "dc")

cc <- readr::read_csv("cc_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "cc")

isr <- readr::read_csv("isr_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "isr")

isc <- readr::read_csv("isc_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "isc")

dsr <- readr::read_csv("dsr_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "dsr")

dsc <- readr::read_csv("dsc_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "dsc")

csr <- readr::read_csv("csr_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "csr")

csc <- readr::read_csv("csc_simulation.csv") |> 
  dplyr::filter(!is.na(mean)) |> 
  tidyr::separate(param, into = c("param", "junk"), sep = "\\[") |> 
  dplyr::select(simrep, param, sp, truth:Rhat) |> 
  tibble::add_column(model = "csc")

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(isr) |> 
  dplyr::full_join(isc) |> 
  dplyr::full_join(dsr) |> 
  dplyr::full_join(dsc) |> 
  dplyr::full_join(csr) |> 
  dplyr::full_join(csc)

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
