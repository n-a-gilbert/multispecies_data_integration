library(here)
library(tidyverse)
library(MetBrewer)

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
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) 

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
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) 

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
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`)

load("is.RData")
is <- is |> 
  dplyr::mutate(type = ifelse(model == "isr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

load("ds.RData")
ds <- ds |> 
  dplyr::mutate(type = ifelse(model == "dsr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

load("cs.RData")
cs <- cs  |> 
  dplyr::mutate(type = ifelse(model == "csr", "rare", "common")) |>
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(is) |> 
  dplyr::full_join(ds) |> 
  dplyr::full_join(cs) 

all |> 
  dplyr::group_by( model, type, param) |>
  dplyr::summarise( rb = 100*sum( mean - truth) / sum(abs(truth))) |>
  dplyr::filter(param %in% c(
    "gamma0_ds",
    "gamma0_c", 
    "alpha0", 
    "alpha1")) |> 
  dplyr::mutate( model_name = ifelse( model == "csr", "Single species, counts",
                                      ifelse( model == "csc", "Single species, counts",
                                              ifelse(model == "dsr", "Single species, distance sampling",
                                                     ifelse(model == "dsc", "Single species, distance sampling",
                                                            ifelse(model == "isc", "Single species, integrated", 
                                                                   ifelse(model == "isr", "Single species, integrated",
                                                                          ifelse( model == "cc", "Community, counts", 
                                                                                  ifelse(model == "dc", "Community, distance sampling", 
                                                                                         ifelse(model == "ic", "Community, integrated", NA)))))))))) |> 
  
  
  
  dplyr::mutate( model_name = factor(model_name, 
                                     levels = c(
                                       "Single species, counts",
                                       "Single species, distance sampling",
                                       "Single species, integrated", 
                                       "Community, counts",
                                       "Community, distance sampling",
                                       "Community, integrated"))) |> 
  dplyr::mutate(type = paste(type, "species")) |> 
  dplyr::mutate(type = str_to_sentence(type)) |>
  dplyr::mutate(type = factor(type, levels = c("Rare species", "Common species"))) |> 
  dplyr::mutate(param = ifelse(param == "alpha0", "Abundance intercept", 
                               ifelse(param == "alpha1", "Abundance covariate coefficient",
                                      ifelse(param == "gamma0_ds", "Detection (distance sampling)", 
                                             "Detection (counts)")))) |> 
  dplyr::mutate( param = factor(param, levels = c(
    "Abundance intercept", 
    "Abundance covariate coefficient",
    "Detection (distance sampling)",
    "Detection (counts)"))) |> 
  dplyr::group_by(type, param) |> 
  dplyr::mutate(
    rank_rb = rank( abs(rb)),
    rb_lab = paste0(  sprintf("%.2f", round(rb, 2)), "%")) |> 
  dplyr::mutate(rb_lab_x = ifelse(rank_rb > 4, rank_rb - 1.25, rank_rb + 1.25)) |> 
  
  ggplot2::ggplot(aes(x = rank_rb, y = model_name, color = type)) +
  ggplot2::facet_grid(type~param, scales = "free_x") +
  ggplot2::geom_point(size = 3) + 
  ggplot2::geom_text(aes(x = rb_lab_x, label = rb_lab), size = 3, color = "black") + 
  ggplot2::xlab( "Rank relative bias") +
  ggplot2::scale_x_continuous(limits = c(0, 7),
                              breaks = c(1:6)) +
  ggplot2::scale_color_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  ggplot2::theme(panel.grid.minor = element_blank(),
                 axis.line = element_line(color = "black", size = 0.2),
                 axis.ticks = element_line(color = "black", size = 0.2),
                 axis.text.x = element_text(size = 9, color = "black"),
                 axis.title.x = element_text(size = 10, color = "black"),
                 legend.text = element_text(size = 10, color = "black"),
                 strip.text = element_text(size = 9, color = "black"), 
                 legend.title = element_blank(), 
                 legend.position = "bottom",
                 axis.title.y = element_blank(),
                 axis.text.y = element_text( size = c(9,9,9,9,9,12), color = "black",
                                             face = c("plain", "plain", "plain", "plain", "plain", "bold")))

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_s6.png", 
  width = 10, 
  height = 5, 
  units = "in", 
  dpi = 300
)