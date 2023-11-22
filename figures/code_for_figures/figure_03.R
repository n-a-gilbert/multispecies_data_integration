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
  dplyr::mutate(type = ifelse(model == "isr", "rare",
                              ifelse(model == "isc", "common", NA))) |> 
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

load("ds.RData")
ds <- ds |> 
  dplyr::mutate(type = ifelse(model == "dsr", "rare",
                              ifelse(model == "dsc", "common", NA))) |> 
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

load("cs.RData")
cs <- cs  |> 
  dplyr::mutate(type = ifelse(model == "csr", "rare",
                              ifelse(model == "csc", "common", NA))) |>
  dplyr::select(model, simrep, param, type, truth, mean, `2.5%`, `97.5%`) |> 
  dplyr::filter(!is.na(mean)) 

all <- dplyr::full_join(ic, dc) |> 
  dplyr::full_join(cc) |> 
  dplyr::full_join(is) |> 
  dplyr::full_join(ds) |> 
  dplyr::full_join(cs) 

all |> 
  dplyr::filter(param %in% c(
    "gamma0_ds", 
    "gamma0_c", 
    "alpha0", 
    "alpha1")) |> 
  dplyr::mutate(diff = mean - truth) |> 
  dplyr::mutate(type = paste(type, "species")) |> 
  dplyr::mutate(type = str_to_sentence(type)) |> 
  dplyr::mutate(param = ifelse(param == "alpha0", "Abundance intercept", 
                               ifelse(param == "alpha1", "Abundance covariate coefficient",
                                      ifelse(param == "gamma0_ds", "Detection intercept (distance sampling)", 
                                             "Detection intercept (counts)")))) |> 
  dplyr::mutate( param = factor(param, levels = c(
    "Abundance intercept", 
    "Abundance covariate coefficient",
    "Detection intercept (distance sampling)",
    "Detection intercept (counts)"))) |> 
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
  
  
  ggplot2::ggplot(aes(x = diff, y = model_name, fill = type)) +
  ggplot2::geom_vline(xintercept = 0,
                      color = MetPalettes$Hiroshige[[1]][c(1)],
                      linetype = "dashed") +
  ggplot2::facet_wrap(~param, scales = "free_x") +
  ggplot2::geom_boxplot( outlier.alpha = 0.2,
                         outlier.size = 0.75, 
                         size = 0.25) +
  ggplot2::scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  ggplot2::guides(fill = guide_legend(reverse = T)) +
  ggplot2::theme_classic() +
  ggplot2::labs(x = "Posterior mean - true value",
                y = "Model") +
  ggplot2::theme( strip.background = element_blank(), 
                  axis.line = element_line(color = "black", size = 0.2),
                  axis.ticks = element_line(color = "black", size = 0.2),
                  legend.title = element_blank(), 
                  legend.position = "bottom",
                  axis.title.y = element_blank(),
                  axis.text.y = element_text( size = c(9,9,9,9,9,12), color = "black",
                                              face = c("plain", "plain", "plain", "plain", "plain", "bold")),
                  axis.title.x = element_text(size = 11, color = "black"),
                  axis.text.x = element_text(size = 9, color = "black"),
                  legend.text = element_text(size = 10, color = "black"),
                  strip.text = element_text(size = 8, color = "black"),
                  panel.background = element_rect(fill = "white", color = NA), 
                  plot.background = element_rect(fill = "white", color = NA))

setwd(here::here("figures"))
ggplot2::ggsave(
  "figure_03.png", 
  width = 6, 
  height = 6, 
  units = "in", 
  dpi = 300
)