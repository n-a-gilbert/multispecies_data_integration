library(here)
library(tidyverse)
library(MetBrewer)

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
  dplyr::mutate( model_name = ifelse( model == "cs", "Single species, counts",
                                      ifelse(model == "ds", "Single species, distance sampling",
                                             ifelse(model == "is", "Single species, integrated", 
                                                    ifelse( model == "cc", "Community, counts", 
                                                            ifelse(model == "dc", "Community, distance sampling", 
                                                                   ifelse(model == "ic", "Community, integrated", NA))))))) |> 
  
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
  ggplot2::labs(x = "Estimate - truth",
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