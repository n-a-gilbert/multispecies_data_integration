# 16 March 2023
# Figure 5
# Visualize differences in number-of-groups and group-size between regions

library(here)
library(tidyverse)
library(MetBrewer)
library(MCMCvis)
library(patchwork)
library(reshape)

setwd(here::here("results"))
load("herbivore_case_study_results_v04.RData")

sp_key <- tibble::tibble(
  sp = 1:11, 
  sp_name =  c( "buffalo", "eland", "elephant",
                "giraffe", "Grant's", "hartebeest",
                "impala",  "Thomson's", "topi",
                "warthog", "waterbuck" ) )

sp_n <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("sp", "region", "iter")) %>% 
  bind_rows(.id = "param") %>% 
  mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  mutate( value = exp(value)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate( N = `Group size` * `Number of groups`) %>%
  dplyr::select(sp, region, iter, value = N) %>% 
  mutate(value = log(value)) %>% 
  full_join(sp_key) %>% 
  group_by(sp_name, region) %>% 
  dplyr::summarise( mean = mean(value), 
                    l95 = quantile(value, c(0.025)), 
                    u95 = quantile(value, c(0.975))) %>% 
  mutate(sp_name = factor(sp_name, 
                          levels = c("Grant's", "hartebeest", "eland", "giraffe", "elephant", "waterbuck", 
                                     "warthog", "buffalo", "topi", "Thomson's", "impala")))

sp_sig <-
  lapply(
    X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
    FUN = reshape::melt, 
    varnames = c("sp", "region", "iter")) %>% 
  bind_rows(.id = "param") %>% 
  mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  mutate( value = exp(value)) %>%
  as_tibble() %>% 
  pivot_wider(names_from = param, values_from = value) %>% 
  mutate( N = `Group size` * `Number of groups`) %>%
  dplyr::select(sp, region, iter, N) %>%
  pivot_wider(names_from = region, values_from = N) %>%
  mutate(diff = Talek - Mara) %>%
  group_by(sp) %>%
  summarise( mean = mean(diff),
             l95 = quantile(diff, c(0.025)),
             u95 = quantile(diff, c(0.975))) %>% 
  full_join(sp_key) %>% 
  mutate(region = ifelse( mean > 0, "Talek", "Mara")) %>% 
  mutate(sig = ifelse( l95 < 0 & u95 > 0, NA, "*")) %>% 
  dplyr::select(sp_name, sig, region) %>% 
  add_column( mean = -3.45)

com_n <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("mu_beta0", "mu_alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("region", "iter")) %>% 
  bind_rows(.id = "param") %>% 
  as_tibble() %>% 
  separate( region, into = c("junk", "region"), sep = "\\[") %>% 
  mutate(region = parse_number(region)) %>% 
  mutate( param = ifelse(param == "mu_alpha0", "Number of groups", "Group size")) %>%
  mutate( value = exp(value)) %>% 
  dplyr::select( param, region, iter, value ) %>%
  pivot_wider(names_from = param, values_from = value) %>%
  mutate( value = `Group size` * `Number of groups`) %>%
  mutate(value = log(value)) %>% 
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  
  group_by( region) %>% 
  dplyr::summarise( mean = mean(value), 
                    l95 = quantile(value, c(0.025)), 
                    u95 = quantile(value, c(0.975))) 

( nplot <- ggplot() +
    geom_rect(
      data = com_n,
      aes( xmin = l95,
           xmax = u95,
           fill = region),
      ymin = -Inf,
      ymax = Inf,
      alpha = 0.1) +
    geom_vline(
      data = com_n,
      aes(xintercept = mean,
          color = region),
      size = 1,
      alpha = 0.2) +
    geom_errorbar(
      data = sp_n, 
      aes(y = sp_name, 
          xmin = l95, 
          xmax = u95, 
          color = region),
      width = 0,
      size = 1,
      position = position_dodge(width = 0.5)) +
    geom_point(
      data = sp_n, 
      aes(y = sp_name,
          x = mean, 
          color = region),
      size = 3,
      position = position_dodge(width = 0.5)) +
    geom_text( data = sp_sig,
               aes(x = mean,
                   y = sp_name,
                   label = sig,
                   color = region),
               size = 6,
               vjust = 0.75,
               fontface = "bold",
               show.legend = FALSE) +
    scale_fill_manual("Region",
                      values = MetBrewer::MetPalettes$Greek[[1]][c(2, 5)]) +
    scale_color_manual( "Region",
                        values = MetBrewer::MetPalettes$Greek[[1]][c(2, 5)]) +
    scale_x_continuous( limits = c(-3.5, 5.5),
                        expand = c(0.02, 0)) +
    theme_minimal() +
    labs( x = " Log ( Abundance ) ",
          title = "(a)") + 
    theme( panel.grid = element_blank(),
           axis.title.y = element_blank(),
           legend.text = element_text(color = "black", size = 10), 
           legend.title = element_text(color = "black", size = 11),
           legend.position = c(0.8, 0.1),
           plot.title = element_text(color = "black", size = 10),
           axis.text = element_text(color = "black", size = 10),
           axis.title.x = element_text(color = "black", size = 11),
           axis.line = element_line(color = "black", size = 0.1),
           plot.background = element_rect(color = NA, fill = "white"), 
           panel.background = element_rect(color = NA, fill = "white")) +
    guides(color = guide_legend(reverse = TRUE), 
           fill = guide_legend(reverse = TRUE)) )

species_ng_gs <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("beta0", "alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("sp", "region", "iter")) %>% 
  bind_rows(.id = "param") %>% 
  mutate( param = ifelse(param == "alpha0", "Number of groups", "Group size")) %>%
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>%
  mutate( value = exp(value)) %>% 
  pivot_wider(names_from = region, values_from = value) %>% 
  mutate( diff = Talek - Mara) %>% 
  group_by(sp, param) %>% 
  summarise( mean = mean(diff), 
             l95 = quantile(diff, c(0.025)), 
             u95 = quantile(diff, c(0.975))) %>% 
  full_join(sp_key) %>% 
  mutate( param = factor(param, 
                         levels = c("Number of groups", "Group size")), 
          sp_name = factor(sp_name,
                           levels = c("Grant's", "hartebeest", "eland", "giraffe", "elephant", "waterbuck",
                                      "warthog", "buffalo", "topi", "Thomson's", "impala")))

com_ng_gs <- lapply(
  X = MCMCvis::MCMCpstr( out, params = c("mu_beta0", "mu_alpha0"), type = "chains"),
  FUN = reshape::melt, 
  varnames = c("region", "iter")) %>% 
  bind_rows() %>%
  as_tibble() %>% 
  mutate( param = ifelse( grepl("alpha0", region), "Number of groups", "Group size")) %>% 
  separate(region, into = c("junk", "region"), sep = "\\[") %>% 
  mutate(region = parse_number(region)) %>% 
  dplyr::mutate(region = ifelse(region == 1, "Mara", "Talek")) %>% 
  mutate( value = exp(value)) %>% 
  pivot_wider(names_from = region, values_from = value) %>% 
  mutate( diff = Talek - Mara)%>% 
  group_by(param) %>% 
  summarise( mean = mean(diff), 
             l95 = quantile(diff, c(0.025)), 
             u95 = quantile(diff, c(0.975))) %>% 
  mutate( param = factor(param, 
                         levels = c("Number of groups", "Group size")))


labs <- tibble(
  param = c("Number of groups", "Group size"),
  sp_name = c("hartebeest", "hartebeest"),
  mean = c(10, 10),
  label = c("More groups in Talek", "Larger groups in Talek")) %>% 
  mutate(sp_name = factor(sp_name,
                          levels = c("Grant's", "hartebeest", "eland", "giraffe", "elephant", "waterbuck",
                                     "warthog", "buffalo", "topi", "Thomson's", "impala"))) %>% 
  mutate( param = factor(param, 
                         levels = c("Number of groups", "Group size")))

( diffplot <- ggplot(  ) +
    geom_rect( data = com_ng_gs, 
               aes( xmin = l95, 
                    xmax = u95), 
               ymin = -Inf, 
               ymax = Inf,
               alpha = 0.2,
               fill = MetBrewer::MetPalettes$Hiroshige[[1]][7]) +
    geom_vline(data = com_ng_gs, 
               aes( xintercept = mean),
                color = MetBrewer::MetPalettes$Hiroshige[[1]][7],
               alpha = 0.2, 
               size = 1) +
    geom_vline( xintercept = 0,
                color = MetBrewer::MetPalettes$Hiroshige[[1]][1],
                linetype = "dashed" ) +
    geom_errorbar( data = species_ng_gs, 
                   aes( y = sp_name, xmin = l95, xmax = u95),
                   width = 0, size = 1,
                   color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    geom_point( data = species_ng_gs, 
                aes(y = sp_name, x = mean), 
                size  = 3, 
                color = MetBrewer::MetPalettes$Hiroshige[[1]][9] ) +
    geom_text(
      data = labs, 
      aes( y = sp_name, 
           x = mean, 
           label = label),
      size = 2.5,
      color = MetBrewer::MetPalettes$Hiroshige[[1]][1]) +
  facet_wrap(~param, ncol = 1)  +
    theme_minimal() +
    scale_x_continuous(breaks = c(-5, 0, 5, 10, 15)) +
    labs( x = " Talek - Mara ",
          title = "(b)") + 
    theme( panel.grid = element_blank(),
           axis.title.y = element_blank(),
           plot.title = element_text(color = "black", size = 10),
           axis.text = element_text(color = "black", size = 10),
           axis.title.x = element_text(color = "black", size = 11),
           axis.line = element_line(color = "black", size = 0.1),
           plot.background = element_rect(color = NA, fill = "white"), 
           panel.background = element_rect(color = NA, fill = "white")) )


nplot + diffplot + plot_layout(widths = c(1, 1))

setwd(here::here("figures"))
ggsave(
  "region_comparsion_v01.png", 
  width = 6.5,
  height = 5, 
  units = "in", 
  dpi = 300)