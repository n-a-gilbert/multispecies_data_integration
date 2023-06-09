library(here)
library(tidyverse)
library(MetBrewer)

setwd(here::here("results"))
load("main_simulation_results_v01.RData")

prop_labs <- icm_n_ds_rb %>% 
  add_column( source = "Distance sampling") %>% 
  full_join(icm_info) %>% 
  full_join( add_column( full_join(icm_n_tc_rb, icm_info), source = "Counts" ) ) %>% 
  mutate(tot = ifelse(source == "Counts", totTC, totDS)) %>%
  dplyr::select( simrep, sp, mean, sd, l95, u95, contain0, tot, nsites_tc_fact, p_bias, source) %>% 
  mutate(p_bias = ifelse(p_bias == 0, "Same    \ndetectability", "Count    \ndetectability\n20% lower")) %>% 
  mutate(p_bias = factor(p_bias, levels = c(
    "Count    \ndetectability\n20% lower",
    "Same    \ndetectability"))) %>% 
  mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) %>% 
  mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
    "4x more count data",
    "Same amount of both data"))) %>% 
  mutate( contain0 = ifelse(contain0 == 0, "False", "True")) %>% 
  ungroup() %>% 
  group_by( source, nsites_tc_fact, p_bias ) %>% 
  mutate(n = n()) %>% 
  summarise( prop = round (sum( contain0 == "False", na.rm = TRUE) / n, 2 )) %>% 
  mutate(prop = paste0( (100 - (prop * 100)), "% unbiased")) %>% 
  distinct(.) %>% 
  tibble::add_column( mean = c( 1.65, 1.65, 
                                1.65, 1.65, 
                                1.65, 1.65, 
                                1.65, 1.65))

( panela <- icm_n_ds_rb %>% 
    add_column( source = "Distance sampling") %>% 
    full_join(icm_info) %>% 
    full_join( add_column( full_join(icm_n_tc_rb, icm_info), source = "Counts" ) ) %>% 
    mutate(tot = ifelse(source == "Counts", totTC, totDS)) %>%
    dplyr::select( simrep, sp, mean, sd, l95, u95, contain0, tot, nsites_tc_fact, p_bias, source) %>% 
    mutate(p_bias = ifelse(p_bias == 0, "Same    \ndetectability", "Count    \ndetectability\n20% lower")) %>% 
    mutate(p_bias = factor(p_bias, levels = c(
      "Count    \ndetectability\n20% lower",
      "Same    \ndetectability"))) %>% 
    mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) %>% 
    mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
      "4x more count data",
      "Same amount of both data"))) %>% 
    mutate( contain0 = ifelse(contain0 == 0, "False", "True")) %>% 
    ggplot( aes( x = mean, y = p_bias, fill = nsites_tc_fact )) +
    facet_wrap(~source, scales = "free_x") +
    geom_vline(xintercept = 0,
               color = "gray60",
               linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.2,
                 outlier.size = 0.75, 
                 size = 0.25) +
    geom_label(data = prop_labs,
               aes(label = prop,
                   color = nsites_tc_fact),
               position = position_dodge(width = 1.2),
               size = 3,
               fill = "white",
               show.legend = FALSE,
               label.size = NA) +
    scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    scale_color_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    theme_classic() +
    xlim(c(-1, 3.55)) +
    labs(x = "Relative bias (%)",
         title = "(a)                    Abundance") +
    theme(legend.position = "bottom",
          legend.title = element_blank(), 
          axis.title.y =element_blank(),
          axis.text = element_text(size = 10, color = "black"), 
          axis.title = element_text(size = 10, color = "black"), 
          strip.text = element_text(size = 11, color = "black"),
          legend.text = element_text(size = 10, color = "black"),
          plot.title = element_text(size = 11, color = "black"),
          panel.background = element_rect(fill = "white", color = NA), 
          plot.background = element_rect(fill = "white", color = NA),
          legend.margin = margin(0, 0, 0, 0), 
          legend.box.margin = margin(-5, 0, 0, 0),
          strip.background = element_rect(color = NA),
          axis.line = element_line(size = 0.1, color = "black"),
          axis.ticks = element_line(size = 0.1, color = "black")) +
    guides( fill = guide_legend(nrow = 1,
                                reverse = TRUE)))

( panelb <- icm_a1 %>% 
    mutate(p_bias = ifelse(p_bias == 0, 
                           "Same     \ndetectability", 
                           "Count   \ndetectability\n20% lower")) %>% 
    mutate(p_bias = factor(p_bias, 
                           levels = c(
                             "Count   \ndetectability\n20% lower",
                             "Same     \ndetectability"))) %>% 
    mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, "Same amount of both data", "4x more count data")) %>% 
    mutate(nsites_tc_fact = factor(nsites_tc_fact, levels = c(
      "4x more count data",
      "Same amount of both data"))) %>% 
    tibble::add_column( type = "Covariate effect") %>% 
    
    ggplot( aes( x = mean, y = p_bias, fill = nsites_tc_fact)) +
    facet_wrap(~type) +
    geom_vline(xintercept = 0,
               color = "gray60",
               linetype = "dashed") +
    geom_boxplot(outlier.alpha = 0.2,
                 outlier.size = 0.75, 
                 size = 0.25) +
    labs( x = "Truth - estimate") +
    scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(1,3)])+
    ggtitle("(b)    Covariate effect") +
    theme_classic() +
    theme(legend.position = "bottom",
          legend.title = element_blank(), 
          axis.title.y =element_blank(),
          axis.text.x = element_text(size = 10, color = "black"),
          axis.text.y = element_blank(), 
          axis.title.x = element_text(size = 10, color = "black"), 
          strip.text = element_text(size = 11, color = "white"),
          legend.text = element_text(size = 10, color = "black"),
          panel.background = element_rect(fill = "white", color = NA), 
          plot.background = element_rect(fill = "white", color = NA),
          legend.margin = margin(0, 18, 0, 0), 
          plot.title = element_text(size = 11, color = "black"),
          legend.box.margin = margin(-5, 18, 0, 0),
          strip.background = element_rect(color = NA),
          axis.line = element_line(size = 0.1, color = "black"),
          axis.ticks.x = element_line(size = 0.1, color = "black"),
          axis.ticks.y = element_blank()) +
    guides(fill = guide_legend(nrow = 1, 
                               reverse = TRUE)))


library(patchwork)

panela + panelb + plot_layout( nrow = 1, widths = c(2, 1), guides = "collect") &
  theme(legend.position='bottom',
        legend.margin = margin(0, 0, 0, 0), 
        legend.box.margin = margin(t = -5, r = 0, b = 0, l = 0)) 

setwd(here::here("figures"))
ggsave(
  "figure_02.png",
  width = 6, 
  height = 4, 
  units = "in", 
  dpi = 300
)

# Table S2 - accuracy and precision of ICM's covariate estimates
icm_a1 %>% 
  mutate(p_bias = ifelse(p_bias == 0, 
                         "Same detectability", 
                         "Count detectability 20% lower")) %>% 
  mutate(p_bias = factor(p_bias, 
                         levels = c(
                           "Same detectability",
                           "Count detectability 20% lower"
                         ))) %>% 
  mutate(nsites_tc_fact = ifelse(nsites_tc_fact == 1, 
                                 "Same amount of both data", 
                                 "4x more count data")) %>% 
  mutate(nsites_tc_fact = factor(nsites_tc_fact, 
                                 levels = c(
                                   "Same amount of both data",
                                   "4x more count data"
                                 ))) %>% 
  group_by( nsites_tc_fact, p_bias ) %>%
  mutate(n = n()) %>% 
  summarise( unbiased = round( 100 * ( sum(contain0) / n), 1),
             biased = round( 100 - unbiased, 1),
             mean_sd = round( mean(sd), 2)) %>% 
  distinct() %>% 
  arrange(nsites_tc_fact, p_bias) %>% 
  mutate(unbiased = paste0(unbiased, "%"), 
         biased = paste0(biased, "%")) %>% 
  rename( `Data amount` = nsites_tc_fact, 
          `Detectability` = p_bias, 
          Unbiased = unbiased, 
          Biased = biased, 
          SD = mean_sd )


setwd(here::here("results"))
load("simulation_alternative_model_results_v01.RData")

icm_for_comparison <- icm_info %>%
  group_by( simrep ) %>% 
  arrange(simrep, totDS) %>% 
  filter(! totDS == max(totDS) ) %>% 
  slice(1, n()) %>% 
  mutate(species = ifelse( totDS == min(totDS), "rare", "common")) %>% 
  full_join( distinct( dplyr::select(icm_n_ds_rb, simrep, nsites_tc_fact, p_bias))) %>% 
  filter(nsites_tc_fact == 1 & p_bias == 0) %>% 
  add_column(model = "ICM") %>% 
  dplyr::select( model, species, sp, simrep ) %>% 
  left_join(icm_a1) %>% 
  ungroup() %>% 
  mutate(contain0 = ifelse(contain0 == 1, "contains 0", "does not contain 0")) %>% 
  dplyr::select(model, species, mean, sd, l95, u95, contain0 )

key <- tribble(
  ~model, ~name,
  "ICM", "Community, integrated",
  "CC", "Community, counts",
  "CDS", "Community, distance sampling",
  "ISS", "Single species, integrated", 
  "SSC", "Single species, counts", 
  "SSDS", "Single species, distance sampling")

alpha1_truth_minus_estimate %>% 
  full_join(icm_for_comparison) %>% 
  full_join(key) %>% 
  mutate(species = ifelse(species == "common", "Common species", "Rare species")) %>% 
  mutate(species = factor(species, levels = c("Rare species", "Common species"))) %>% 
  mutate( name = factor(name, 
                        levels = c(
                          "Community, integrated",
                          "Community, distance sampling", 
                          "Community, counts",
                          "Single species, integrated", 
                          "Single species, distance sampling",
                          "Single species, counts"))) %>% 
  ggplot(
    aes( x = mean, y = rev(name), fill = species)) +
  scale_y_discrete(expand = c(0, 0)) +
  geom_vline(xintercept = 0,
             color = MetPalettes$Hiroshige[[1]][c(1)]) +
  geom_boxplot( outlier.alpha = 0.2,
                outlier.size = 0.75, 
                size = 0.25) +
  scale_color_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  scale_fill_manual(values = MetPalettes$Hiroshige[[1]][c(7,9)])+
  theme_classic() +
  labs(x = "Covariate effect: truth - estimate") +
  theme(
    legend.position = "bottom",
    legend.title = element_blank(), 
    axis.title.y = element_blank(),
    axis.text.y = element_text( size = c(10, 10, 10, 10, 10, 14), color = "black",
                                face = c("plain", "plain", "plain", "plain", "plain", "bold")),
    axis.title = element_text(size = 11, color = "black"), 
    legend.text = element_text(size = 10, color = "black"),
    panel.background = element_rect(fill = "white", color = NA), 
    plot.background = element_rect(fill = "white", color = NA),
    legend.margin = margin(0, 12.5, 0, 0), 
    legend.box.margin = margin(-10, 12.5, 0, 0),
    axis.line = element_line(size = 0.1, color = "black"), 
    axis.ticks = element_line(size = 0.1, color = "black")) + 
  guides(color = guide_legend(reverse = TRUE)) +
  scale_y_discrete(limits = rev)

setwd(here::here("figures"))
ggsave(
  "figure_03.png",
  width = 4.75, 
  height = 4, 
  units = "in", 
  dpi = 300
)

# relative precision of models
# Table S3
alpha1_truth_minus_estimate %>% 
  full_join(icm_for_comparison) %>% 
  full_join(key) %>% 
  mutate(species = ifelse(species == "common", "Common species", "Rare species")) %>% 
  mutate(species = factor(species, levels = c("Rare species", "Common species"))) %>% 
  mutate( name = factor(name, 
                        levels = c(
                          "Community, integrated",
                          "Community, distance sampling", 
                          "Community, counts",
                          "Single species, integrated", 
                          "Single species, distance sampling",
                          "Single species, counts"))) %>% 
  group_by(name, species) %>% 
  summarise(mean_sd = mean(sd)) %>% 
  ungroup() %>% 
  group_by(species) %>% 
  mutate( rel_sd = mean_sd / min(mean_sd) ) %>% 
  dplyr::select(name, species, rel_sd) %>% 
  pivot_wider(names_from = species, values_from = rel_sd) %>% 
  rename(Model = name)
