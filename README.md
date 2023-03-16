# A multispecies hierarchical model to integrate count and distance sampling data

### [Neil A. Gilbert](https://gilbertecology.com), [Caroline M. Blommel](https://www.researchgate.net/profile/Caroline-Blommel), [Matthew T. Farr](https://farrmt.github.io/), [David S. Green](https://scholar.google.com/citations?user=zZf1ct0AAAAJ), [Kay E. Holekamp](https://www.holekamplab.org/), [Elise F. Zipkin](https://zipkinlab.org/)

### Data/code DOI: To be added upon acceptance

#### Please contact the first author for questions about the code or data: Neil A. Gilbert (neil.allen.gilbert@gmail.com)
__________________________________________________________________________________________________________________________________________

## Abstract

To be added upon submission

## Repository Directory

### [code](./code): Contains code for preparing case study data, running case study model, and simulations
*  [case_study_analysis](./code/case_study_analysis)
   * [herbivore_case_study_analaysis_v01.R](./code/case_study_analysis/herbivore_case_study_analysis_v01.R) - code to run case study model
*  [data_processing](./code/data_processing)
   * [prepare_distance_sampling_data_v01.R](./code/data_processing/prepare_distance_sampling_data_v01.R) - format case study distance sampling data
   * [prepare_count_data_v01.R](./code/data_processing/prepare_count_data_v01.R) - format case study count data
* [simulations](./code/simulations)
   * [alternative_model_comparison](./code/simulations/alternative_model_comparison) - folder containing scripts to run simulations for alternative single datastream / single-species models
      * [community_count_v01.R](./code/simulations/alternative_model_comparison/community_count_v01.R) - community count-only model
      * [community_distance_sampling_v01.R](./code/simulations/alternative_model_comparison/community_distance_sampling_v01.R) - community distance sampling-only model
      * [single_species_common_count_v01.R](./code/simulations/alternative_model_comparison/single_species_common_count_v01.R) - single species (common) count only model
      * [single_species_common_distance_sampling_v01.R](./code/simulations/alternative_model_comparison/single_species_common_distance_sampling_v01.R) - single species (common) distance sampling only model
      * [single_species_common_integrated_v01.R](./code/simulations/alternative_model_comparison/single_species_common_integrated_v01.R) - single species (common) integrated model
      * [single_species_rare_count_v01.R](./code/simulations/alternative_model_comparison/single_species_rare_count_v01.R) - single species (rare) count only model
      * [single_species_rare_distance_sampling_v01.R](./code/simulations/alternative_model_comparison/single_species_common_distance_sampling_v01.R) - single species (rare) distance sampling only model
      * [single_species_rare_integrated_v01.R](./code/simulations/alternative_model_comparison/single_species_rare_integrated_v01.R) - single species (rare) integrated model
   * [main_simulation_v01.R](./code/simulations/main_simulation_v01.R) - script to run the main simulation

### [data](./data): Contains data for case study
* [Shapefiles](./data/Shapefiles) - various shapefiles
  * [DS](./data/Shapefiles/DS) - shapefiles for distance sampling transects
  * [Transects](./data/Shapefiles/Transects) - shapefiles for transects where count data was collected
  * [reserve](./data/Shapefiles/reserve) - shapefile for reserve / management zone boundaries
* [Herbivore Utilization Complete.csv](./data/Herbivore%20Utilization%20Complete.csv) - unformatted distance sampling data
* [count_data_v01.RData](./data/count_data_v01.RData) - formatted count data
* [distance_sampling_data_v01.RData](./data/distance_smapling_data_v01.RData) - formatted distance sampling data
* [tblPreyCensus_2012to2014.csv](./data/tblPreyCensus_2012to2014.csv) - unformatted count data 

### [figures](./figures): contains figures, and code to create them
* [density_estimates_mean_v01.png](./figures/density_estimates_mean_v01.png) - Fig. S3 - posterior means of abundance estimates
* [density_estimates_sd_v01.png](./figures/density_estimates_sd_v01.png) - Fig. S4 - posterior standard deviation of abundance estimates
* [figure_01.png](./figures/figure_01.png) - Figure 1 - conceptual overview of model
* [figure_05_contrasts_v01.png](./figures/figure_05_contrasts_v01.png) - Figure 5 - group size & number of group differences between regions
* [main_simulation_icm_relative_bias_v01.png](./figures/main_simulation_icm_relative_bias_v01.png) - Figure 3 - relative bias of abundance estimates from main simulation
* [simulation_model_comparison_covariate_v01.png](./figures/simulation_model_comparison_covariate_v01.png) - Figure 4 - accuracy/precision of covariate estimates from ICM & alternative models
* [study_area_map_v02.png](./figures/study_area_map_v02.png) - map of study area
* [code_for_figures](./figures/code_for_figures) - folder with scripts to create figures
   * [figure_02_study_area_map_v01.R](./figures/code_for_figures/figure_02_study_area_map_v01.R) - make study area map
   * [figure_5_region_comparison_v01.R](./figures/code_for_figures/figure_5_region_comparison_v01.R) - code for fig. 5 (differences between regions)
   * [figure_s1_s2_simulated_communit_v01.R](./figures/code_for_figures/figure_s1_s2_simulated_communit_v01.R) - plot simulated community example 
   * [figures_3_4_tables_s2_s3_v01.R](./figures/code_for_figures/figures_3_4_tables_s2_s3_v01.R) - create figures 3 & 4, plus tables s2 and s3
   * [figures_s3_s4_plot_abundance_estimates_v01.R](./figures/code_for_figures/figures_s3_s4_plot_abundance_estimates_v01.R) - creates Figs. S3, S4
   
### [results](./results): contains results files
* [herbivore_case_study_results_v01.RData](./results/herbivore_case_study_results_v01.RData) - Model output for Mara herbivores case study
* [main_simulation_results_v01.RData](./results/main_simulation_results_v01.RData) - Summarized results from main simulation
* [simulation_alternative_model_results_v01.RData](./results/simulation_alternative_model_results_v01.RData) - Summarised results for alternative model simulations
