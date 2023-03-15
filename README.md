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
* [main_simulation_icm_relative_bias_v01.png](./figures/main_simulation_icm_relative_bias_v01.png) - Figure 3 - relative bias of abundance estimates from main simulation
* [study_area_map_v02.png](./figures/study_area_map_v02.png) - map of study area
* [code_for_figures](./figures/code_for_figures) - folder with scripts to create figures
   * [figure_02_study_area_map_v01.R](./figures/code_for_figures/figure_02_study_area_map_v01.R) - make study area map
   * [figure_s1_s2_simulated_communit_v01.R](./figures/code_for_figures/figure_s1_s2_simulated_communit_v01.R) - plot simulated community example 
