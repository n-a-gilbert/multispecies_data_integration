# A multispecies hierarchical model to integrate count and distance sampling data

### [Neil A. Gilbert](https://gilbertecology.com), [Caroline M. Blommel](https://www.researchgate.net/profile/Caroline-Blommel), [Matthew T. Farr](https://farrmt.github.io/), [David S. Green](https://scholar.google.com/citations?user=zZf1ct0AAAAJ), [Kay E. Holekamp](https://www.holekamplab.org/), [Elise F. Zipkin](https://zipkinlab.org/)

### Data/code DOI: To be added upon acceptance

#### Please contact the first author for questions about the code or data: Neil A. Gilbert (neil.allen.gilbert@gmail.com)
__________________________________________________________________________________________________________________________________________

## Abstract

Integrated community models—an emerging framework in which multiple data sources for multiple species are analyzed simultaneously—offer opportunities to expand inferences beyond the single-species and single-data source approaches common in ecology. We developed a novel integrated community model combining distance sampling and single-visit count data; within the model, information is shared among data sources (via a joint likelihood) and species (via a random effects structure) to estimate abundance patterns across a community. Simulations showed that the integrated community model produced more precise estimates of ecological quantities such as covariate effects than alternative single-species and single-data source models. The model provided unbiased estimates of abundance—even for locations with single-visit count data that contain no information about the observation process—assuming comparable detection probabilities between data sources. When detection probabilities for simulated count data were different from distance sampling, however, abundance estimates were systematically biased. We applied the model to datasets on 11 herbivore species from the Masai Mara National Reserve, Kenya, and found considerable interspecific variation in response to local wildlife management practices: four species showed higher abundances in a region with passive conservation enforcement, four species showed higher abundances in a region with active conservation enforcement, and the remaining three species showed no abundance differences between the two regions. Furthermore, given the hierarchical structure of the model, we identified several species that showed between-region differences in group size and number of groups that were of greater magnitude than the community average. Future applications of this modeling framework should consider the circumstances under which data integration is appropriate given assumptions about shared abundance patterns and detection probabilities between data sources.

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
* [distance_sampling_data_v01.RData](./data/distance_smapling_data_v01.RData) - formatted distance sampling data. This .RData file contains 3 objects:
  * final2: a dataframe containing the following columns: 
  | Variable name | Meaning |
  |---------------|---------|
  | sp | Species id |
  | site | Site (transect) id |
  | rep | Visit id |
  | gs | Observed group size |
  | dclass | Distance class (1 through 40) of observed group |
  | ng | Observed number of groups for species x site x rep combo |
  | area | Area offset for transect |
  | region | Binary variable indicating Mara (0) or Talek (1) region |
  | date | Date of survey |
  | sp_name | common name of species |

* [tblPreyCensus_2012to2014.csv](./data/tblPreyCensus_2012to2014.csv) - unformatted count data 

### [figures](./figures): contains figures, and code to create them
* [code_for_figures](./figures/code_for_figures) - folder with scripts to create figures
   * [figure_02_03.R](./figures/code_for_figures/figure_02_03.R) - create figures 2 & 3, plus tables s2 and s3
   * [figure_04.R](./figures/code_for_figures/figure_04.R) - code for fig. 4 (differences between regions)
   * [figure_s1_s2.R](./figures/code_for_figures/figure_s1_s2.R) - plot simulated community example 
   * [figure_s3.R](./figures/code_for_figures/figure_s3.R) - make study area map
   * [figure_s4_s5.R](./figures/code_for_figures/figure_s4_s5.R) - creates Figs. S4, S5
* [figure_01.png](./figures/figure_01.png) - Figure 1 - conceptual overview of model
* [figure_02.png](./figures/figure_02.png) - Figure 2 - main simulation results
* [figure_03.png](./figures/figure_03.png) - Figure 3 - simulation - comparison to alternative models
* [figure_04.png](./figures/figure_04.png) - Figure 4 - case study results
* [figure_s1.png](./figures/figure_s1.png) - Figure S1 - simulated covariate effect
* [figure_s2.png](./figures/figure_s2.png) - Figure S2 - simulated detection function
* [figure_s3.png](./figures/figure_s3.png) - Figure S3 - case study map
* [figure_s4.png](./figures/figure_s4.png) - Figure S4 - case study - transect-level density estimates
* [figure_s5.png](./figures/figure_s5.png) - Figure S5 - case study - transect-level density estimate uncertainty

### [results](./results): contains results files
* [herbivore_case_study_results_v04.RData](./resultsherbivore_case_study_results_v04.RData) - Model output for Mara herbivores case study
* [main_simulation_results_v01.RData](./results/main_simulation_results_v01.RData) - Summarized results from main simulation
* [simulation_alternative_model_results_v01.RData](./results/simulation_alternative_model_results_v01.RData) - Summarised results for alternative model simulations

