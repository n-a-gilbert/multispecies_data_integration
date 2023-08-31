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
* [count_data_v01.RData](./data/count_data_v01.RData) - formatted count data. This .RData file contains 1 object:
  * **transect_data**: a dataframe containing the following columns:
 
    | Variable name | Meaning |
    |---------------|---------|
    | transect | Transect name |
    | sp_name | Common name of species |
    | date | Date of survey |
    | sp | Species id |
    | site | Site (transect} id |
    | rep | Visit id |
    | count | Count of the total number of individuals of a species observed on a survey |
    | area | Area offset for transect |
    | region | Binary variable indicating Mara (0) or Talek (1) region |

* [distance_sampling_data_v01.RData](./data/distance_smapling_data_v01.RData) - formatted distance sampling data. This .RData file contains 3 objects:
  * **b**: A scalar, the maximum distance to which animals are counted (1000 m)
  * **mdpt**: A vector, the distance (in m) to the midpoint of each distance bin from the transect line
  * **v**: A scalar, the width (in m) of the distance bins
  * **final2**: a dataframe containing the following columns:
    
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
  | sp_name | Common name of species |

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
* [herbivore_case_study_results_v04.RData](./resultsherbivore_case_study_results_v04.RData) - Model output for Mara herbivores case study. This .RData contains 4 objects
  * **constants** A list of constants used in Nimble model

    | Variable name | Meaning |
    |---------------|---------|
    | NSPECIES | Number of species |
    | NBINS | Number of distance bins |
    | NDISTANCES | Number of distance observations |
    | NSURVEYS | Number of distance sampling surveys |
    | NCOUNTS | Number of count surveys |
    | SP_GS | Species index for the group size data |
    | SP_NG | Species index for the number-of-group data |
    | SP_TC | Species index for the count data |
    | REGION_GS | Region index for the group size data |
    | REGION_NG | Region index for the number-of-group data |
    | NREGION | Number of regions |

  * **data** A list of data used in the Nimble model
 
    | Variable name | Meaning |
    |---------------|---------|
    | MIDPOINT | Distance to the midpoint of each distance bin |
    | yGS | Observed group sizes |
    | DCLASS | Observed distance class |
    | B_DS | Maximum distance to which animals are counted for distance sampling |
    | B_TC | Maximum distance to which animals are counted for counts |
    | V | Width of distance bins
    | yNG | Observed number of groups (distance sampling) |
    | yN_DS | Observed total number of animals (distance sampling |
    | yN | Observed total number of animals (counts) |
    | OFFSET_DS | Area offset for distance sampling transects |
    | OFFSET_TC | Area offset for count transects |
    | MASS | Body mass of each species |

  * **out** A list of the MCMC chains with the posterior samples for model parameters
  * **model.code** Code for the Nimble model

* [main_simulation_results_v01.RData](./results/main_simulation_results_v01.RData) - Summarized results from main simulation. This .RData file contains 4 dataframes:
  * **icm_a1**. This summarizes model estimates of a covariate effect.

    | Variable name | Meaning |
    |---------------|---------|
    | simrep | Which replicate simulation (1 through 100) |
    | sp | Species id |
    | mean | Posterior mean of difference between estimate and true value |
    | sd | Posterior standard deviation of difference between estimate and true value |
    | l95 | Lower bound of 95% credible interval of difference between estimate and true value |
    | u95 | Upper bound of 95% credible interval of difference between estimate and true value |
    | contain0 | Does the 95% credible interval of the difference between the estimate and true value contain 0? |
    | nsites_tc_fact | Which relative data amount scenario? 1 = same number of distance sampling and count transects, 4 = 4x more count transects |
    | p_bias | Which relative detectability scenario?  -1 = count detectability 10% lower, 0 = same detectability |
    
  * **icm_info**. True values of parameters used to simulate data

    | Variable name | Meaning |
    |---------------|---------|
    | sp | Species id |
    | totDS | Total number of individuals at distance sampling sites |
    | totTC | Total number of individuals at count sites |
    | alpha0 | Number-of-groups intercept |
    | alpha1 | Effect of covariate on number-of-groups |
    | beta0 | Group size intercept |
    | gamma0 | Intercept for scale parameter in detection function |
    | zeta | Hyperparameter for number-of-group overdispersion |
    | xi | Hyperparameter for group size overdispersion |
    | simrep | Which replicate simulation (1 through 100) |
    
  * **icm_n_ds_rb*.* Relative bias of estimated abundance underlying distance sampling data.

    | Variable name | Meaning |
    |---------------|---------|
    | simrep | Which replicate simulation (1 through 100) |
    | sp | Species id |
    | mean | Posterior mean of relative bias |
    | sd | Posterior sd of relative bias |
    | l95 | Lower bound of 95% credible interval for relative bias |
    | u95 | Upper bound of 95% credible interval for relative bias |
    | contain0 | Does the 95% credible interval for relative bias contain 0 (no bias) |
    | nsites_tc_fact | Which relative data amount scenario? 1 = same number of distance sampling and count transects, 4 = 4x more count transects |
    | p_bias | Which relative detectability scenario?  -1 = count detectability 10% lower, 0 = same detectability |

  * **icm_n_tc_rb**. Relative bias of estimated abundance underlying count data. Same variable definitions as above.   
* [simulation_alternative_model_results_v01.RData](./results/simulation_alternative_model_results_v01.RData) - Summarised results for alternative model simulations

