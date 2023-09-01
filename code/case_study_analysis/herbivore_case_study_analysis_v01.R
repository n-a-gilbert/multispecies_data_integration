library(here)
library(nimble)
library(parallel)
library(tidyverse)

setwd(here::here("data"))

load("distance_sampling_data_v01.RData")
load("count_data_v01.RData")

# format data to get number-of-group data
ng_data <- final2 %>%
  dplyr::select( sp, site, rep, ng, gs, area, region ) |>
  dplyr::filter( !base::is.na( ng ) ) |>
  dplyr::group_by( sp, site, rep, area, region ) |>
  dplyr::summarise( ng = base::unique( ng ),
                    yN_DS = base::sum(gs, na.rm = TRUE ) ) |>
  dplyr::arrange( sp, site, rep )

# format data to get group size data
gs_data <- final2 |>
  dplyr::filter( !base::is.na( gs ) ) |>
  dplyr::arrange( sp, site, rep )

# from PanTHERIA - species mass
sp_mass <- tibble::tibble(
  sp_name = c("buffalo", "eland", "elephant", "giraffe", "grants",
              "hartebeest", "impala",  "thomsons", "topi",
              "warthog", "waterbuck"),
  mass = c(592666, 562592.7, 3824540, 964654.7, 55464.46,
           160937.9,  52591.69, 22907.43, 136000.3,
           82499.99, 204393.5))

det_attributes <- final2 %>%
  dplyr::select( sp, sp_name ) |>
  dplyr::distinct() |>
  dplyr::full_join( sp_mass ) |>
  dplyr::mutate( mass = base::as.numeric( base::scale( base::log( mass ) ) ) )

#package data up for nimble model
data <- list(
  MIDPOINT = mdpt,  # midpoint of distance class bins
  yGS = gs_data$gs, # group size counts (vector)
  DCLASS = gs_data$dclass, # distance class of group size counts (vector)
  B_DS = b,  # geometry for distance function
  B_TC = 100, # width of count transect
  V = v,   # geometry for distance function 
  yNG = ng_data$ng, # number of groups per site x rep x species
  yN_DS = ng_data$yN_DS, # total abundance on DS transect
  yN = transect_data$count, # transect counts
  OFFSET_DS = ng_data$area, # area offset for distance sampling transects
  OFFSET_TC = transect_data$area, # area offset for count transects
  MASS = det_attributes$mass ) # species body mass

constants <- list(
  NSPECIES = length(unique(final2$sp)), # number of species (11)
  NBINS = length(mdpt), # number of distance bins 
  NDISTANCES = nrow(gs_data), # number of group size observations
  NSURVEYS = nrow(ng_data),  # number of number-of-group observations
  NCOUNTS = nrow(transect_data), # number of transect counts
  SP_GS = gs_data$sp, # species ID, used for indexing group size data
  SP_NG = ng_data$sp,# species ID, used for indexing number-of-groups data
  SP_TC = transect_data$sp, # species ID, used for indexing transect counts
  REGION_GS = gs_data$region + 1, # Region = Mara (1) or Talek (2)
  REGION_NG = ng_data$region + 1, # Region - Mara (1) or Talek (2)
  REGION_TC = transect_data$region + 1, # Region  - Mara (1) or Talek (2)
  NREGION = 2) # number of regions

# NIMBLE code
model.code <- nimble::nimbleCode({
  mu_gamma0 ~ dnorm(5.5, sd = 0.5) # pretty informative prior for the scale parameter intercept
  # weakly informative priors for everything else
  sd_gamma0 ~ dexp(1)
  gamma1 ~ dnorm(0, sd = 2)
  for( j in 1:NREGION) {
    mu_alpha0[j] ~ dnorm(0, sd = 2)
    sd_alpha0[j] ~ dexp(1)
    mu_beta0[j]  ~ dnorm(0, sd = 2)
    sd_beta0[j]  ~ dexp(1)
  }
  # compute contrasts in intercepts for number-of-groups and group size between regions, community level
  mu_alpha0_contrast <- mu_alpha0[2] - mu_alpha0[1]
  sd_alpha0_contrast <- sd_alpha0[2] - sd_alpha0[1]
  mu_beta0_contrast <- mu_beta0[2] - mu_beta0[1]
  sd_beta0_contrast <- sd_beta0[2] - sd_beta0[1]
  # intercepts for species x region
  for (s in 1:NSPECIES) {
    for (j in 1:NREGION) {
      alpha0[s, j] ~ dnorm(mu_alpha0[j], sd = sd_alpha0[j])
      beta0[s, j] ~ dnorm(mu_beta0[j], sd = sd_beta0[j])
    }
    alpha0_contrast[s] <- alpha0[s, 2] - alpha0[s, 1]
    beta0_contrast[s] <- beta0[s, 2] - beta0[s, 1]
    # scale parameter intercept for detection function
    gamma0[s] ~ dnorm(mu_gamma0, sd = sd_gamma0)
    omega[s] <- exp(gamma0[s] + gamma1 * MASS[s]) # scale parameter for detection function
    pie_sp[s] <- sum(pie_ds[1:NBINS, s]) # transect-level deteciton probability for distance sampling
    pie_sp_tc[s] <- sum(pie_tc[1:4, s]) # transect-level detection probability for counts (only 4 first bins b/c only surveyed out to 100 m)
    zeta[s] ~ dgamma(4, 2) # number-of-group overdispersion hyperparameter
    xi[s] ~ dgamma(4, 2) # group size hyperparameter
    for (k in 1:NBINS) { # loop through distance classes
      log(g[k, s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s] * omega[s])
      pie_ds[k, s] <- g[k, s] * (V/B_DS)
      pie_tc[k, s] <- g[k, s] * (V/B_TC)
      pie_cell[k, s] <- pie_ds[k, s]/pie_sp[s] # cell probabilties, standardized to sum to 1
    }
  }
  # number-of-groups submodel
  for (i in 1:NSURVEYS) {
    log(lambda[i]) <- alpha0[SP_NG[i], REGION_NG[i]] + log(OFFSET_DS[i]) # expected number of groups
    rho[i] ~ dgamma(zeta[SP_NG[i]], zeta[SP_NG[i]]) # overdispersion
    NG[i] ~ dpois(lambda[i] * rho[i]) # latent number of groups
    yNG[i] ~ dbin(pie_sp[SP_NG[i]], NG[i]) # observed number of groups modeled w/ binomial distribution
    NG_missed[i] <- NG[i] - yNG[i] # how many groups were missed?
    log( kappa_ng[i] ) <- beta0[SP_NG[i], REGION_NG[i]] # average group size
    GS_missed[i] ~ T( dpois( kappa_ng[i] ), 1, ) # size of missed groups
    N_missed[i] <- NG_missed[i] * GS_missed[i] # how many total animals were missed
    N_ds[i] <- yN_DS[i] + N_missed[i] # total abundance as derived variable
  }
  # loop through the distance observations
  for (i in 1:NDISTANCES) {
    log(kappa[i]) <- beta0[SP_GS[i], REGION_GS[i]] # expected group size
    phi[i] ~ dgamma(xi[SP_GS[i]], xi[SP_GS[i]]) # group size overdispersion parameter
    yGS[i] ~ T(dpois(kappa[i] * phi[i]), 1, ) # model group size as draw from zero-truncated poisson
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i]]) # distance class modeled as draw from categorical distribution
  }
  # count submodel
  for (i in 1:NCOUNTS) {
    rho_tc[i] ~ dgamma(zeta[SP_TC[i]], zeta[SP_TC[i]]) # number-of-groups overdispersion
    phi_tc[i] ~ dgamma(xi[SP_TC[i]], xi[SP_TC[i]]) # group size overdispersion
    NG_tc[i] <- exp(alpha0[SP_TC[i], REGION_TC[i]] + log(OFFSET_TC[i])) * rho_tc[i] # number of groups
    GS_tc[i] <- exp(beta0[SP_TC[i], REGION_TC[i]]) * phi_tc[i] # group size
    N[i] ~ dpois(NG_tc[i] * GS_tc[i]) # total abundance is just number-of-groups times group size
    yN[i] ~ dbin(pie_sp_tc[SP_TC[i]], N[i]) # observed count is modeled with binomial distribution
  }
})

# parameters to monitor
params <- c("mu_gamma0", "sd_gamma0",
            "mu_alpha0", "sd_alpha0", "mu_alpha0_contrast", "sd_alpha0_contrast",
            "mu_beta0", "sd_beta0", "mu_beta0_contrast", "sd_beta0_contrast",
            "zeta", "xi", "gamma0", "gamma1", "alpha0", "alpha0_contrast", 
            "beta0_contrast", "beta0", "pie_sp", "pie_sp_tc", "N_ds", "N")

# make initial values for MCMC chains
inits <- function(){
  list(
    NG = data$yNG + 1,
    DCLASS = replace(data$DCLASS, is.na(data$DCLASS), 5),
    mu_gamma0 = 5.5, 
    sd_gamma0 = 0.5, 
    gamma1 = rnorm(1, 0, 1),
    mu_alpha0 = rnorm(constants$NREGION, 0, 2), 
    sd_alpha0 = rexp(constants$NREGION, 1), 
    mu_beta0 = rnorm(constants$NREGION, 0, 2), 
    sd_beta0 = rexp(constants$NREGION, 1), 
    gamma0 = rep(5.5, constants$NSPECIES),
    alpha0 = array(rnorm(constants$NSPECIES * constants$NREGION, 0, 2),
                   dim = c(constants$NSPECIES, constants$NREGION)),
    beta0 = array(rnorm(constants$NSPECIES * constants$NREGION, 0, 2),
                  dim = c(constants$NSPECIES, constants$NREGION)),
    zeta =  rgamma(constants$NSPECIES, 4, 2),
    xi = rgamma(constants$NSPECIES, 4, 2), 
    rho_tc = rgamma(constants$NCOUNTS, 2, 2), 
    phi_tc = rgamma(constants$NCOUNTS, 2, 2), 
    N =  data$yN + round(runif(length(data$yN), 1,  6)),
    rho = rgamma(constants$NSURVEYS, 2, 2), 
    phi = rgamma(constants$NDISTANCES, 2, 2))
}

# following code sets up to run the model in parallel, 3 chains
# this is a pretty decent-sized model (took ~ 12 hours on supercomputer)
# so if running on a laptop or regular desktop, you might want to run just 1 chain
# for fewer iterations for it to finish in a reasonable amount of time
nc <- 3
nburn <- 100000
ni <- nburn + 200000
nt <- 200

start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model.code",
                              "inits", 
                              "data", 
                              "constants", 
                              "params", 
                              "nburn", 
                              "ni", 
                              "nt"))

for(j in seq_along(cl)) {
  set.seed(j)
  init <- inits()
  parallel::clusterExport(cl[j], "init")
}

out <- parallel::clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model <- nimble::nimbleModel(code = model.code, 
                               name = "model.code", 
                               constants = constants, 
                               data = data, 
                               inits = init)
  
  Cmodel <- nimble::compileNimble(model)
  modelConf <- nimble::configureMCMC(model)
  modelConf$nimble::addMonitors(params)
  modelMCMC <- nimble::buildMCMC(modelConf)
  CmodelMCMC <- nimble::compileNimble(modelMCMC, project = model)
  out1 <- nimble::runMCMC(CmodelMCMC, 
                          nburnin = nburn, 
                          niter = ni, 
                          thin = nt)
  
  return(as.mcmc(out1))
})
end <- Sys.time()
end - start
parallel::stopCluster(cl)

# package up and save results
setwd(here::here("results"))
save(
  out, 
  model.code, 
  data, 
  constants,
  file = "herbivore_case_study_results_v02.RData"
)