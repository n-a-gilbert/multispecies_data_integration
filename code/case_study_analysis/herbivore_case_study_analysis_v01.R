library(here)
library(nimble)
library(parallel)
library(tidyverse)

setwd(here::here("data"))

load("distance_sampling_data_v01.RData")
load("count_data_v01.RData")

# format data to get ds counts
ng_data <- final2 |>
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

det_attributes <- final2 |>
  dplyr::select( sp, sp_name ) |>
  dplyr::distinct() |>
  dplyr::full_join( sp_mass ) |>
  dplyr::mutate( mass = base::as.numeric( base::scale( base::log( mass ) ) ) )

#package data up for nimble model
data <- list(
  MIDPOINT = mdpt,  # midpoint of distance class bins
  DCLASS = gs_data$dclass, # distance class of group size counts (vector)
  B_DS = b,  # geometry for distance function
  B_TC = 100, # width of count transect
  V = v,   # geometry for distance function 
  yN_DS = ng_data$yN_DS, # total abundance on DS transect
  yN_TC = transect_data$count, # transect counts
  OFFSET_DS = ng_data$area, # area offset for distance sampling transects
  OFFSET_TC = transect_data$area, # area offset for count transects
  MASS = det_attributes$mass ) # species body mass

constants <- list(
  NSPECIES = length(unique(final2$sp)), # number of species (11)
  NBINS = length(mdpt), # number of distance bins 
  NBINS_C = 4L,
  NDISTANCES = nrow(gs_data), # number of distance observations
  NSURVEYS = nrow(ng_data),  # number of ds surveys
  NCOUNTS = nrow(transect_data), # number of transect counts
  SP_GS = gs_data$sp, # species ID, used for indexing distance data
  SP_NG = ng_data$sp,# species ID, used for indexing counts from distance sampling data
  SP_TC = transect_data$sp, # species ID, used for indexing transect counts
  REGION_NG = ng_data$region + 1, # Region - Mara (1) or Talek (2)
  REGION_TC = transect_data$region + 1, # Region  - Mara (1) or Talek (2)
  REGION_GS = gs_data$region + 1,
  NREGION = 2) # number of regions

# NIMBLE code
model.code <- nimble::nimbleCode({
  gamma1 ~ dnorm(0, sd = 0.1) # effect of mass on detectabilitiy
  for( j in 1:NREGION) {
    mu_alpha0[j] ~ dnorm(0, sd = 2)       # community average of abundance intercept
    sd_alpha0[j] ~ dexp(1)               # sd among species in abundance intercept
    mu_gamma0[j] ~ dnorm(5.5, sd = 0.5)  # community average of detection function scale parameter intercept
    sd_gamma0[j] ~ dexp(1)               # sd among species for detection function scale parameter intercept
  }
  for (s in 1:NSPECIES) {
    zeta[s] ~ dgamma(4, 2) # overdispersion hyperparameter
    for (j in 1:NREGION) {
      alpha0[s,j] ~ dnorm( mu_alpha0[j], sd = sd_alpha0[j] )  # species x region abundance intercept
      gamma0[s,j] ~ dnorm( mu_gamma0[j], sd = sd_gamma0[j] )  # species x region scale parameter intercept
      omega[s,j] <- exp( gamma0[s,j] + gamma1 * MASS[s] )    # scale parameter
      pie_sp[s,j] <- sum(pie_ds[1:NBINS, s, j])              # transect-level detection probability (DS)
      pie_sp_tc[s,j] <- sum(pie_c[1:NBINS_C, s, j])          # transect-level detection probability (counts)
      for (k in 1:NBINS) {
        log(g[k, s, j]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s,j] * omega[s,j]) # half-normal detection function
        pie_ds[k, s, j] <- g[k, s, j] * (V/B_DS)                                    # bin-level detection probability (DS)
        pie_cell[k, s, j] <- pie_ds[k, s, j]/pie_sp[s,j]                            # cell probabilities
      }
      for(k in 1:NBINS_C) {
        log( g_c[k,s,j] ) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s,j] * omega[s,j]) # half-normal detection function (counts)
        pie_c[k, s,j] <- g_c[k,s,j] * (V/B_TC)                                        # bin-level detection probability (counts)
      }
    }
  }
  for (i in 1:NSURVEYS) {
    log(lambda[i]) <- alpha0[SP_NG[i], REGION_NG[i]] + log(OFFSET_DS[i]) # expected abundance
    rho[i] ~ dgamma(zeta[SP_NG[i]], zeta[SP_NG[i]])                      # overdispersion parameter
    N_DS[i] ~ dpois( lambda[i] * rho[i] )                               # latent abundance
    yN_DS[i] ~ dbin( pie_sp[SP_NG[i], REGION_NG[i]], N_DS[i] )          # observed count
  }
  
  # loop through the distance observations
  for (i in 1:NDISTANCES) {
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i], REGION_GS[i]] ) # distance class modeled as draw from categorical distribution
  }
  
  # count submodel
  for (i in 1:NCOUNTS) {
    log(lambda_tc[i] ) <- alpha0[SP_TC[i], REGION_TC[i]] + log(OFFSET_TC[i]) # expected abundance
    rho_tc[i] ~ dgamma(zeta[SP_TC[i]], zeta[SP_TC[i]])                       # overdispersion
    N_TC[i] ~ dpois( lambda_tc[i] * rho_tc[i] )                               # latent abundance
    yN_TC[i] ~ dbin(pie_sp_tc[SP_TC[i], REGION_TC[i]], N_TC[i]) # observed count is modeled with binomial distribution
  }
})

# parameters to monitor
params <- c("mu_gamma0", "sd_gamma0",
            "mu_alpha0", "sd_alpha0", 
            "zeta",
            "gamma0",
            "gamma1",
            "alpha0",
            "pie_sp",
            "pie_sp_tc",
            "N_DS",
            "N_TC")

make_inits <- function(data, constants){
  mu_gamma0_st <- rnorm(2, 5.5, 0.1)
  sd_gamma0_st <- runif(2, 0.5, 1.25)
  gamma1_st <- rnorm(1, 0, 0.1)
  mu_alpha0_st <- rnorm(constants$NREGION, 0,  2)
  sd_alpha0_st <- runif(constants$NREGION, 1, 2)
  zeta_st <- rgamma(constants$NSPECIES, 4, 2)
  zalpha0_st <- array( rnorm (constants$NSPECIES * constants$NREGION, 0, 0.5), 
                       dim = c(constants$NSPECIES, constants$NREGION))
  gamma0_st <-  array( rnorm (constants$NSPECIES * constants$NREGION, 5.5, 0.25), 
                       dim = c(constants$NSPECIES, constants$NREGION))
  alpha0_st <- array(NA, dim = c(constants$NSPECIES, constants$NREGION))
  for( s in 1:constants$NSPECIES ){
    for( j in 1:constants$NREGION){
      alpha0_st[s, j] <- mu_alpha0_st[j] + sd_alpha0_st[j]*zalpha0_st[s,j]
    }
  }
  
  inits <- list(
    mu_gamma0 = mu_gamma0_st,
    sd_gamma0 = sd_gamma0_st,
    gamma1 = gamma1_st,
    mu_alpha0 = mu_alpha0_st, 
    sd_alpha0 = sd_alpha0_st,
    alpha0 = alpha0_st,
    gamma0 = gamma0_st,
    zeta = zeta_st,
    rho = rgamma(constants$NSURVEYS, 2, 2), 
    N_DS = data$yN_DS + 1,
    DCLASS = replace(data$DCLASS, is.na(data$DCLASS), 5),
    rho_tc = rgamma(constants$NCOUNTS, 2, 2), 
    N_TC =  data$yN_TC + round(runif(length(data$yN_TC), 1,  6)))
  return(inits)
}

# following code sets up to run the model in parallel, 3 chains
# this is a pretty decent-sized model (took ~ 18 hours on supercomputer)
# so if running on a laptop or regular desktop, you might want to run just 1 chain
# for fewer iterations for it to finish in a reasonable amount of time
nc <- 3
nburn <- 300000
ni <- nburn + 600000
nt <- 200


start <- Sys.time()
cl <- parallel::makeCluster(nc)

parallel::clusterExport(cl, c("model.code",
                              "make_inits", 
                              "data", 
                              "constants", 
                              "params", 
                              "nburn", 
                              "ni", 
                              "nt"))

for(j in seq_along(cl)) {
  set.seed(j)
  init <- make_inits(data, constants)
  parallel::clusterExport(cl[j], "init")
}

out <- parallel::clusterEvalQ(cl, {
  library(nimble)
  library(coda)
  
  model <- nimbleModel(code = model.code, 
                       name = "model.code", 
                       constants = constants, 
                       data = data, 
                       inits = init)
  
  Cmodel <- compileNimble(model)
  modelConf <- configureMCMC(model)
  modelConf$addMonitors(params)
  modelMCMC <- buildMCMC(modelConf)
  CmodelMCMC <- compileNimble(modelMCMC, project = model)
  out1 <- runMCMC(CmodelMCMC, 
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
  file = "herbivore_case_study_results_v01.RData"
)