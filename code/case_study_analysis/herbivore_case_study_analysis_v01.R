# library(here)
library(nimble)
library(parallel)
library(tidyverse)

# setwd(here::here("data"))

load("distance_sampling_data_v01.RData")
load("count_data_v01.RData")

# format data to get number-of-group data
ng_data <- final2 %>%
  dplyr::select( sp, site, rep, ng, gs, area, region ) %>%
  dplyr::filter( !base::is.na( ng ) ) %>%
  dplyr::group_by( sp, site, rep, area, region ) %>%
  dplyr::summarise( ng = base::unique( ng ),
                    yN_DS = base::sum(gs, na.rm = TRUE ) ) %>%
  dplyr::arrange( sp, site, rep )

# format data to get group size data
gs_data <- final2 %>%
  dplyr::filter( !base::is.na( gs ) ) %>%
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
  dplyr::select( sp, sp_name ) %>%
  dplyr::distinct() %>%
  dplyr::full_join( sp_mass ) %>%
  dplyr::mutate( mass = base::as.numeric( base::scale( base::log( mass ) ) ) )

#package data up for nimble model
data <- list(
  MIDPOINT = mdpt,  # midpoint of distance class bins
  # yGS = gs_data$gs, # group size counts (vector)
  DCLASS = gs_data$dclass, # distance class of group size counts (vector)
  B_DS = b,  # geometry for distance function
  B_TC = 100, # width of count transect
  V = v,   # geometry for distance function 
  # yNG = ng_data$ng, # number of groups per site x rep x species
  yN_DS = ng_data$yN_DS, # total abundance on DS transect
  yN_TC = transect_data$count, # transect counts
  OFFSET_DS = ng_data$area, # area offset for distance sampling transects
  OFFSET_TC = transect_data$area, # area offset for count transects
  MASS = det_attributes$mass ) # species body mass

constants <- list(
  NSPECIES = length(unique(final2$sp)), # number of species (11)
  NBINS = length(mdpt), # number of distance bins 
  NBINS_C = 4L,
  NDISTANCES = nrow(gs_data), # number of group size observations
  NSURVEYS = nrow(ng_data),  # number of number-of-group observations
  NCOUNTS = nrow(transect_data), # number of transect counts
  SP_GS = gs_data$sp, # species ID, used for indexing group size data
  SP_NG = ng_data$sp,# species ID, used for indexing number-of-groups data
  SP_TC = transect_data$sp, # species ID, used for indexing transect counts
  # REGION_GS = gs_data$region + 1, # Region = Mara (1) or Talek (2)
  REGION_NG = ng_data$region + 1, # Region - Mara (1) or Talek (2)
  REGION_TC = transect_data$region + 1, # Region  - Mara (1) or Talek (2)
  REGION_GS = gs_data$region + 1,
  NREGION = 2) # number of regions

# NIMBLE code
model.code <- nimble::nimbleCode({
  gamma1 ~ dnorm(0, sd = 0.1)
  for( j in 1:NREGION) {
    mu_alpha0[j] ~ dnorm(0, sd = 2)
    sd_alpha0[j] ~ dexp(1)
    mu_gamma0[j] ~ dnorm(5.5, sd = 0.5) 
    sd_gamma0[j] ~ dexp(1)
  }
  for (s in 1:NSPECIES) {
    zeta[s] ~ dgamma(4, 2)
    for (j in 1:NREGION) {
      alpha0[s,j] ~ dnorm( mu_alpha0[j], sd = sd_alpha0[j] )
      gamma0[s,j] ~ dnorm( mu_gamma0[j], sd = sd_gamma0[j] ) 
      omega[s,j] <- exp( gamma0[s,j] + gamma1 * MASS[s] )
      pie_sp[s,j] <- sum(pie_ds[1:NBINS, s, j])
      pie_sp_tc[s,j] <- sum(pie_c[1:NBINS_C, s, j]) 
      for (k in 1:NBINS) {
        log(g[k, s, j]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s,j] * omega[s,j])
        pie_ds[k, s, j] <- g[k, s, j] * (V/B_DS)
        pie_cell[k, s, j] <- pie_ds[k, s, j]/pie_sp[s,j] 
      }
      for(k in 1:NBINS_C) {
        log( g_c[k,s,j] ) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s,j] * omega[s,j])
        pie_c[k, s,j] <- g_c[k,s,j] * (V/B_TC)
      }
    }
  }
  for (i in 1:NSURVEYS) {
    log(lambda[i]) <- alpha0[SP_NG[i], REGION_NG[i]] + log(OFFSET_DS[i])
    rho[i] ~ dgamma(zeta[SP_NG[i]], zeta[SP_NG[i]]) 
    N_DS[i] ~ dpois( lambda[i] * rho[i] )
    yN_DS[i] ~ dbin( pie_sp[SP_NG[i], REGION_NG[i]], N_DS[i] )
  }
  
  # loop through the distance observations
  for (i in 1:NDISTANCES) {
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i], REGION_GS[i]] ) # distance class modeled as draw from categorical distribution
  }
  
  # count submodel
  for (i in 1:NCOUNTS) {
    log(lambda_tc[i] ) <- alpha0[SP_TC[i], REGION_TC[i]] + log(OFFSET_TC[i])
    rho_tc[i] ~ dgamma(zeta[SP_TC[i]], zeta[SP_TC[i]]) 
    N_TC[i] ~ dpois( lambda_tc[i] * rho_tc[i] ) 
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
  
  # mu_gamma0C_st <- rnorm(1, 5.5, 0.1)
  # sd_gamma0C_st <- runif(1, 0.5, 1.25)
  
  gamma1_st <- rnorm(1, 0, 0.1)
  
  mu_alpha0_st <- rnorm(constants$NREGION, 0,  2)
  sd_alpha0_st <- runif(constants$NREGION, 1, 2)
  
  # mu_beta0_st <- rnorm(constants$NREGION, 0, 2)
  # sd_beta0_st <- runif(constants$NREGION, 0.5, 1)
  
  zeta_st <- rgamma(constants$NSPECIES, 4, 2)
  # xi_st <- rgamma(constants$NSPECIES, 4, 2)
  
  zalpha0_st <- array( rnorm (constants$NSPECIES * constants$NREGION, 0, 0.5), 
                       dim = c(constants$NSPECIES, constants$NREGION))
  
  # zbeta0_st <- array( rnorm (constants$NSPECIES * constants$NREGION, 0, 0.5), 
                      # dim = c(constants$NSPECIES, constants$NREGION))
  
  zgamma0_st <- rnorm(constants$NSPECIES, 0, 0.5)
  # zgamma0C_st <- rnorm(constants$NSPECIES, 0, 0.5)
  gamma0_st <-  array( rnorm (constants$NSPECIES * constants$NREGION, 5.5, 0.25), 
                       dim = c(constants$NSPECIES, constants$NREGION))
  
  # gamma0_st <- gamma0C_st <- numeric(length = constants$NSPECIES)
  # alpha0_st <- beta0_st <- array(NA, dim = c(constants$NSPECIES, constants$NREGION))
  alpha0_st <- array(NA, dim = c(constants$NSPECIES, constants$NREGION))
  
  
  for( s in 1:constants$NSPECIES ){
    # gamma0_st[s] <- mu_gamma0_st + sd_gamma0_st * zgamma0_st[s]
    # gamma0C_st[s] <- mu_gamma0C_st + sd_gamma0C_st*zgamma0C_st[s]
    for( j in 1:constants$NREGION){
      alpha0_st[s, j] <- mu_alpha0_st[j] + sd_alpha0_st[j]*zalpha0_st[s,j]
      # beta0_st[s, j] <- mu_beta0_st[j] + sd_beta0_st[j]*zbeta0_st[s,j]
    }
  }
  
  inits <- list(
    mu_gamma0 = mu_gamma0_st,
    sd_gamma0 = sd_gamma0_st,
    # mu_gamma0C = mu_gamma0C_st,
    # sd_gamma0C = sd_gamma0C_st,
    gamma1 = gamma1_st,
    mu_alpha0 = mu_alpha0_st, 
    sd_alpha0 = sd_alpha0_st,
    # mu_beta0 = mu_beta0_st, 
    # sd_beta0 = sd_beta0_st,
    # zalpha0 = zalpha0_st,
    alpha0 = alpha0_st,
    # zbeta0 = zbeta0_st,
    # beta0 = beta0_st,
    # zgamma0 = zgamma0_st,
    gamma0 = gamma0_st,
    # zgamma0C = zgamma0C_st,
    # gamma0C = gamma0C_st,
    zeta = zeta_st,
    # xi = xi_st,
    rho = rgamma(constants$NSURVEYS, 2, 2), 
    N_DS = data$yN_DS + 1,
    # N_missed = ( data$yNG + 1 )  * sample(1:5, constants$NSURVEYS, replace = TRUE),
    # phi = rgamma(constants$NDISTANCES, 2, 2),
    DCLASS = replace(data$DCLASS, is.na(data$DCLASS), 5),
    rho_tc = rgamma(constants$NCOUNTS, 2, 2), 
    # phi_tc = rgamma(constants$NCOUNTS, 2, 2), 
    N_TC =  data$yN_TC + round(runif(length(data$yN_TC), 1,  6)))
  
  return(inits)
  
}

str( make_inits(data, constants))

model <- nimbleModel(code = model.code,
                     name = "model.code",
                     constants = constants,
                     data = data,
                     inits = make_inits(data, constants))

# following code sets up to run the model in parallel, 3 chains
# this is a pretty decent-sized model (took ~ 12 hours on supercomputer)
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
# setwd(here::here("results"))
save(
  out, 
  model.code, 
  data, 
  constants,
  file = "herbivore_case_study_results_noGS_v03.RData"
)