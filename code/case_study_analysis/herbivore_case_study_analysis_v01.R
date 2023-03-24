library(here)
library(nimble)
library(parallel)
library(tidyverse)

setwd(here::here("data"))

load("distance_sampling_data_v01.RData")
load("count_data_v01.RData")

ng_data <- final2 %>%
  dplyr::select( sp, site, rep, ng, gs, area, region ) %>% 
  dplyr::filter( !base::is.na( ng ) ) %>% 
  dplyr::group_by( sp, site, rep, area, region ) %>% 
  dplyr::summarise( ng = base::unique( ng ),
                    yN_DS = base::sum(gs, na.rm = TRUE ) ) %>% 
  dplyr::arrange( sp, site, rep )

gs_data <- final2 %>% 
  dplyr::filter( !base::is.na( gs ) ) %>% 
  dplyr::arrange( sp, site, rep )

# from PanTHERIA
sp_mass <- tibble(
  sp_name = c("buffalo", "eland", "elephant", "giraffe", "grants",
              "hartebeest", "impala",  "thomsons", "topi",
              "warthog", "waterbuck"),
  mass = c(592666, 562592.7, 3824540, 964654.7, 55464.46,
           160937.9,  52591.69, 22907.43, 136000.3,
           82499.99, 204393.5))

det_attributes <- final2 %>%
  dplyr::select( sp, sp_name ) %>% 
  dplyr::distinct(.) %>% 
  dplyr::full_join( sp_mass ) %>% 
  dplyr::mutate( mass = base::as.numeric( base::scale( base::log( mass ) ) ) )

data <- list(
  MIDPOINT = mdpt,  # midpoint of distance class bins
  yGS = gs_data$gs, # group size counts (vector)
  DCLASS = gs_data$dclass, # distance class of group size counts (vector)
  B_DS = b,  # geometry for distance function
  B_TC = 100, # width of count transect
  V = v,   # geometry for distance function 
  yNG = ng_data$ng, # number of groups per site x rep x species
  yN_DS = ng_data$yN_DS,
  yN = transect_data$count, # transect counts
  OFFSET_DS = ng_data$area, 
  OFFSET_TC = transect_data$area,
  MASS = det_attributes$mass )

constants <- list(
  NSPECIES = length(unique(final2$sp)), # number of species (11)
  NBINS = length(mdpt), # number of distance bins 
  NDISTANCES = nrow(gs_data), # number of group size observations
  NSURVEYS = nrow(ng_data),  # number of number-of-group observations
  NCOUNTS = nrow(transect_data), # number of transect counts
  SP_GS = gs_data$sp, # species ID, used for indexing group size data
  SP_NG = ng_data$sp,# species ID, used for indexing number-of-groups data
  SP_TC = transect_data$sp, # species ID, used for indexing transect counts
  REGION_GS = gs_data$region + 1,
  REGION_NG = ng_data$region + 1,
  REGION_TC = transect_data$region + 1,
  NREGION = 2)

model.code <- nimble::nimbleCode({
  mu_gamma0 ~ dnorm(5.5, sd = 0.5)
  sd_gamma0 ~ dexp(1)
  gamma1 ~ dnorm(0, sd = 2)
  for( j in 1:NREGION) {
    mu_alpha0[j] ~ dnorm(0, sd = 2)
    sd_alpha0[j] ~ dexp(1)
    mu_beta0[j]  ~ dnorm(0, sd = 2)
    sd_beta0[j]  ~ dexp(1)
  }
  mu_alpha0_contrast <- mu_alpha0[2] - mu_alpha0[1]
  sd_alpha0_contrast <- sd_alpha0[2] - sd_alpha0[1]
  mu_beta0_contrast <- mu_beta0[2] - mu_beta0[1]
  sd_beta0_contrast <- sd_beta0[2] - sd_beta0[1]
  for (s in 1:NSPECIES) {
    for (j in 1:NREGION) {
      alpha0[s, j] ~ dnorm(mu_alpha0[j], sd = sd_alpha0[j])
      beta0[s, j] ~ dnorm(mu_beta0[j], sd = sd_beta0[j])
    }
    alpha0_contrast[s] <- alpha0[s, 2] - alpha0[s, 1]
    beta0_contrast[s] <- beta0[s, 2] - beta0[s, 1]
    gamma0[s] ~ dnorm(mu_gamma0, sd = sd_gamma0)
    omega[s] <- exp(gamma0[s] + gamma1 * MASS[s])
    pie_sp[s] <- sum(pie_ds[1:NBINS, s])
    pie_sp_tc[s] <- sum(pie_tc[1:4, s])
    zeta[s] ~ dgamma(4, 2)
    xi[s] ~ dgamma(4, 2)
    for (k in 1:NBINS) {
      log(g[k, s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s] * omega[s])
      pie_ds[k, s] <- g[k, s] * (V/B_DS)
      pie_tc[k, s] <- g[k, s] * (V/B_TC)
      pie_cell[k, s] <- pie_ds[k, s]/pie_sp[s]
    }
  }
  for (i in 1:NSURVEYS) {
    log(lambda[i]) <- alpha0[SP_NG[i], REGION_NG[i]] + log(OFFSET_DS[i])
    rho[i] ~ dgamma(zeta[SP_NG[i]], zeta[SP_NG[i]])
    NG[i] ~ dpois(lambda[i] * rho[i])
    yNG[i] ~ dbin(pie_sp[SP_NG[i]], NG[i])
    NG_missed[i] <- NG[i] - yNG[i]
    log( kappa_ng[i] ) <- beta0[SP_NG[i], REGION_NG[i]]
    GS_missed[i] ~ T( dpois( kappa_ng[i] ), 1, )
    N_missed[i] <- NG_missed[i] * GS_missed[i]
    N_ds[i] <- yN_DS[i] + N_missed[i]
  }
  for (i in 1:NDISTANCES) {
    log(kappa[i]) <- beta0[SP_GS[i], REGION_GS[i]]
    phi[i] ~ dgamma(xi[SP_GS[i]], xi[SP_GS[i]])
    yGS[i] ~ T(dpois(kappa[i] * phi[i]), 1, )
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i]])
  }
  for (i in 1:NCOUNTS) {
    rho_tc[i] ~ dgamma(zeta[SP_TC[i]], zeta[SP_TC[i]])
    phi_tc[i] ~ dgamma(xi[SP_TC[i]], xi[SP_TC[i]])
    NG_tc[i] <- exp(alpha0[SP_TC[i], REGION_TC[i]] + log(OFFSET_TC[i])) * rho_tc[i]
    GS_tc[i] <- exp(beta0[SP_TC[i], REGION_TC[i]]) * phi_tc[i]
    N[i] ~ dpois(NG_tc[i] * GS_tc[i])
    yN[i] ~ dbin(pie_sp_tc[SP_TC[i]], N[i])
  }
})

params <- c("mu_gamma0", "sd_gamma0",
            "mu_alpha0", "sd_alpha0", "mu_alpha0_contrast", "sd_alpha0_contrast",
            "mu_beta0", "sd_beta0", "mu_beta0_contrast", "sd_beta0_contrast",
            "zeta", "xi", "gamma0", "gamma1", "alpha0", "alpha0_contrast", 
            "beta0_contrast", "beta0", "pie_sp", "pie_sp_tc", "N_ds", "N")

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

nc <- 3
nburn <- 100000
ni <- nburn + 200000
nt <- 200

start <- Sys.time()
cl <- makeCluster(nc)

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
  clusterExport(cl[j], "init")
}

out <- clusterEvalQ(cl, {
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
stopCluster(cl)

setwd(here::here("results"))
save(
  out, 
  model.code, 
  data, 
  constants,
  file = "herbivore_case_study_results_v02.RData"
)