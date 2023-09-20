# 19 September 2023
# Author: Neil Gilbert
# Script to simulate data for non-group-living species
# And updated model code showing a form of the model for non-group-living species

library(tidyverse)
library(nimble)
library(parallel)

sim_icm <- function(
    nsp = 20,          # number of species
    mu_alpha0 = 1.0,   # community average abundance intercept (log scale)
    sigma_alpha0 = 1.5, # standard deviation of abundance intercept among species
    mu_alpha1 = 0.05,    # community average for covariate effect
    sigma_alpha1 = 0.25,   # standard deviation among species for covariate effect
    mu_gamma0 = 5.5,     # community average of scale parameter (detection function) intercept: distance sampling data
    mu_gamma0C = 4.8, # community average of scale parameter (detection function) intercept: count data
    sigma_gamma0 = 0.25, # standard deviation among species for scale parameter intercept: distance sampling data
    sigma_gamma0C = 0.3, # standard deviation among species for scale parameter intercept: count data
    min_zeta = 0.2,      # minimum of uniform distribution from which to draw abundance overdispersion
    max_zeta = 1,      # maximum of uniform distribution from which to draw abundance overdispersion
    nsites = 50,       # number of sites for distance sampling
    nrep = 1,          # number of temporal replicates
    b = 1000,          # distance to which animals are counted
    width = 25,        # width of distance classes
    nsites_tc_fact = 1 # multiplication factor of how much more count data sites there are
){
  
  nsites_tc <- nsites * nsites_tc_fact
  
  # abundance - intercept & covariate coefficient
  alpha0 <- rnorm( nsp, mean = mu_alpha0, sd = sigma_alpha0 )
  alpha1 <- rnorm( nsp, mean = mu_alpha1, sd = sigma_alpha1 )
  
  # intercept for scale parameter for distance sampling
  gamma0 <- rnorm( nsp, mean = mu_gamma0, sd = sigma_gamma0 )
  
  # intercept for scale paramter for count data
  gamma0C <- rnorm( nsp, mean = mu_gamma0C, sd = sigma_gamma0C )
  
  # overdispersion hyperparameter for abundance
  zeta <- runif( nsp, min_zeta, max_zeta )
  
  sp_df <- tibble::tibble(
    sp = 1:nsp,
    alpha0 = alpha0,
    alpha1 = alpha1, 
    gamma0 = gamma0, 
    gamma0C = gamma0C,
    zeta = zeta)
  
  com_truth <- tibble::tribble(
    ~param, ~truth,
    "mu_gamma0", mu_gamma0,
    "sd_gamma0", sigma_gamma0,
    "mu_gamma0C", mu_gamma0C,
    "sd_gamma0C", sigma_gamma0C,
    "mu_alpha0", mu_alpha0,
    "sd_alpha0", sigma_alpha0,
    "mu_alpha1", mu_alpha1,
    "sd_alpha1", sigma_alpha1)
  
  site_covs_ds <- tibble::tibble(
    site = 1:nsites,
    x = runif(nsites, -2, 2)) |> 
    mutate(x = as.numeric(scale(x)))
  
  n_df_ds <- expand.grid(sp = 1:nsp,
                         site =  1:nsites,
                         rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df) |>
    dplyr::full_join(site_covs_ds) |>
    dplyr::rename( xvar = x) |> 
    ( function(x) dplyr::mutate(x,
                                rho = rgamma(nrow(x), zeta, zeta),
                                en = exp( alpha0 + alpha1 * xvar) * rho,
                                n = rpois(nrow(x), en)))()
  
  
  n_vector <- c()
  site_vector <- c()
  rep_vector <- c()
  sp_vector <- c()
  for(i in 1 : nrow( n_df_ds )) {
    if( n_df_ds[[i, "n"]] == 0){
      n_vector <- c(n_vector, 0)
      site_vector <- c(site_vector, n_df_ds[[i, "site"]])
      rep_vector <- c(rep_vector, n_df_ds[[i, "rep"]])
      sp_vector <- c(sp_vector, n_df_ds[[i, "sp"]])
    } else {
      n_vector <- c(n_vector, rep(1, n_df_ds[[i, "n"]]))
      site_vector <- c(site_vector, rep(n_df_ds[[i, "site"]], n_df_ds[[i, "n"]]))
      rep_vector <- c(rep_vector, rep(n_df_ds[[i, "rep"]], n_df_ds[[i, "n"]]))
      sp_vector <- c(sp_vector, rep(n_df_ds[[i, "sp"]], n_df_ds[[i, "n"]]))
    }
  }
  
  # expanded df so each individual can have a distance observation
  n_df_ds_expanded <- tibble::tibble(
    site = site_vector, 
    rep = rep_vector, 
    sp = sp_vector, 
    ind = n_vector) |> # placeholder for saying yes, there is an individual
    dplyr::full_join(n_df_ds)
  
  # assign distances to each individual and simulate observation process, based on distance
  sigma <- exp(sp_df$gamma0)
  data_ds <- NULL
  for( i in 1 : nrow(n_df_ds_expanded) ) {
    if(n_df_ds_expanded[[i, "n"]] == 0){
      data_ds <- tibble::as_tibble(
        rbind(data_ds,
              cbind(
                site = n_df_ds_expanded[[i, "site"]],
                rep = n_df_ds_expanded[[i, "rep"]],
                sp = n_df_ds_expanded[[i, "sp"]],
                ind = n_df_ds_expanded[[i, "ind"]],
                en = n_df_ds_expanded[[i, "en"]],
                n = n_df_ds_expanded[[i, "n"]],
                ind_obs = 0, 
                dclass = NA)))
    } else {
      d <- runif( 1, 0, b) # animals distributed uniformly
      dclass <- d %/% width + 1 # grab the dclass that it falls into
      # detection probability is a function of distance and the scale parameter
      p <- exp( -d * d / (2 * sigma[n_df_ds_expanded[[i, "sp"]]] ^ 2))
      # was or was not the individual observed?
      ind_obs <- rbinom(n_df_ds_expanded[[i, "ind"]], 1, p)
      
      data_ds <- tibble::as_tibble(
        rbind(data_ds,
              cbind(
                site = n_df_ds_expanded[[i, "site"]],
                rep = n_df_ds_expanded[[i, "rep"]],
                sp = n_df_ds_expanded[[i, "sp"]],
                ind = n_df_ds_expanded[[i, "ind"]],
                en = n_df_ds_expanded[[i, "en"]],
                n = n_df_ds_expanded[[i, "n"]],
                ind_obs = ind_obs, 
                dclass = dclass)))
    }
  }
  
  site_covs_c <- tibble::tibble(
    site = 1:nsites_tc,
    x = runif(nsites_tc, -2, 2)) |> 
    mutate(x = as.numeric(scale(x)))
  
  n_df_c <- expand.grid(sp = 1:nsp,
                        site =  1:nsites_tc,
                        rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df) |>
    dplyr::full_join(site_covs_c) |>
    dplyr::rename( xvar = x) |> 
    ( function(x) dplyr::mutate(x,
                                rho = rgamma(nrow(x), zeta, zeta),
                                en = exp( alpha0 + alpha1 * xvar) * rho,
                                n = rpois(nrow(x), en)))()
  
  
  n_vector_c <- c()
  site_vector_c <- c()
  rep_vector_c <- c()
  sp_vector_c <- c()
  for(i in 1 : nrow( n_df_c )) {
    if( n_df_c[[i, "n"]] == 0){
      n_vector_c <- c(n_vector_c, 0)
      site_vector_c <- c(site_vector_c, n_df_c[[i, "site"]])
      rep_vector_c <- c(rep_vector_c, n_df_c[[i, "rep"]])
      sp_vector_c <- c(sp_vector_c, n_df_c[[i, "sp"]])
    } else {
      n_vector_c <- c(n_vector_c, rep(1, n_df_c[[i, "n"]]))
      site_vector_c <- c(site_vector_c, rep(n_df_c[[i, "site"]], n_df_c[[i, "n"]]))
      rep_vector_c <- c(rep_vector_c, rep(n_df_c[[i, "rep"]], n_df_c[[i, "n"]]))
      sp_vector_c <- c(sp_vector_c, rep(n_df_c[[i, "sp"]], n_df_c[[i, "n"]]))
    }
  }
  
  n_df_c_expanded <- tibble::tibble(
    site = site_vector_c, 
    rep = rep_vector_c, 
    sp = sp_vector_c, 
    ind = n_vector_c) |> # ind is just a placeholder - means yes, there is an individual
    dplyr::full_join(n_df_c)
  
  # assign distances to each individual and simulate observation process, based on distance
  sigmaC <- exp(sp_df$gamma0C)
  data_c <- NULL
  for( i in 1 : nrow(n_df_c_expanded) ) {
    if(n_df_c_expanded[[i, "n"]] == 0){
      data_c <- tibble::as_tibble(
        rbind(data_c,
              cbind(
                site = n_df_c_expanded[[i, "site"]],
                rep = n_df_c_expanded[[i, "rep"]],
                sp = n_df_c_expanded[[i, "sp"]],
                ind = n_df_c_expanded[[i, "ind"]],
                en = n_df_c_expanded[[i, "en"]],
                n = n_df_c_expanded[[i, "n"]],
                ind_obs = 0, 
                dclass = NA)))
    } else {
      d <- runif( 1, 0, b) # animals distributed uniformly
      dclass <- d %/% width + 1 # grab the dclass that it falls into
      # detection probability is a function of distance and the scale parameter
      p <- exp( -d * d / (2 * sigmaC[n_df_c_expanded[[i, "sp"]]] ^ 2))
      # was or was not the individual observed?
      ind_obs <- rbinom(n_df_c_expanded[[i, "ind"]], 1, p)
      
      data_c <- tibble::as_tibble(
        rbind(data_c,
              cbind(
                site = n_df_c_expanded[[i, "site"]],
                rep = n_df_c_expanded[[i, "rep"]],
                sp = n_df_c_expanded[[i, "sp"]],
                ind = n_df_c_expanded[[i, "ind"]],
                en = n_df_c_expanded[[i, "en"]],
                n = n_df_c_expanded[[i, "n"]],
                ind_obs = ind_obs, 
                dclass = dclass)))
    }
  }
  
  transect_counts <- sp_df |> 
    dplyr::full_join(data_c) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(sp, site, rep) |> 
    dplyr::mutate(n_obs = sum(ind_obs)) |> 
    dplyr::select(sp, site, rep, dclass, ind_obs, en, n, n_obs) |> 
    dplyr::distinct() |> 
    dplyr::group_by(sp, site, rep) |> 
    dplyr::summarise( 
      true_n = unique(n),
      count = unique(n_obs)) |> 
    dplyr::full_join(site_covs_c) |> 
    dplyr::select(sp, site, rep, true_n, count, x_tc = x)
  
  ds_data_n <- data_ds |> 
    group_by( sp, site, rep ) |> 
    dplyr::summarise( true_n = unique(n), 
                      count = sum(ind_obs)) |> 
    full_join(site_covs_ds)
  
  ds_data_ds <- data_ds |> 
    filter(ind_obs == 1) |> 
    dplyr::select(site, rep, sp, dclass)
  
  data <- list(
    MIDPOINT_DS = seq(from = 12.5, to = 987.5, by = 25),
    MIDPOINT_C = seq(from = 12.5, to = 987.5, by = 25),
    DCLASS = ds_data_ds$dclass,  
    V_DS = 25, 
    V_C = 25,
    B_DS = 1000,
    B_C = 1000,
    yN_DS = ds_data_n$count,
    HAB_DS = ds_data_n$x,
    HAB_C = transect_counts$x_tc,
    yN_C = transect_counts$count,
    true_n_ds = ds_data_n$true_n,
    true_n_c = transect_counts$true_n)
  
  constants <- list(
    N_SPECIES = length(unique(transect_counts$sp)),
    N_BINS_DS = length(data$MIDPOINT_DS),
    N_BINS_C = length(data$MIDPOINT_C),
    N_DS_D = nrow(ds_data_ds), 
    SP_DS_D = ds_data_ds$sp, 
    SP_DS_N = ds_data_n$sp,
    N_DS_N = nrow(ds_data_n),  
    N_COUNTS = nrow(transect_counts), 
    SP_C = transect_counts$sp)
  
  sp_info <- ds_data_n |> 
    dplyr::group_by(sp) |> 
    dplyr::summarise(totDS = sum(true_n)) |> 
    dplyr::full_join(dplyr::summarise(dplyr::group_by(transect_counts, sp), totTC = sum(true_n))) |> 
    dplyr::full_join(sp_df)
  
  return(list(data = data,
              constants = constants,
              sp_info = sp_info,
              com_truth = com_truth))
}

#### Model code ####
model.code <- nimble::nimbleCode({
  
  mu_gamma0_ds ~ dnorm(5.5, sd = 0.5) 
  sd_gamma0_ds ~ dexp(1)
  
  mu_gamma0_c ~ dnorm(5.5, sd = 0.5)
  sd_gamma0_c ~ dexp(1)
  
  mu_alpha0 ~ dnorm(0, sd = 2) 
  sd_alpha0 ~ dexp(1)          
  
  mu_alpha1 ~ dnorm(0, sd = 2)
  sd_alpha1 ~ dexp(1)
  
  for ( s in 1:N_SPECIES ) {
    gamma0_ds[s] ~ dnorm( mu_gamma0_ds, sd = sd_gamma0_ds )
    gamma0_c[s] ~ dnorm( mu_gamma0_c, sd = sd_gamma0_c )
    omega_ds[s] <- exp( gamma0_ds[s] )
    omega_c[s] <- exp( gamma0_c[s] )
    alpha0[s] ~ dnorm( mu_alpha0, sd = sd_alpha0 )
    alpha1[s] ~ dnorm( mu_alpha1, sd = sd_alpha1 )
    zeta[s] ~ dgamma(4, 2)    
    pie_sp_ds[s] <- sum( pie_ds[1:N_BINS_DS, s] )
    for( k in 1:N_BINS_DS){
      log(g_ds[k, s]) <- -MIDPOINT_DS[k] * MIDPOINT_DS[k]/(2 * omega_ds[s] * omega_ds[s])
      pie_ds[k, s] <- g_ds[k,s] * (V_DS/B_DS)
      pie_cell_ds[k, s] <- pie_ds[k, s] / pie_sp_ds[s]
    }
    pie_sp_c[s] <- sum(pie_c[1:N_BINS_C, s])
    for (k in 1:N_BINS_C ) {
      log(g_c[k,s]) <- -MIDPOINT_C[k] * MIDPOINT_C[k]/(2 * omega_c[s] * omega_c[s])
      pie_c[k,s] <- g_c[k,s] * (V_C/B_C)
    }
  }
  for( i in 1:N_DS_N ) {
    log( lambda_ds[i] ) <- alpha0[ SP_DS_N[i]] + alpha1[ SP_DS_N[i]] * HAB_DS[i] 
    rho_ds[i] ~ dgamma( zeta[SP_DS_N[i]], zeta[SP_DS_N[i]] )
    N_DS[i] ~ dpois( lambda_ds[i] * rho_ds[i] )
    yN_DS[i] ~ dbin( pie_sp_ds[SP_DS_N[i]], N_DS[i] )
  }
  for (i in 1:N_DS_D ) {
    DCLASS[i] ~ dcat(pie_cell_ds[1:N_BINS_DS, SP_DS_D[i] ] )
  }
  for(i in 1:N_COUNTS) {
    log(lambda_c[i]) <- alpha0[SP_C[i]] + alpha1[ SP_C[i]] * HAB_C[i]
    rho_c[i] ~ dgamma(zeta[SP_C[i]], zeta[SP_C[i]])
    N_C[i] ~ dpois( lambda_c[i] * rho_c[i] )
    yN_C[i] ~ dbin( pie_sp_c[SP_C[i]], N_C[i] )
  }
})

#### MCMC settings & simulation scenarios ####

params <- c(
  "mu_gamma0_ds", 
  "sd_gamma0_ds", 
  "mu_gamma0_c", 
  "sd_gamma0_c", 
  "mu_alpha0", 
  "sd_alpha0", 
  "mu_alpha1", 
  "sd_alpha1",
  "gamma0_ds", 
  "gamma0_c",
  "alpha0",
  "alpha1",
  "zeta",
  "pie_sp_ds",
  "pie_sp_c",
  "N_DS",
  "N_C")

simdat <- sim_icm()
attach(simdat)
data <- simdat$data
constants <- simdat$constants

make_inits <- function(data, constants){
  
  mu_gamma0_ds_st <- rnorm(1, 5.5, 0.1)
  sd_gamma0_ds_st <- runif(1, 0.1, 0.5)
  
  mu_gamma0_c_st <- rnorm(1, 5.5, 0.1)
  sd_gamma0_c_st <- runif(1, 0.1, 0.5)
  
  mu_alpha0_st <- rnorm(1, 0,  2)
  sd_alpha0_st <- runif(1, 1, 2)
  
  mu_alpha1_st <- rnorm(1, 0, 2)
  sd_alpha1_st <- runif(1, 0.3, 0.75)
  
  zeta_st <- rgamma(constants$N_SPECIES, 4, 2)
  
  alpha0_st <- rnorm( constants$N_SPECIES, mu_alpha0_st, sd_alpha0_st )
  alpha1_st <- rnorm( constants$N_SPECIES, mu_alpha1_st, sd_alpha1_st )
  gamma0_ds_st <- rnorm( constants$N_SPECIES, mu_gamma0_ds_st, sd_gamma0_ds_st )
  gamma0_c_st <- rnorm( constants$N_SPECIES, mu_gamma0_c_st, sd_gamma0_c_st )
  
  inits <- list(
    mu_gamma0_ds = mu_gamma0_ds_st,
    sd_gamma0_ds = sd_gamma0_ds_st,
    mu_gamma0_c = mu_gamma0_c_st,
    sd_gamma0_c = sd_gamma0_c_st,
    mu_alpha0 = mu_alpha0_st, 
    sd_alpha0 = sd_alpha0_st,
    mu_alpha1 = mu_alpha1_st,
    sd_alpha1 = sd_alpha1_st,
    alpha0 = alpha0_st,
    alpha1 = alpha1_st,
    gamma0_ds = gamma0_ds_st,
    gamma0_c = gamma0_c_st,
    zeta = zeta_st,
    rho_ds = rgamma(constants$N_DS_N, 2, 2), 
    N_DS = data$yN_DS + 1,
    DCLASS = replace(data$DCLASS, is.na(data$DCLASS), 5),
    rho_c = rgamma(constants$N_COUNTS, 2, 2), 
    N_C =  data$yN_C + round(runif(length(data$yN_C), 1,  6)))
  return(inits)
}

# Ran in several hours on HPCC and converged
nc <- 3
nburn <- 100000
ni <- nburn + 200000
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
  set.seed(NULL)
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
time <- difftime(end, start, units = "hours")
parallel::stopCluster(cl)

save( 
  model.code, 
  data, 
  constants, 
  sp_info, 
  file = "simulated_no_gs_v01.RData")