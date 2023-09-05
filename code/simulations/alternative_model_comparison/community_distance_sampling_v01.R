# 15 March 2023
# Neil Gilbert
# This script runs simulations for a hierarchical community model for distance sampling data only

library(nimble)
library(parallel)
library(tidyverse)
library(extraDistr)
library(MCMCvis)
library(here)

sim_icm <- function(
    nsp = 20,          # number of species
    mu_alpha0 = -0.65,
    sigma_alpha0 = 1.7,
    mu_alpha1 = 0.05,    # community average for covariate effect
    sigma_alpha1 = 0.25,   # standard deviation among species for covariate effect
    mu_beta0 = 1.9,         # average group size for the community
    sigma_beta0 = 1.25,      # standard deviation among species for group size intercept
    mu_gamma0 = 5.5,
    sigma_gamma0 = 0.25,
    min_zeta = 0.2,      # minimum of uniform distribution from which to draw number of groups overdispersion
    max_zeta = 1,      # maximum of uniform distribution from which to draw number of groups overdispersion
    min_xi = 0.25,      # minimum of uniform distribution from which to draw group size overdispersion
    max_xi = 1.8,      # maximum of uniform distribution from which to draw group size overdispersion
    nsites = 50,       # number of sites for distance sampling
    nrep = 1,          # number of temporal replicates
    b = 1000,          # distance to which animals are counted
    width = 25,        # width of distance classes
    nsites_tc_fact = 1,# multiplication factor of how much more count data sites there are
    p_bias = 0        # should count detection probability be same (0), 20% lower (-1), or 20% higher (1) than DS
){
  
  nsites_tc <- nsites * nsites_tc_fact
  
  # number of groups - intercept & covariate coefficient
  alpha0 <- rnorm( nsp, mean = mu_alpha0, sd = sigma_alpha0 )
  alpha1 <- rnorm( nsp, mean = mu_alpha1, sd = sigma_alpha1 )
  
  # group size - intercept
  beta0 <- rnorm( nsp, mean = mu_beta0, sd = sigma_beta0 )
  
  # intercept for scale parameter
  gamma0 <- rnorm( nsp, mean = mu_gamma0, sd = sigma_gamma0 )
  
  # overdispersion hyperparameter for number-of-groups
  zeta <- runif( nsp, min_zeta, max_zeta )
  
  # overdispersion hyperparameter for group size
  xi <- runif( nsp, min_xi, max_xi )
  
  sp_df <- tibble::tibble(
    sp = 1:nsp,
    alpha0 = alpha0,
    alpha1 = alpha1, 
    beta0 = beta0,
    gamma0 = gamma0, 
    zeta = zeta, 
    xi = xi)
  
  site_covs <- tibble::tibble(
    site = 1:nsites,
    x = runif(nsites, -2, 2)) |> 
    mutate(x = as.numeric(scale(x)))
  
  ng_df <- expand.grid(sp = 1:nsp,
                       site =  1:nsites,
                       rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df) |>
    dplyr::full_join(site_covs) |>
    dplyr::rename( xvar = x) |> 
    ( function(x) dplyr::mutate(x,
                                rho = rgamma(nrow(x), zeta, zeta),
                                eng = exp( alpha0 + alpha1 * xvar) * rho,
                                ng = rpois(nrow(x), eng)))()
  
  
  ng_vector <- c()
  site_vector <- c()
  rep_vector <- c()
  sp_vector <- c()
  for(i in 1 : nrow( ng_df )) {
    if( ng_df[[i, "ng"]] == 0){
      ng_vector <- c(ng_vector, 0)
      site_vector <- c(site_vector, ng_df[[i, "site"]])
      rep_vector <- c(rep_vector, ng_df[[i, "rep"]])
      sp_vector <- c(sp_vector, ng_df[[i, "sp"]])
    } else {
      ng_vector <- c(ng_vector, rep(1, ng_df[[i, "ng"]]))
      site_vector <- c(site_vector, rep(ng_df[[i, "site"]], ng_df[[i, "ng"]]))
      rep_vector <- c(rep_vector, rep(ng_df[[i, "rep"]], ng_df[[i, "ng"]]))
      sp_vector <- c(sp_vector, rep(ng_df[[i, "sp"]], ng_df[[i, "ng"]]))
    }
  }
  
  # expanded df so each group can have a dclass :)
  ng_df_expanded <- tibble::tibble(
    site = site_vector, 
    rep = rep_vector, 
    sp = sp_vector, 
    group = ng_vector) |> # group is just a placeholder - means yes, there is a group
    dplyr::full_join(ng_df)
  
  # assign distances to each group and simulate observation process, based on distance
  sigma <- exp(sp_df$gamma0)
  data <- NULL
  for( i in 1 : nrow(ng_df_expanded) ) {
    if(ng_df_expanded[[i, "ng"]] == 0){
      data <- tibble::as_tibble(
        rbind(data,
              cbind(
                site = ng_df_expanded[[i, "site"]],
                rep = ng_df_expanded[[i, "rep"]],
                sp = ng_df_expanded[[i, "sp"]],
                group = ng_df_expanded[[i, "group"]],
                eng = ng_df_expanded[[i, "eng"]],
                ng = ng_df_expanded[[i, "ng"]],
                group_obs = 0, 
                dclass = NA)))
    } else {
      d <- runif( 1, 0, b) # animals distributed uniformly
      dclass <- d %/% width + 1 # grab the dclass that it falls into
      # detection probability is a function of distance and the scale parameter
      p <- exp( -d * d / (2 * sigma[ng_df_expanded[[i, "sp"]]] ^ 2))
      # was or was not the group observed?
      group_obs <- rbinom(ng_df_expanded[[i, "group"]], 1, p)
      
      data <- tibble::as_tibble(
        rbind(data,
              cbind(
                site = ng_df_expanded[[i, "site"]],
                rep = ng_df_expanded[[i, "rep"]],
                sp = ng_df_expanded[[i, "sp"]],
                group = ng_df_expanded[[i, "group"]],
                eng = ng_df_expanded[[i, "eng"]],
                ng = ng_df_expanded[[i, "ng"]],
                group_obs = group_obs, 
                dclass = dclass)))
    }
  }
  
  ds_data_all <- sp_df |> 
    dplyr::full_join(data) |>
    ( function(x) dplyr::mutate(x,
                                phi = rgamma(nrow(x), xi, xi),
                                egs = exp(beta0) * phi))() |>
    ( function(x) dplyr::mutate(x,
                                gs = ifelse(group == 0, 0, # if there is not a group, observed gs must be 0
                                            extraDistr::rtpois( nrow(x),
                                                                lambda = egs,
                                                                a = 0,
                                                                b = Inf))))() |> 
    dplyr::mutate(gs_obs = ifelse(group_obs == 1, gs, 0)) |> 
    dplyr::mutate(dclass = ifelse(gs_obs == 0, NA, dclass)) |> 
    dplyr::ungroup() |> 
    dplyr::group_by(sp, site, rep) |> 
    dplyr::mutate(ng_obs = sum(group_obs)) |> 
    dplyr::select(sp, site, rep, 
                  dclass,
                  eng, ng, ng_obs,
                  egs, gs, gs_obs) |> 
    dplyr::distinct()
  
  ds_data <- ds_data_all |>  
    dplyr::filter(! (ng_obs > 0 & is.na(dclass))) |> 
    dplyr::ungroup() |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::mutate(rowid = row_number()) |> 
    dplyr::group_by(sp) |> 
    dplyr::mutate(nobs_start = min(rowid), 
                  nobs_end = max(rowid)) 
  
  g <- function(x, sig) exp(- x * x /(2 * sig^2))
  dist_breaks <- seq(0, b, by = width)
  p <- array(NA, 
             dim = c(length(dist_breaks) - 1, nsp))
  
  for(j in 1:nrow(p)) {
    for(s in 1:nsp) {
      p[j, s] <- ( integrate( g, 
                              dist_breaks[j], 
                              dist_breaks[j + 1], 
                              sig = sigma[s])$value / (dist_breaks[j + 1] - dist_breaks[j])) * (width / b)
    }
  }
  
  p_tc <- c()
  for(i in 1:ncol(p)) {
    p_tc <- c(p_tc, sum(p[, i]))
  }
  
  if(p_bias == -1){
    p_tc <- p_tc - (p_tc * 0.20)
  } else if(p_bias == 1){
    p_tc <- p_tc + (p_tc * 0.20)
  }
  
  sp_df_tc <-
    cbind(sp_df, p_tc) |> 
    tibble::as_tibble()
  
  tc_site_covs <- tibble::tibble(
    site = 1:nsites_tc,
    x_tc = runif(nsites_tc, -2, 2)) |> 
    dplyr::mutate(x_tc = as.numeric(scale(x_tc)))
  
  transect_counts <- expand.grid(sp = 1:nsp,
                                 site =  1:nsites_tc,
                                 rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df_tc) |>
    dplyr::full_join(tc_site_covs) |> 
    ( function(x) dplyr::mutate(x,
                                rho = rgamma(nrow(x), zeta, zeta),
                                phi = rgamma(nrow(x), xi, xi),
                                eng = exp( alpha0 + alpha1 * x_tc) * rho,
                                egs = exp( beta0 ) * phi,
                                lambda_ltn = eng * egs,
                                totalN = rpois(nrow(x), lambda_ltn),
                                count = rbinom(nrow(x), totalN, p_tc)))() |> 
    dplyr::arrange(sp) |> 
    dplyr::mutate(rowid = row_number()) |> 
    dplyr::group_by(sp) |> 
    dplyr::mutate(ntc_start = min(rowid), 
                  ntc_end = max(rowid))
  
  ng_data <- ds_data_all |>
    dplyr::group_by(sp, site, rep) |> 
    dplyr::summarise( ng = unique(ng), 
                      ng_obs = unique(ng_obs), 
                      N = sum(gs), 
                      yN_DS = sum(gs_obs) ) |>
    dplyr::ungroup() |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::mutate(rowid = row_number()) |> 
    dplyr::group_by(sp) |> 
    dplyr::mutate(ngobs_start = min(rowid), 
                  ngobs_end = max(rowid)) |> 
    dplyr::full_join(site_covs)
  
  ng_indices <- ng_data |> 
    dplyr::select(sp, ngobs_start, ngobs_end) |> 
    dplyr::distinct()
  
  transect_indices <- transect_counts |> 
    dplyr::select(sp, ntc_start, ntc_end) |> 
    dplyr::distinct()
  
  ds_data_final <- ds_data |> 
    dplyr::select(sp, site, rep, dclass, gs_obs) |> 
    dplyr::filter(!is.na(dclass)) |> 
    dplyr::arrange(sp, site, rep)
  
  focal_sp <- ds_data_final |>
    dplyr::count(sp) |>
    dplyr::ungroup() |>
    dplyr::full_join(tibble::tibble(sp = 1:20)) |>
    dplyr::mutate(n = replace_na(n, 0)) |>
    dplyr::arrange(n) |>
    dplyr::filter(n > 0) |>
    dplyr::filter(n < max(n)) |>
    dplyr::slice(c(1, n())) |>
    tibble::add_column(stat = c("rare", "common"))
  
  data <- list(
    MIDPOINT = seq(from = 12.5, to = 987.5, by = 25),
    yGS = ds_data_final$gs_obs, 
    DCLASS = ds_data_final$dclass, 
    V = 25, 
    B = 1000,
    yNG = ng_data$ng_obs,
    # yN_DS = ng_data$yN_DS,
    true_N_DS = ng_data$N,
    HAB_DS = ng_data$x,
    HAB_TC = transect_counts$x_tc,
    yN = transect_counts$count,
    true_N = transect_counts$totalN,
    true_alpha0 = alpha0, 
    true_alpha1 = alpha1, 
    true_beta0 = beta0)
  
  constants <- list(
    NSPECIES = length(unique(transect_counts$sp)),
    NBINS = length(data$MIDPOINT), 
    NDISTANCES = nrow(ds_data_final), 
    SP_GS = ds_data_final$sp, 
    SP_NG = ng_data$sp,
    NSURVEYS = nrow(ng_data),  
    NCOUNTS = nrow(transect_counts), 
    SP_TC = transect_counts$sp,
    COUNT_START = transect_indices$ntc_start,
    COUNT_END = transect_indices$ntc_end,
    NG_START = ng_indices$ngobs_start,
    NG_END = ng_indices$ngobs_end)
  
  sp_info <- ng_data |> 
    dplyr::group_by(sp) |> 
    dplyr::summarise(totDS = sum(N)) |> 
    dplyr::full_join(dplyr::summarise(dplyr::group_by(transect_counts, sp), totTC = sum(totalN))) |> 
    dplyr::full_join(sp_df) |> 
    dplyr::full_join(focal_sp)
  
  return(list(data = data,
              constants = constants,
              sp_info = sp_info))
}

code_cds <- nimble::nimbleCode({
  
  mu_gamma0    ~ dnorm(5.5, sd = 0.5)
  sigma_gamma0 ~ dexp(1)
  mu_alpha0    ~ dnorm(0, sd = 2) 
  sigma_alpha0 ~ dexp(1)          
  mu_alpha1    ~ dnorm(0, sd = 2)
  sigma_alpha1 ~ dexp(1)
  mu_beta0     ~ dnorm(0, sd = 2)  
  sigma_beta0  ~ dexp(1)        
  
  for ( s in 1:NSPECIES ) {
    gamma0[s] ~ dnorm(mu_gamma0, sd = sigma_gamma0)
    omega[s] <- exp(gamma0[s])
    pie_sp[s] <- sum( pie[1:NBINS, s] )
    N_RB_DS[s] <- sum ( N_diff_ds[ NG_START[s]:NG_END[s]]) / sum(true_N_DS[NG_START[s]:NG_END[s]])
    alpha0[s] ~ dnorm(mu_alpha0, sd = sigma_alpha0)
    alpha1[s] ~ dnorm(mu_alpha1, sd = sigma_alpha1)
    alpha1_diff[s] <- true_alpha1[s] - alpha1[s]
    beta0[s]  ~ dnorm(mu_beta0,  sd = sigma_beta0) 
    zeta[s] ~ dgamma(4, 2)    
    xi[s] ~ dgamma(4, 2) 
    
    for (k in 1:NBINS ) {
      log(g[k, s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega[s] * omega[s])
      pie[k, s] <- g[k,s] * (V/B)
      pie_cell[k, s] <- pie[k, s] / pie_sp[s]
    }
  }
  
  for( i in 1:NSURVEYS ) {
    log( lambda[i] ) <- alpha0[ SP_NG[i]] + alpha1[ SP_NG[i]] * HAB_DS[i] 
    rho[i] ~ dgamma( zeta[SP_NG[i]], zeta[SP_NG[i]] )
    NG[i] ~ dpois( lambda[i] * rho[i] )
    yNG[i] ~ dbin( pie_sp[SP_NG[i]], NG[i] )
    phi_ng[i] ~ dgamma( xi[SP_NG[i]], xi[SP_NG[i]] )
    log(kappa_ng[i]) <- beta0[SP_NG[i]]
    GS[i] ~ T(dpois( kappa_ng[i] * phi_ng[i]), 1, )
    N_ds[i] <- NG[i] * GS[i]
    N_diff_ds[i] <- N_ds[i] - true_N_DS[i]
  }
  
  for (i in 1:NDISTANCES ) {
    log( kappa[i] ) <- beta0[SP_GS[i]]
    phi[i] ~ dgamma( xi[SP_GS[i]], xi[SP_GS[i]])
    yGS[i] ~ T(dpois( kappa[i] * phi[i] ), 1, )
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i] ] )
  }
  
})

params_cds <- c(
  "alpha1_diff", 
  "N_RB_DS")

inits_cds <- function(){ 
  list(
    NG = data_run$yNG + 1,
    mu_gamma0 = rnorm(1, 5.5, 0.5), 
    sigma_gamma0 = rexp(1, 1),
    mu_alpha0 = runif(1, 1.5, 2.5),
    sigma_alpha0 = rexp(1, 1),
    mu_alpha1 = runif(1, 0, 1), 
    sigma_alpha1 = runif(1, 0.5, 1),
    mu_beta0 = runif(1, 2, 3),
    sigma_beta0 = rexp(1, 1), 
    gamma0 = rnorm(constants_run$NSPECIES, 5.5, 0.5),
    alpha0 = runif(constants_run$NSPECIES, 1, 2), 
    alpha1 = rnorm(constants_run$NSPECIES, 0, 1),
    beta0 = runif(constants_run$NSPECIES, 1, 2),
    zeta = runif(constants_run$NSPECIES, 0.2, 0.8),
    xi = runif(constants_run$NSPECIES, 0.2, 0.8),
    rho = rgamma(constants_run$NSURVEYS, 2, 2),
    phi = rgamma(constants_run$NDISTANCES, 2, 2),
    phi_ng = rgamma(constants_run$NSURVEYS, 2, 2))
}

nc <- 3
nburn <- 50000
ni <- nburn + 100000
nt <- 100

nsimreps <- 100

outs <- list(list())
alpha1_diff <- list(list())
N_RB_DS <- list(list())
sp_info <- list(list())
simdat <- list(list())

for( i in 1: nsimreps ) {
  
  simdat[[i]] <- sim_icm()
  data_run <- simdat[[i]]$data
  constants_run <- simdat[[i]]$constants
  model.code <- code_cds
  inits <- inits_cds
  params <- params_cds
  
  start <- Sys.time()
  cl <- parallel::makeCluster(nc)
  parallel::clusterExport(cl, c("model.code",
                                "inits",
                                "data_run",
                                "constants_run",
                                "params",
                                "nburn",
                                "ni",
                                "nt"))
  
  for(j in seq_along(cl)) {
    set.seed(j)
    init <- inits()
    set.seed(NULL)
    parallel::clusterExport(cl[j], "init")
  }
  
  out <- parallel::clusterEvalQ(cl, {
    library(nimble)
    library(coda)
    
    model <- nimble::nimbleModel(code = model.code,
                                 name = "model.code",
                                 constants = constants_run,
                                 data = data_run,
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
    
    return(coda::as.mcmc(out1))
  })
  end <- Sys.time()
  stopCluster(cl)
  
  outs[[i]] <- out
  
  sp_info[[i]] <- simdat[[i]]$sp_info |>  
    add_column(simrep = i)
  
  alpha1_diff[[i]] <- MCMCvis::MCMCpstr(out, params = "alpha1_diff", type = "chains" )[[1]] |> 
    tibble::as_tibble(rownames = "sp") |> 
    ( function(x) tidyr::pivot_longer(x,
                                      2:ncol(x),
                                      names_to = "iter",
                                      values_to = "value"))() |> 
    dplyr::mutate(sp = stringr::str_replace(sp, "alpha1_diff", "")) |> 
    dplyr::mutate(sp = readr::parse_number(sp), 
                  iter = readr::parse_number(iter)) |> 
    tibble::add_column(param = "alpha1_diff", 
                       simrep = i) |> 
    dplyr::select(simrep, sp, iter, param, value) 
  
  N_RB_DS[[i]] <- MCMCvis::MCMCpstr(out, params = "N_RB_DS", type = "chains" )[[1]] |> 
    tibble::as_tibble(rownames = "sp") |> 
    ( function(x) tidyr::pivot_longer(x,
                                      2:ncol(x),
                                      names_to = "iter",
                                      values_to = "value"))() |> 
    dplyr::mutate(sp = stringr::str_replace(sp, "N_RB_DS", "")) |> 
    dplyr::mutate(sp = readr::parse_number(sp), 
                  iter = readr::parse_number(iter)) |> 
    tibble::add_column(param = "N_RB_DS", 
                       simrep = i) |> 
    dplyr::select(simrep, sp, iter, param, value) 
  
  rm(data_run, constants_run)
  print(paste("Finished", i, "of", nsimreps ))
  
}

setwd(here::here("results"))
save( outs, sp_info, alpha1_diff, N_RB_DS, file = paste0("simulation_comparison_cds_", Sys.Date(), ".RData"))