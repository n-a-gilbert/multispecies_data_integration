# simulation script for single-species, distance-sampling-only model for rare species
library(tidyverse)
library(nimble)
library(parallel)
library(MCMCvis)
library(here)

# see main simulation script for more comments / details
sim_icm <- function(
    nsp = 15,          # number of species
    mu_alpha0 = 0.87,
    sigma_alpha0 = 1.95,
    mu_alpha1 = 0.05,    # community average for covariate effect
    sigma_alpha1 = 0.25,   # standard deviation among species for covariate effect
    mu_gamma0_ds = 5.5,
    mu_gamma0_c = 5.0,
    sigma_gamma0_ds = 0.25,
    sigma_gamma0_c = 0.25,
    nsites = 50,       # number of sites for distance sampling
    nrep = 1,          # number of temporal replicates
    b = 1000,          # distance to which animals are counted
    width = 25,        # width of distance classes
    nsites_tc_fact = 2 # multiplication factor of how much more count data sites there are
){
  
  nsites_tc <- nsites * nsites_tc_fact
  
  # abundance - intercept & covariate coefficient
  alpha0 <- rnorm( nsp, mean = mu_alpha0, sd = sigma_alpha0 )
  alpha1 <- rnorm( nsp, mean = mu_alpha1, sd = sigma_alpha1 )
  
  # intercept for scale parameter
  gamma0_ds <- rnorm( nsp, mean = mu_gamma0_ds, sd = sigma_gamma0_ds )

  sp_df <- tibble::tibble(
    sp = 1:nsp,
    alpha0 = alpha0,
    alpha1 = alpha1, 
    gamma0_ds = gamma0_ds)
  
  com_truth <- tibble::tribble(
    ~param, ~truth,
    "mu_gamma0", mu_gamma0_ds,
    "sd_gamma0", sigma_gamma0_ds,
    "mu_alpha0", mu_alpha0,
    "sd_alpha0", sigma_alpha0,
    "mu_alpha1", mu_alpha1,
    "sd_alpha1", sigma_alpha1) 
  
  site_covs <- tibble::tibble(
    site = 1:nsites,
    x = runif(nsites, -2, 2)) |> 
    dplyr::mutate(x = as.numeric(scale(x)))
  
  n_df <- expand.grid(sp = 1:nsp,
                      site =  1:nsites,
                      rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df) |>
    dplyr::full_join(site_covs) |>
    dplyr::rename( xvar = x) |> 
    ( function(x) dplyr::mutate(x,
                                en = exp( alpha0 + alpha1 * xvar),
                                n = rpois(nrow(x), en)))()  |> 
    dplyr::rowwise()  |>  
    # how many groups were there? (For assigning distance measurements)
    dplyr::mutate( ng = ifelse(n > 0, sample(1:n, 1), 0))  |> 
    dplyr::ungroup()
  
  n_vector <- c()
  site_vector <- c()
  rep_vector <- c()
  sp_vector <- c()
  for(i in 1 : nrow( n_df )) {
    if( n_df[[i, "n"]] == 0){
      n_vector <- c(n_vector, 0)
      site_vector <- c(site_vector, n_df[[i, "site"]])
      rep_vector <- c(rep_vector, n_df[[i, "rep"]])
      sp_vector <- c(sp_vector, n_df[[i, "sp"]])
    } else {
      n_vector <- c(n_vector, rep(1, n_df[[i, "ng"]]))
      site_vector <- c(site_vector, rep(n_df[[i, "site"]], n_df[[i, "ng"]]))
      rep_vector <- c(rep_vector, rep(n_df[[i, "rep"]], n_df[[i, "ng"]]))
      sp_vector <- c(sp_vector, rep(n_df[[i, "sp"]], n_df[[i, "ng"]]))
    }
  }
  
  get_unique_integers <- function(n, ng){
    mat <- rmultinom(n, size = 1, prob = c(runif(ng, 0, 0.5)))
    rows <- apply(mat, 1, sum)
    return( rows )
  }
  
  # expanded df so each group can have a dclass :)
  n_df_expanded <- tibble::tibble(
    site = site_vector, 
    rep = rep_vector, 
    sp = sp_vector, 
    group = n_vector) |> # group is just a placeholder - means yes, there is a group
    dplyr::full_join(n_df) |> 
    dplyr::group_by(sp, site, rep) |> 
    dplyr::mutate(gs = ifelse(ng == 0, 0, 
                       get_unique_integers(n = n, ng = ng)))  |> 
    dplyr::ungroup()
  
  
  # assign distances to each group and simulate observation process, based on distance
  sigma <- exp(sp_df$gamma0_ds)
  data <- NULL
  for( i in 1 : nrow(n_df_expanded) ) {
    if(n_df_expanded[[i, "ng"]] == 0){
      data <- tibble::as_tibble(
        rbind(data,
              cbind(
                site = n_df_expanded[[i, "site"]],
                rep = n_df_expanded[[i, "rep"]],
                sp = n_df_expanded[[i, "sp"]],
                group = n_df_expanded[[i, "group"]],
                eng = n_df_expanded[[i, "eng"]],
                n = n_df_expanded[[i, "n"]],
                ng = n_df_expanded[[i, "ng"]],
                gs = n_df_expanded[[i, "gs"]],
                group_obs = 0, 
                dclass = NA)))
    } else {
      d <- runif( 1, 0, b) # animals distributed uniformly
      dclass <- d %/% width + 1 # grab the dclass that it falls into
      # detection probability is a function of distance and the scale parameter
      p <- exp( -d * d / (2 * sigma[n_df_expanded[[i, "sp"]]] ^ 2))
      # was or was not the group observed?
      group_obs <- rbinom(n_df_expanded[[i, "group"]], 1, p)
      
      data <- tibble::as_tibble(
        rbind(data,
              cbind(
                site = n_df_expanded[[i, "site"]],
                rep = n_df_expanded[[i, "rep"]],
                sp = n_df_expanded[[i, "sp"]],
                group = n_df_expanded[[i, "group"]],
                eng = n_df_expanded[[i, "eng"]],
                n = n_df_expanded[[i, "n"]],
                ng = n_df_expanded[[i, "ng"]],
                gs = n_df_expanded[[i, "gs"]],
                group_obs = group_obs, 
                dclass = dclass)))
    }
  }
  
  # filter out only the most rare species
  # defined as the species with at least one distance observation with the smallest sum of observed counts
  ng_data <- data |> 
    dplyr::filter(gs > 0) |> 
    dplyr::filter(group_obs == 1) |> 
    dplyr::group_by(sp, site, rep) |> 
    dplyr::summarise( true_n = unique(n), 
                      count = sum(gs),
                      ng = sum(gs > 0)) |> 
    dplyr::full_join(
      dplyr::select(n_df, sp, site, rep, true_n = n)
    ) |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::mutate(count = tidyr::replace_na(count, 0),
                  ng = tidyr::replace_na(ng, 0)) |> 
    dplyr::full_join(site_covs) |>
    dplyr::group_by(sp) |>
    dplyr::mutate(totDS = sum(true_n), 
                  totDS_obs = sum(count),
                  ndistances = sum(ng)) |>
    dplyr::ungroup() |>
    dplyr::filter( totDS_obs > 1 & ndistances > 1 ) |>
    dplyr::filter(totDS_obs == min(totDS_obs)) |>
    dplyr::arrange( sp, site, rep) |>
    dplyr::slice(1:nsites)
  
  ds_data_final <- data |> 
    dplyr::filter(gs > 0) |> 
    dplyr::filter(group_obs == 1) |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::select(sp, site, rep, gs, dclass) |>
    dplyr::filter( sp == unique(ng_data$sp))
  
  data <- list(
    MIDPOINT = seq(from = 12.5, to = 987.5, by = 25),
    DCLASS = ds_data_final$dclass, 
    V = 25, 
    B = 1000,
    yN_DS = ng_data$count,
    HAB_DS = ng_data$x,
    true_n_ds = ng_data$true_n)
  
  constants <- list(
    NSPECIES = length(unique(ng_data$sp)),
    NBINS = length(data$MIDPOINT), 
    NDISTANCES = nrow(ds_data_final), 
    SP_GS = ds_data_final$sp - (unique(ds_data_final$sp) - 1),
    SP_NG = ng_data$sp - (unique(ng_data$sp) - 1),
    NSURVEYS = nrow(ng_data))
  
  sp_info <- ng_data |>
    dplyr::group_by(sp) |>
    dplyr::summarise(totDS = unique(totDS)) |>
    dplyr::left_join(sp_df)
  
  return(list(data = data,
              constants = constants,
              sp_info = sp_info,
              com_truth = com_truth))
}

# see main simulation for more comments
# here, community-level parameters and count submodel are omitted
model.code <- nimble::nimbleCode({
  
  for(s in 1:NSPECIES){
    gamma0_ds[s] ~ dunif(0, 10)
    alpha0[s] ~ dnorm(0, sd = 2)
    alpha1[s] ~ dnorm(0, sd = 2)
    omega_ds[s] <- exp(gamma0_ds[s])
    pie_sp[s] <- sum( pie[1:NBINS, s])
    for (k in 1:NBINS ) {
      log(g[k,s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega_ds[s] * omega_ds[s] )
      pie[k,s] <- g[k,s] * (V/B)
      pie_cell[k,s] <- pie[k,s] / pie_sp[s]
    }
  }
  for( i in 1:NSURVEYS ) {
    log( lambda[i] ) <- alpha0[ SP_NG[i]] + alpha1[ SP_NG[i]] * HAB_DS[i] 
    N_DS[i] ~ dpois( lambda[i] )
    yN_DS[i] ~ dbin( pie_sp[SP_NG[i]], N_DS[i] )
  }
  
  for (i in 1:NDISTANCES ) {
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i] ] )
  }
})

params <- c(
  "gamma0_ds", 
  "alpha0",
  "alpha1",
  "pie_sp",
  "N_DS")

make_inits <- function(data, constants) { 
  inits <- list(
    N_DS = data$yN_DS + 1,
    gamma0_ds = rnorm(1, 5.5, 0.5),
    alpha0 = rnorm(1, 0, 1), 
    alpha1 = rnorm(1, 0, 1))
  return(inits)
}

nburn <- 100000
ni <- nburn + 100000
nt <- 100
nc <- 3

min_simrep <- 1
max_simrep <- 1000

simrep_rank <- rank(min_simrep:max_simrep)
simrep_raw <- min_simrep:max_simrep

for( i in min(simrep_rank):max(simrep_rank)){
  
  simdat <- sim_icm()
  data <- simdat$data
  constants <- simdat$constants
  sp_info <- simdat$sp_info
  com_truth <- simdat$com_truth
  print(paste( "Starting rep", simrep_rank[i], "of", max(simrep_rank))) 
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
    
    model <- nimble::nimbleModel(code = model.code,
                         name = "model.code",
                         constants = constants,
                         data = data,
                         inits = init)
    
    Cmodel <- nimble::compileNimble(model)
    modelConf <- nimble::configureMCMC(model)
    modelConf$addMonitors(params)
    modelMCMC <- nimble::buildMCMC(modelConf)
    CmodelMCMC <- nimble::compileNimble(modelMCMC, project = model)
    out1 <- nimble::runMCMC(CmodelMCMC, 
                    nburnin = nburn, 
                    niter = ni, 
                    thin = nt)
    
    return(coda::as.mcmc(out1))
  })
  end <- Sys.time()
  time <- difftime(end, start, units = "hours")
  parallel::stopCluster(cl)
  
  outsum <- MCMCvis::MCMCsummary( out ) |> 
    tibble::as_tibble(rownames = "param")
  
  res <- sp_info |> 
    tidyr::pivot_longer(c("gamma0_ds", 
                   "alpha0",
                   "alpha1"), 
                 names_to = "param", values_to = "truth") |>
    dplyr::mutate(sp = sp - ( min(sp) - 1)) |> 
    dplyr::mutate(param = paste0(param, '[', sp, ']')) |> 
    dplyr::select(param, totDS, truth) |> 
    dplyr::full_join(
      dplyr::full_join( dplyr::mutate( dplyr::select( sp_info, sp, totDS), sp = sp - (min(sp) - 1)),
                 tibble::tibble(
                   sp = constants$SP_NG,
                   param = paste0("N_DS[", 1:length(data$true_n_ds), "]"),
                   truth = data$true_n_ds)
      )
    ) |> 
    dplyr::left_join(outsum) |> 
    tibble::add_column(simrep = simrep_raw[i])
  
  readr::write_csv(res, paste0("ssds_rare_no_od_simrep_", formatC(simrep_raw[i], width = 4, format = "d", flag = "0"), "_results.csv"))
  print(paste("Rep", simrep_rank[i], "took", round(time[[1]], 3), "hours"))
  rm( cl, com_truth, constants, data, init, out, outsum, res, simdat, sp_info, end, start, time)
}
