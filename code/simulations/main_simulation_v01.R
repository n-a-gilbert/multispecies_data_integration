library(tidyverse)
library(nimble)
library(parallel)
library(MCMCvis)

# functtion to simulate distance sampling and count data for a community of species
sim_icm <- function(
    nsp = 15,               # number of species
    mu_alpha0 = 0.87,       # community average for abundance intercept
    sigma_alpha0 = 1.95,    # variation in abundance intercepts among species
    mu_alpha1 = 0.05,       # community average for covariate effect
    sigma_alpha1 = 0.25,    # standard deviation among species for covariate effect
    mu_gamma0_ds = 5.5,     # community average of intercept for scale parameter in DS detection function
    mu_gamma0_c = 5.0,      # community average of intercept for scale parameter in count detection function
    sigma_gamma0_ds = 0.25, # standard deviation among species for intercept for scale parameter in DS detection function
    sigma_gamma0_c = 0.25,  # standard deviation among species for intercept for scale parameter in count detection function
    nsites = 50,            # number of sites for distance sampling
    nrep = 1,               # number of temporal replicates
    b = 1000,               # maximum distance to which animals are counted
    width = 25,             # width of distance classes
    nsites_tc_fact = 2      # multiplication factor of how much more count data sites there are
){
  
  # how many sites with count data 
  nsites_tc <- nsites * nsites_tc_fact
  
  # simulate species-level abundance intercept
  alpha0 <- rnorm( nsp, mean = mu_alpha0, sd = sigma_alpha0 )
  
  # simulate species-level covariate effect
  alpha1 <- rnorm( nsp, mean = mu_alpha1, sd = sigma_alpha1 )
  
  # simulate species-level intercept for scale parameter for DS detection function
  gamma0_ds <- rnorm( nsp, mean = mu_gamma0_ds, sd = sigma_gamma0_ds )
  
  # simulate species-level intercept for scale parameter for count detection function
  gamma0_c <- rnorm( nsp, mean = mu_gamma0_c, sd = sigma_gamma0_c )
  
  # put the true values for the species-level parameters in a dataframe
  sp_df <- tibble::tibble(
    sp = 1:nsp,
    alpha0 = alpha0,
    alpha1 = alpha1, 
    gamma0_ds = gamma0_ds, 
    gamma0_c = gamma0_c)
  
  # put the true values for the community-level parameters in a dataframe
  com_truth <- tibble::tribble(
    ~param, ~truth,
    "mu_gamma0", mu_gamma0_ds,
    "sd_gamma0", sigma_gamma0_ds,
    "mu_gamma0_c", mu_gamma0_c,
    "sd_gamma0_c", sigma_gamma0_c,
    "mu_alpha0", mu_alpha0,
    "sd_alpha0", sigma_alpha0,
    "mu_alpha1", mu_alpha1,
    "sd_alpha1", sigma_alpha1) 
  
  # simulate a coveriate associated with abundance for distance sampling sites
  site_covs <- tibble::tibble(
    site = 1:nsites,
    x = runif(nsites, -2, 2)) |> 
    mutate(x = as.numeric(scale(x)))
  
  # simulate abundance for distance sampling sites
  n_df <- expand.grid(sp = 1:nsp,
                      site =  1:nsites,
                      rep =  1:nrep) |> 
    tibble::as_tibble() |> 
    dplyr::full_join(sp_df) |>
    dplyr::full_join(site_covs) |>
    dplyr::rename( xvar = x) |>  # covariate
    ( function(x) dplyr::mutate(x,
                                en = exp( alpha0 + alpha1 * xvar), # expected abundance - function of intercept, slope, and covariate
                                n = rpois(nrow(x), en)))()  |>    # latent abundance - simulate with Poisson
    dplyr::rowwise()  |>  
    # how many groups were there? (For assigning distance measurements)
    dplyr::mutate( ng = ifelse(n > 0, sample(1:n, 1), 0))  |> 
    dplyr::ungroup()
  
  # expand the abundance dataframe so each group of animals has its own rown
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
    group_by(sp, site, rep) |> 
    mutate(gs = ifelse(ng == 0, 0, 
                       get_unique_integers(n = n, ng = ng))) |> # this gives us group sizes that sums up to the true abundance 
    ungroup()
  
  
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
      p <- exp( -d * d / (2 * sigma[n_df_expanded[[i, "sp"]]] ^ 2)) # half-normal detection function for detection probability
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
  
  # Now, we do the same thing, but for the count data
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
                                en = exp( alpha0 + alpha1 * xvar),
                                n = rpois(nrow(x), en)))() |> 
    dplyr::rowwise()  |>  
    # how many groups were there? (For assigning distance measurements)
    dplyr::mutate( ng = ifelse(n > 0, sample(1:n, 1), 0))  |> 
    dplyr::ungroup()
  
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
      n_vector_c <- c(n_vector_c, rep(1, n_df_c[[i, "ng"]]))
      site_vector_c <- c(site_vector_c, rep(n_df_c[[i, "site"]], n_df_c[[i, "ng"]]))
      rep_vector_c <- c(rep_vector_c, rep(n_df_c[[i, "rep"]], n_df_c[[i, "ng"]]))
      sp_vector_c <- c(sp_vector_c, rep(n_df_c[[i, "sp"]], n_df_c[[i, "ng"]]))
    }
  }
  
  # expanded df so each group can have a dclass :)
  n_df_expanded_c <- tibble::tibble(
    site = site_vector_c, 
    rep = rep_vector_c, 
    sp = sp_vector_c, 
    group = n_vector_c) |> # group is just a placeholder - means yes, there is a group
    dplyr::full_join(n_df_c) |> 
    group_by(sp, site, rep) |> 
    mutate(gs = ifelse(ng == 0, 0, 
                       get_unique_integers(n = n, ng = ng)))  |> 
    ungroup()
  
  # assign distances to each group and simulate observation process, based on distance
  sigmaC <- exp(sp_df$gamma0_c)
  data_c <- NULL
  for( i in 1 : nrow(n_df_expanded_c) ) {
    if(n_df_expanded_c[[i, "ng"]] == 0){
      data_c <- tibble::as_tibble(
        rbind(data_c,
              cbind(
                site = n_df_expanded_c[[i, "site"]],
                rep = n_df_expanded_c[[i, "rep"]],
                sp = n_df_expanded_c[[i, "sp"]],
                group = n_df_expanded_c[[i, "group"]],
                eng = n_df_expanded_c[[i, "eng"]],
                n = n_df_expanded_c[[i, "n"]],
                ng = n_df_expanded_c[[i, "ng"]],
                gs = n_df_expanded_c[[i, "gs"]],
                group_obs = 0, 
                dclass = NA)))
    } else {
      d <- runif( 1, 0, b) # animals distributed uniformly
      dclass <- d %/% width + 1 # grab the dclass that it falls into
      # detection probability is a function of distance and the scale parameter
      p <- exp( -d * d / (2 * sigmaC[n_df_expanded_c[[i, "sp"]]] ^ 2))
      # was or was not the group observed?
      group_obs <- rbinom(n_df_expanded_c[[i, "group"]], 1, p)
      
      data_c <- tibble::as_tibble(
        rbind(data_c,
              cbind(
                site = n_df_expanded_c[[i, "site"]],
                rep = n_df_expanded_c[[i, "rep"]],
                sp = n_df_expanded_c[[i, "sp"]],
                group = n_df_expanded_c[[i, "group"]],
                eng = n_df_expanded_c[[i, "eng"]],
                n = n_df_expanded_c[[i, "n"]],
                ng = n_df_expanded_c[[i, "ng"]],
                gs = n_df_expanded_c[[i, "gs"]],
                group_obs = group_obs, 
                dclass = dclass)))
    }
  }
  
  # this is the final product for the count data
  # we discard the distance observations and just keep the number of animals detected
  transect_counts <- data_c |> 
    dplyr::filter( gs > 0) |> 
    dplyr::filter(group_obs == 1) |> 
    dplyr::group_by(sp, site, rep) |> 
    summarise( count = sum(gs)) |> 
    ungroup() |> 
    full_join(
      dplyr::select( n_df_c, sp, site, rep, true_n = n)
    ) |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::mutate(count = tidyr::replace_na(count, 0)) |> 
    dplyr::full_join(site_covs_c) |> 
    dplyr::select(sp, site, rep, true_n, count, x_tc = x)
  
  # distance sampling data - generate counts for each survey
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
    dplyr::mutate(totDS = sum(true_n),      # total abundance (across sites)
                  totDS_obs = sum(count),    # total count (across sites)
                  ndistances = sum(ng)) |>  # number of distance observations
    dplyr::ungroup()
  
  # the distance observations
  ds_data_final <- data |> 
    dplyr::filter(gs > 0) |> 
    dplyr::filter(group_obs == 1) |> 
    dplyr::arrange(sp, site, rep) |> 
    dplyr::select(sp, site, rep, gs, dclass)
  
  # package up data for NIMBLE
  data <- list(
    MIDPOINT = seq(from = 12.5, to = 987.5, by = 25), # midpoint of each distance band
    DCLASS = ds_data_final$dclass,                    # distance class observations
    V = 25,                                           # distance band width
    B = 1000,                                         # maximum distance to which animals are counted
    yN_DS = ng_data$count,                            # distance sampling count
    HAB_DS = ng_data$x,                               # distance sampling abundance coveriate
    HAB_TC = transect_counts$x_tc,                    # count abundance covariate
    yN_TC = transect_counts$count,                    # count count
    true_n_ds = ng_data$true_n,                       # latent abundance at distance sampling sites
    true_n_tc = transect_counts$true_n)               # latent abundance at count sites
  
  # constants (control loops, etc) for NIMBLE model
  constants <- list(
    NSPECIES = length(unique(transect_counts$sp)), # number of species
    NBINS = length(data$MIDPOINT),                 # number of distance bands
    NDISTANCES = nrow(ds_data_final),              # number of distance observations
    SP_GS = ds_data_final$sp,                      # species index for the distance observations 
    SP_NG = ng_data$sp,                            # species index for distance sampling counts
    NSURVEYS = nrow(ng_data),                      # number of distance sampling surveys
    NCOUNTS = nrow(transect_counts),               # number of count surveys
    SP_TC = transect_counts$sp)                    # species index for count data
  
  # summarise species-level information
  sp_info <- ng_data |>
    dplyr::select(sp, totDS, totDS_obs, ndistances) |> 
    dplyr::distinct() |> 
    dplyr::full_join(dplyr::summarise(dplyr::group_by(transect_counts, sp),
                                      totTC = sum(true_n),              # total abudance across count sites
                                      totTC_obs = sum(count))) |>       # total number of individuals counted across count sites
    dplyr::full_join(sp_df) 
  
  return(list(data = data,
              constants = constants,
              sp_info = sp_info,
              com_truth = com_truth))
}

#### Model code ####
model.code <- nimble::nimbleCode({
  # uninformative prior for community average of scale parameter intercept for DS detection function
  # this implies that detection probability can be anywhere from 0 to 1 across the width of the surveyed area
  mu_gamma0 ~ dunif(0, 10)
  # typical prior for standard deviation - this is the variation among species for DS detection function scale param intercept
  sd_gamma0 ~ dexp(1)
  # uninformative prior for community average of scale parameter intercept for DS detection function
  # this implies that detection probability can be anywhere from 0 to 1 across the width of the surveyed area
  mu_gamma0_c ~ dunif(0, 10)
  # typical prior for standard deviation - this is the variation among species for DS detection function scale param intercept
  sd_gamma0_c ~ dexp(1)
  # typical weakly informative prior for abundance intercept
  mu_alpha0 ~ dnorm(0, sd = 2) 
  # prior for standard devation among species-level abundance intercepts
  sd_alpha0 ~ dexp(1)          
  # typical weakly informative prior for abundance covariate effect
  mu_alpha1 ~ dnorm(0, sd = 2)
  # prior for standard devation among species-level covariate effects
  sd_alpha1 ~ dexp(1)
  # loop through species
  for ( s in 1:NSPECIES ) {
    # species-level intercept for scale param for DS detection function
    gamma0_ds[s] ~ dnorm( mu_gamma0, sd = sd_gamma0 )
    # species-level intercept for scale param for count detection function
    gamma0_c[s] ~ dnorm( mu_gamma0_c, sd = sd_gamma0_c )
    # scale param for DS detection function
    omega_ds[s] <- exp( gamma0_ds[s] )
    # scale param for count detection function
    omega_c[s] <- exp( gamma0_c[s] )
    # species-level abundance intercept
    alpha0[s] ~ dnorm( mu_alpha0, sd = sd_alpha0 )
    # species-level abundance covariate effect
    alpha1[s] ~ dnorm( mu_alpha1, sd = sd_alpha1 )
    # species-level overall detection probability for DS data
    pie_sp[s] <- sum( pie[1:NBINS, s] )
    # species-level overall detection probability for count data
    pie_sp_c[s] <- sum(pie_c[1:NBINS, s])
    # loop through distance bins - if DS and count data have different number, just create separate loops for each
    for (k in 1:NBINS ) {
      # half normal detection function - DS
      log(g[k, s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega_ds[s] * omega_ds[s])
      # half normal detection function - counts
      log(g_c[k,s]) <- -MIDPOINT[k] * MIDPOINT[k]/(2 * omega_c[s] * omega_c[s])
      # bin-level detection probability for each species - DS
      pie[k, s] <- g[k,s] * (V/B)
      # bin-level detection probability for each species - counts
      pie_c[k,s] <- g_c[k,s] * (V/B)
      # cell probabilities for DS data
      pie_cell[k, s] <- pie[k, s] / pie_sp[s]
    }
  }
  # loop through counts for distance sampling data
  for( i in 1:NSURVEYS ) {
    # expected abundance - log-linear regression
    log( lambda[i] ) <- alpha0[ SP_NG[i]] + alpha1[ SP_NG[i]] * HAB_DS[i] 
    # latent abundance at distance sampling sites!
    N_DS[i] ~ dpois( lambda[i] )
    # observed count at distance sampling sites
    yN_DS[i] ~ dbin( pie_sp[SP_NG[i]], N_DS[i])
  }
  # loop through distance observations
  for (i in 1:NDISTANCES ) {
    # distance class observations modeled with categorical distribution w/ cell probabilities
    DCLASS[i] ~ dcat(pie_cell[1:NBINS, SP_GS[i] ] )
  }
  # loop through count data
  for(i in 1:NCOUNTS) {
    # expected abundance at count sites. Notice that alpha0 and alpha1 appear here again!!
    log(lambda_tc[i]) <- alpha0[SP_TC[i]] + alpha1[SP_TC[i]] * HAB_TC[i]
    # latent abundance at count sites
    N_TC[i] ~ dpois( lambda_tc[i] )
    # observed count at count sites
    yN_TC[i] ~ dbin( pie_sp_c[SP_TC[i]], N_TC[i] )
  }
})

#### MCMC settings & simulation scenarios ####
# parameters to track
params <- c(
  "mu_gamma0", 
  "sd_gamma0",
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
  "N_DS", 
  "N_TC")

#function to generate random initial values for MCMC chains
make_inits <- function(data, constants){
  mu_gamma0_st <- rnorm(1, 5.5, 0.2)
  sd_gamma0_st <- runif(1, 0.1, 0.5)
  mu_gamma0_c_st <- rnorm(1, 5.5, 0.2)
  sd_gamma0_c_st <- runif(1, 0.1, 0.5)
  mu_alpha0_st <- rnorm(1, 0,  2)
  sd_alpha0_st <- runif(1, 1, 2)
  mu_alpha1_st <- rnorm(1, 0, 2)
  sd_alpha1_st <- runif(1, 0.3, 0.75)
  alpha0_st <- alpha1_st  <- gamma0_ds_st <- gamma0_c_st <- numeric(length = constants$NSPECIES)
  gamma0_ds_st <- rnorm( constants$NSPECIES, mu_gamma0_st, sd_gamma0_st )
  gamma0_c_st <- rnorm( constants$NSPECIES, mu_gamma0_c_st, sd_gamma0_c_st )
  alpha0_st <- rnorm( constants$NSPECIES, mu_alpha0_st, sd_alpha0_st )
  alpha1_st <- rnorm( constants$NSPECIES, mu_alpha1_st, sd_alpha1_st )
  inits <- list(
    mu_gamma0 = mu_gamma0_st,
    sd_gamma0 = sd_gamma0_st,
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
    N_DS = data$yN_DS + 1,
    N_TC = data$yN_TC + 1)
  return(inits)
}

# burn-in
nburn <- 100000
# total number of iterations
ni <- nburn + 100000
# thinning interval
nt <- 100
# number of MCMC chains
nc <- 3

# I split this script into multiple jobs on the supercomputer, running 100 at a time
# 1 replicate of the model took ~ 30-40 minutes with the above settings
min_simrep <- 1
max_simrep <- 1000

simrep_rank <- rank(min_simrep:max_simrep)
simrep_raw <- min_simrep:max_simrep

# loop through replicate simulations
for( i in min(simrep_rank):max(simrep_rank)){
  # simulate data
  simdat <- sim_icm()
  data <- simdat$data
  constants <- simdat$constants
  sp_info <- simdat$sp_info
  com_truth <- simdat$com_truth
  print(paste( "Starting rep", simrep_rank[i], "of", max(simrep_rank))) 
  # run NIMBLE model MCMC chains in parallel
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
  
  # summarize replicate output
  outsum <- MCMCvis::MCMCsummary( out ) |> 
    as_tibble(rownames = "param")
  
  # make a dataframe with summary of model estimates and true parameter values
  res <- sp_info |> 
    pivot_longer(c("gamma0_ds", 
                   "gamma0_c",
                   "alpha0",
                   "alpha1"), 
                 names_to = "param", values_to = "truth") |> 
    mutate(param = paste0(param, '[', sp, ']')) |> 
    dplyr::select(param, sp, totDS, totDS_obs, ndistances, totTC, totTC_obs, truth) |> 
    full_join(com_truth) |> 
    full_join(
      full_join( dplyr::select( sp_info, sp, totDS, totDS_obs, ndistances, totTC, totTC_obs),
                 tibble::tibble(
                   sp = constants$SP_NG,
                   param = paste0("N_DS[", 1:length(data$true_n_ds), "]"),
                   truth = data$true_n_ds)
      )
    ) |> 
    full_join(
      full_join( dplyr::select( sp_info, sp, totDS, totDS_obs, ndistances, totTC, totTC_obs),
                 tibble::tibble(
                   sp = constants$SP_TC,
                   param = paste0("N_TC[", 1:length(data$true_n_tc), "]"),
                   truth = data$true_n_tc)
      )
    ) |> 
    left_join(outsum) |> 
    add_column(simrep = simrep_raw[i])
  
  # write out csv of results for the replicate
  write_csv(res, paste0("simrep_no_overdispersion_", formatC(simrep_raw[i], width = 4, format = "d", flag = "0"), "_results.csv"))
  
  print(paste("Rep", simrep_rank[i], "took", round(time[[1]], 3), "hours"))
  # clean up before it all starts over again
  rm( cl, com_truth, constants, data, init, out, outsum, res, simdat, sp_info, end, start, time)
}
