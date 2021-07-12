library(extraDistr)
# Hyperparameters for the alpha interaction matrix priors
hyperparameters <- list(intercept_median = ifelse(!grepl("2tr", traits_biotic), 2,5),
                        intercept_sd = ifelse(!grepl("2tr", traits_biotic), 2.3,3),
                        
                        mu_mean = 0,
                        mu_sd = 1.8,
                        
                        sigma_median = 1, #Single trait
                        sigma_sd = 3.5,
                        
                        sigma12_median = 1, #Double traits
                        sigma12_sd = 2,
                        
                        rho_mean = 0,
                        rho_sd = 1,
                        rho_max = 0.9,
                        
                        # carrying capacity
                        c2_min = 0.025,
                        c2_max = 1.5,
                        
                        a2_min = 0.01,
                        a2_max = 0.99,
                        # Likelihood
                        phi_sd = 2.5)

density <- function(pars){
  names(pars)<- list_params
  
  if (traits_biotic == "none"){
    dbio <- 0
  }
  if (traits_biotic %in% c("NoTr")){
    d1 <- dlnorm(pars["intercept"],log(hyperparameters$intercept_median), hyperparameters$intercept_sd, log = T)
    dbio <- d1
  }
  if (traits_biotic %in% c("Height", "SLA")){
    d1 <- dlnorm(pars["intercept"],log(hyperparameters$intercept_median), hyperparameters$intercept_sd, log = T)
    d2 <- dnorm(pars["mu"], hyperparameters$mu_mean, hyperparameters$mu_sd, log = T) 
    d3 <- dlnorm(pars["sigma"], meanlog = log(hyperparameters$sigma_median), sdlog = hyperparameters$sigma_sd, log = T)
    dbio <- d1+d2+d3
  }
  if (grepl("2tr", traits_biotic)){
    d1 <- dlnorm(pars["intercept"],log(hyperparameters$intercept_median), hyperparameters$intercept_sd, log = T)
    d2 <- dnorm(pars["mu1"], hyperparameters$mu_mean, hyperparameters$mu_sd, log = T) 
    d3 <- dlnorm(pars["sigma1"], log(hyperparameters$sigma12_median), sdlog = hyperparameters$sigma12_sd, log = T)
    d4 <- dnorm(pars["mu2"], hyperparameters$mu_mean, hyperparameters$mu_sd, log = T) 
    d5 <- dlnorm(pars["sigma2"], meanlog = log(hyperparameters$sigma12_median), sdlog = hyperparameters$sigma12_sd, log = T)
    dbio <- d1+d2+d3+d4+d5
    if (traits_biotic == "2tr_rho"){
      d6 <- dtnorm(pars["rho"], mean = hyperparameters$rho_mean, sd = hyperparameters$rho_sd, a = -hyperparameters$rho_max, b = hyperparameters$rho_max, log = T)
      dbio <- dbio+d6
    }
  }
  #### Carrying capacity params
  da2 <- dunif(pars["a2"],  hyperparameters$a2_min, hyperparameters$a2_max, log = T)
  dphi <- dhcauchy(pars["phi"], hyperparameters$phi_sd, log = TRUE)
  if (abio){
    dc2 <- dunif(pars["c2"], hyperparameters$c2_min, hyperparameters$c2_max, log = T)
    dll <- da2+dc2+dphi
  }else{
    dll <- da2+dphi
  }
  dll + dbio
}

sampler <- function(n=1){
  if (traits_biotic == "none"){
    dbio <- NULL
  }
  if (traits_biotic %in% c("NoTr")){
    d1 <- rlnorm(n, log(hyperparameters$intercept_median), hyperparameters$intercept_sd)
    dbio <- cbind(d1)
  }
  if (traits_biotic %in% c("Height", "SLA")){
    d1 <- rlnorm(n, log(hyperparameters$intercept_median), hyperparameters$intercept_sd)
    d2 <- rnorm(n, hyperparameters$mu_mean, hyperparameters$mu_sd) 
    d3 <- rlnorm(n, meanlog = log(hyperparameters$sigma_median), sdlog = hyperparameters$sigma_sd)
    dbio <- cbind(d1, d2, d3)
  }
  if (grepl("2tr", traits_biotic)){
    d1 <- rlnorm(n, log(hyperparameters$intercept_median),hyperparameters$intercept_sd)
    d2 <- rnorm(n, hyperparameters$mu_mean, hyperparameters$mu_sd) 
    d3 <- rlnorm(n, meanlog = log(hyperparameters$sigma12_median), sdlog = hyperparameters$sigma12_sd)
    d4 <- rnorm(n, hyperparameters$mu_mean, hyperparameters$mu_sd) 
    d5 <- rlnorm(n, meanlog = log(hyperparameters$sigma12_median), sdlog = hyperparameters$sigma12_sd)
    dbio <- cbind(d1, d2, d3, d4, d5)
    if (traits_biotic == "2tr_rho"){
      d6 <- rtnorm(n, mean = hyperparameters$rho_mean, sd = hyperparameters$rho_sd, a = -hyperparameters$rho_max, b = hyperparameters$rho_max)
      dbio <- cbind(dbio, d6)
    }

  }
  #### Carrying capacity params
  da2 <- runif(n, hyperparameters$a2_min, hyperparameters$a2_max)
  dphi <- rhcauchy(n, hyperparameters$phi_sd)
  if (abio){
    dc2 <- runif(n, hyperparameters$c2_min, hyperparameters$c2_max)
    dll <- cbind(dbio, da2,dc2, dphi)
  }else{
    dll <- cbind(dbio, da2, dphi)
  }

  return(dll)
}
## Banquo params
if (traits_biotic %in% c("Height", "SLA")){
  list_params_bio <- c("intercept", "mu", "sigma")
  bounds_bio <- data.frame(lower = c(0.05,-3,0.05), upper = c(Inf,3,10), row.names = list_params_bio)
} else if (grepl("2tr", traits_biotic)){
  list_params_bio <- c("intercept", "mu1", "sigma1", "mu2", "sigma2", 'rho')
  bounds_bio <- data.frame(lower = c(0.05,-3,0.05,-3,0.05,-0.9), upper = c(Inf,3,10,3,10,0.9), row.names = list_params_bio)
  if (traits_biotic != "2tr_rho"){
    list_params_bio <- list_params_bio[-6]
    bounds_bio <- bounds_bio[-6,]
  }
} else if (traits_biotic == "NoTr"){
  list_params_bio <- c("intercept")
  bounds_bio <- data.frame(lower = c(0.05), upper = c(Inf), row.names = list_params_bio)
  
}else{
  list_params_bio <- NULL
  bounds_bio <- NULL
}
if (abio){
  list_params_ll<- c("a2","c2", "phi")
  bounds_ll <- data.frame(lower =  c(hyperparameters$a2_min, hyperparameters$c2_min,0), upper = c(hyperparameters$a2_max,hyperparameters$c2_max,Inf), row.names = list_params_ll)
  }else{
   list_params_ll<- c("a2", "phi")
  bounds_ll <- data.frame(lower =  c(hyperparameters$a2_min,0), upper = c(hyperparameters$a2_max,Inf), row.names = list_params_ll)
}

bounds <- rbind(bounds_bio, bounds_ll)
list_params <- c(list_params_bio, list_params_ll)
rm(list = c("list_params_bio", "list_params_ll", "bounds_bio", "bounds_ll"))