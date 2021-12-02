###############################################################################
#                                                                             #
#   Function to specify the priors of the model                              #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin       #
#   *contact                                                                  #
###############################################################################

## Objects generated
# list_params: List of parameter names
# Hyperparameters: list containing the hyperparameter values to parametrize the prior distributions..
# density(pars): function to calculate the prior likelihood of the parameters
# sampler(n): function to sampler the parameter values (n number of parameter vectors)
# Bounds: data frame containing the boundaries of the parameters.

# The script depends on the two model markers : 
## abio (logical) indicates if the abiotic model is modeled (carrying capacities are not constant across species and plots)
## trait_biotics : indicates which traits are used to calibrate the interaction matrix
library(extraDistr)
library(dplyr)

# Hyperparameters for the alpha interaction matrix priors
hyperparameters <- list(intercept_median = 0.5,
                        intercept_sd = 0.8,
                        
                        mu_mean = 0,
                        mu_sd = 1.8,
                        mu_max = 3,                      

                        sigma_median = 1,#Single trait
                        sigma_sd = 1.4,
                        sigma_min = 0.1, 
                        sigma_max = 10,
                          
                        rho_mean = 0,
                        rho_sd = 1,
                        rho_max = 0.9,
                        
                        # carrying capacity
                        c2_min = 0.1,
                        c2_max = 1.5,
                        
                        a2_min = 0.01,
                        a2_max = 0.99,
                        
                        # Likelihood
                        phi_sd = 2.5,
                        phi_min = 0, 
                        phi_max = Inf)

density <- function(pars){
  names(pars)<- list_params
  
  if (traits_biotic == "none"){
    dbio <- 0
  }else{
    d1 <- dlnorm(pars["intercept"],log(hyperparameters$intercept_median), hyperparameters$intercept_sd, log = T)
    dbio <- d1
  }
  if (grepl("2tr|3tr", traits_biotic) | traits_biotic %in% c("Height", "SLA", "poros")){
    d2 <- dtnorm(pars[grep("mu", names(pars))], hyperparameters$mu_mean, hyperparameters$mu_sd, a = -hyperparameters$mu_max, b = hyperparameters$mu_max, log = T) 
    d3 <- dlnorm(pars[grep("sigma", names(pars))], log(hyperparameters$sigma_median), sdlog = hyperparameters$sigma_sd, log = T)
    dbio <- sum(c(dbio, d2,d3))
  }
  if (grepl("2tr|3tr", traits_biotic)){
      d6 <- dtnorm(pars[grep("rho", names(pars))], mean = hyperparameters$rho_mean, sd = hyperparameters$rho_sd, a = -hyperparameters$rho_max, b = hyperparameters$rho_max, log = T)
      dbio <- dbio+sum(d6)
  }
  #### Carrying capacity params
  da2 <- dunif(pars["a2"],  hyperparameters$a2_min, hyperparameters$a2_max, log = T)
  dphi <- dhcauchy(pars["phi"], hyperparameters$phi_sd, log = TRUE)
  if (any(list_params == "c2")){
    dc2 <- dunif(pars["c2"], hyperparameters$c2_min, hyperparameters$c2_max, log = T)
    dll <- da2+dc2+dphi
  }else{
    dll <- da2+dphi
  }
  dll + dbio
}

sampler <- function(n=1){
  
  pars <- matrix(NA, nrow = n, ncol = length(list_params), dimnames = list(NULL, list_params))
  if (traits_biotic != "none"){
    pars[,"intercept"] <-  rlnorm(n, log(hyperparameters$intercept_median), hyperparameters$intercept_sd)
  }
  if (traits_biotic  %in% c("Height", "SLA", "poros") | grepl("2tr|3tr", traits_biotic)){
    pars[,grep("mu", list_params)] <- rtnorm(n*sum(grepl("mu", list_params)), hyperparameters$mu_mean, hyperparameters$mu_sd, a = -hyperparameters$mu_max, b = hyperparameters$mu_max) 
    pars[,grep("sigma", list_params)] <- rlnorm(n*sum(grepl("sigma", list_params)), meanlog = log(hyperparameters$sigma_median), sdlog = hyperparameters$sigma_sd) 
  }
  if (grepl("2tr|3tr", traits_biotic)){
    pars[,grep("rho", list_params)] <- rtnorm(n*sum(grepl("rho", list_params)), mean = hyperparameters$rho_mean, sd = hyperparameters$rho_sd, a = -hyperparameters$rho_max, b = hyperparameters$rho_max)
  }
  
  #### Carrying capacity params
  pars[,"a2"] <- runif(n, hyperparameters$a2_min, hyperparameters$a2_max)
  if (any(list_params == "c2")){
    pars[,"c2"] <- runif(n, hyperparameters$c2_min, hyperparameters$c2_max)
  } 
  pars[,"phi"] <- rhcauchy(n, hyperparameters$phi_sd)

  return(pars)
}

#### Generation of the object list_params
Ntraits <- case_when(traits_biotic %in% c("Height", "SLA", "poros") ~ 1,
          grepl('2tr', traits_biotic) ~ 2,
          grepl('3tr', traits_biotic) ~ 3,
          TRUE ~ 0)
if (Ntraits > 0){
  list_params_bio <- c("intercept", 
                       paste0(rep("mu", Ntraits), 1:Ntraits), 
                       paste0(rep("sigma", Ntraits), 1:Ntraits))
  if (Ntraits > 1){
    N_rhos <- 0.5*Ntraits*(Ntraits -1)
    list_params_bio <- c(list_params_bio, paste0(rep("rho", N_rhos), 1:N_rhos))
  }
} else if (traits_biotic == "NoTr"){
  list_params_bio <- c("intercept")
} else {
  list_params_bio <- NULL
}

if (abio){
  list_params <- c(list_params_bio, "a2","c2", "phi")
}else{
  list_params <- c(list_params_bio, "a2", "phi")
}
rm(list = c("Ntraits", "list_params_bio"))
#############

### Generation of the object "bounds"
bounds <- as.data.frame(matrix(ncol = 2, nrow = length(list_params), dimnames = list(list_params, c("lower", "upper"))))
if (traits_biotic != "none"){
  bounds["intercept",] <-  c(0.05, 5)
}
if (traits_biotic  %in% c("Height", "SLA", "poros") | grepl("2tr|3tr", traits_biotic)){
  bounds[grep("mu", list_params),"lower"] <- -hyperparameters$mu_max
  bounds[grep("mu", list_params),"upper"] <- hyperparameters$mu_max
  
  bounds[grep("sigma", list_params),"lower"] <- hyperparameters$sigma_min
  bounds[grep("sigma", list_params),"upper"] <- hyperparameters$sigma_max
  
}
if (grepl("2tr|3tr", traits_biotic)){
    bounds[grep("rho", list_params),"lower"] <- -hyperparameters$rho_max
    bounds[grep("rho", list_params),"upper"] <- hyperparameters$rho_max

}
bounds["a2",] <-  c(hyperparameters$a2_min, hyperparameters$a2_max)
if (abio) bounds["c2",] <- c(hyperparameters$c2_min, hyperparameters$c2_max)
bounds["phi",] <- c(hyperparameters$phi_min, hyperparameters$phi_max)
###################
