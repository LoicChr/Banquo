

#Librairies
library(plyr)
library(ade4)
library(extraDistr)
library(BayesianTools)

#Load functions
source("lib/interactionMatrix.R")
source("lib/traitspace_and_banquo.R")

# Retrieve model characteristics
params <- expand.grid(traits_biotic = c("none","NoTr", "Height", "SLA","2tr", "2tr_rho"),
                      abio = c(TRUE, FALSE), stringsAsFactors = F)


id <- as.numeric(commandArgs(trailingOnly = TRUE))

abio = params[id, "abio"]
traits_biotic = params[id, "traits_biotic"]
rho_pars <- params[id, "rho_pars"]

result_file = paste0("results/abio_", abio,"/",traits_biotic)
if (!file.exists(result_file)) dir.create(result_file, recursive = T)

# Data
# Beginning preparation ##############
## Load traitspace object
load(XXXXXX)
load(XXXX) # load biotic trait values
cover_class <- read.table("data/cover_class.txt") # Load cover classes



comm.red <- obs.comm[row.names(P_S_E_tr), colnames(P_S_E_tr)]
comm.red <- comm.red[, colSums(comm.red) > 0]
P_S_E_tr <- P_S_E_tr[, colnames(comm.red)]
row.names(t.avg) <- t.avg[,1]
t.avg <- t.avg[colnames(P_S_E_tr),]

if (traits_biotic == 'Height'){
  trait <- t.avg$maxht
}
if (traits_biotic == 'SLA'){
  trait <- t.avg$sla
}
if (grepl("2tr", traits_biotic)){
  trait1 <- t.avg$Axis1
  trait2 <- t.avg$Axis2
}

## Load priors and likelihoods
source("main/priors.R")
source("lib/Likelihood_Ab.R")

# Setup of the mcmc
prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
bayesianSetup <- createBayesianSetup(likelihoodAb, prior, names = list_params)
# # settings for the sampler
settings <- list(iterations = DDDD, nrChains = 4) #50000 per chains in the one chain triplet

out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings) # Run the mcmc

save(out, file = paste0(result_file, "/chain1.Rdata")) #Save the results
