#Librairies
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

result_file = paste0("results/abio_", abio,"/",traits_biotic)
if (!file.exists(result_file)) dir.create(result_file, recursive = T)

# Data
# Beginning preparation ##############
## Load traitspace object
load("data/traitspace_objs.Rdata")
cover_class <- read.table("data/cover_class.txt") # Load cover classes


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
source("lib/Likelihood.R")

# Setup of the mcmc
prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
bayesianSetup <- createBayesianSetup(likelihoodAb, prior, names = list_params)
# # settings for the sampler
settings <- list(iterations = 5e5, nrChains = 4) #50000 per chains in the one chain triplet
out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings) # Run the mcmc

# Compute posterior and basic statistics
nstep <- nrow(getSample(out, start = 0, parametersOnly = F, thin = 0))/(4*3)
thres <- nstep - 2e5
dis <- getSample(out, parametersOnly = F, start = thres, thin = 50)

dic.val <- DIC(out, start = thres)$DIC
conv.val <- gelmanDiagnostics(out, start = thres)$mpsrf
rm(out)

save(list = ls(), file = paste0(result_file, "/output.Rdata")) #Save the results
