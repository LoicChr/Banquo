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

load(paste0(result_file, "/chains.Rdata"))
for (i in 1:10){
  dis <- getSample(out)
  dis <- getSample(out, start = nrow(dis)/(3*4)-5e3)
  newZ <- do.call(cbind, lapply(1:ncol(dis), function(i){
    rtnorm(500, mean(dis[,i]), sd(dis[,i]), a = bounds[i,1], b = bounds[i,2])
  }))
  # # settings for the sampler
  settings <- list(iterations = 1e4, nrChains = 4, Z = newZ) #50000 per chains in the one chain triplet
  out <- runMCMC(out, settings = settings) # Run the mcmc
}

save(list = ls(), file = paste0(result_file, "/chains_comp.Rdata")) #Save the results
