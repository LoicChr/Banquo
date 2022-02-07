################################################################
#                                                              #
#                   Running of Banquo model                    #
#                                                              #
################################################################

# To specify which assembly model should be run, the script uses 'model markers' (character length one object). 
# The first model markers is abio = c(TRUE, FALSE), a logical. It indicates if the model takes into account the abiotic filtering. 
## Within the code, it serves to specify to the parameter b is not necessary because it is fixed to 0.
# The second model marker is traits_biotic, a length one character vector that serves to specify the parameters to be estimated.
## It serves to specify the traits to be used to compute the pairwise interaction matrix and specify the parameters to be estimated.
## It can take the following values:
### 'none': no interaction (null model and abiotic model)
### 'NoTr': the interspecific coefficients are fixed to a constant.
### 'Height', 'SLA', 'poros': the interactions are estimated with a single trait
### 2tr_HS, 2tr_HP, 2tr_PS: the interactions are estimated with two traits indicated by the capital letters. S 'SLA', P 'root Porosity', H 'Height'
### 3tr: the interactions are estimated with the three traits.

# The script further uses parallelisation to run in parallel all models specify in the params objects. We recommend adaptating this to the architecture of your computing system.


#Librairies
library(BayesianTools)
library(foreach)
library(doParallel)
library(DEoptim)

# Retrieve model characteristics
params <- expand.grid(traits_biotic = c("none","NoTr", "Height","SLA", "poros", "2tr_HS","2tr_HP","2tr_PS","3tr"),
                      abio = c(TRUE, FALSE), stringsAsFactors = F)



wrapper_model <- function(id){
  #Load functions
  source("lib/interactionMatrix.R")
  source("lib/traitspace_and_banquo.R")
  
  # Retrieve model characteristics  
  abio = params[id, "abio"]
  traits_biotic = params[id, "traits_biotic"]
  
  result_file = paste0("results/abio_", abio,"/",traits_biotic)
  if (!file.exists(result_file)) dir.create(result_file, recursive = T)

  # Data
  # Beginning preparation ##############le
  ## Load traitspace object
  load("data/traitspace_objs.Rdata")
  cover_class <- read.table("data/cover_class.txt") # Load cover classes
  
  if (traits_biotic == 'Height'){
    tr <- cbind(t.avg$maxht)
  }
  if (traits_biotic == 'SLA'){
    tr <- cbind(t.avg$sla)
  }
  if (traits_biotic == 'poros'){
    tr <- cbind(t.avg$poros)
  }
  if (traits_biotic == "2tr_HS"){
    tr <- t.avg[,c("maxht", "sla")]
  }
  if (traits_biotic == "2tr_HP"){
    tr <- t.avg[,c("maxht", "poros")]
  }
  if (traits_biotic == "2tr_PS"){
    tr <- t.avg[,c("poros", "sla")]
  }
  if (traits_biotic == "3tr"){
    tr <- t.avg[,c("poros", "sla", "maxht")]
  }
  
  ## Load priors and likelihoods
  source("main/priors.R", local = T)
  source("lib/Likelihood.R", local = T)
  bounds[bounds ==Inf] <- 10

  # Generate the initial Zmatrix using DEoptim
  negLikelihod <- function(pars){
    -1*likelihoodAb(pars)
  }
  itermax <- 1000
  out.DEopt1 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = itermax))
  out.DEopt2 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = itermax))
  out.DEopt3 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = itermax))
  out.DEopt4 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = itermax))
  
  
  # Setup of the mcmc
  prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
  bayesianSetup <- createBayesianSetup(likelihoodAb, prior, names = list_params)
  
  # # settings for the sampler
  source("main/priors.R", local = T)
  newZ <- rbind(out.DEopt1$member$bestmemit, out.DEopt2$member$bestmemit, out.DEopt3$member$bestmemit, out.DEopt4$member$bestmemit)
  newZ <- newZ[sample(1:nrow(newZ)),]
  settings <- list(iterations = 3*4e5, nrChains = 4, Z= newZ) 

  out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings, sampler = 'DEzs') # Run the mcmc

  # # # Compute posterior and basic statistics
  thres <- 3e5
  thin <- 50 
  dis <- getSample(out, parametersOnly = F, start = thres, thin = thin)

  dic.val <- DIC(out, start = thres)
  conv.val <- gelmanDiagnostics(out, start = thres)

  list.objs <- c("dis", "conv.val", "dic.val", "abio", "traits_biotic", "list_params", "tr")
  
  save(list = list.objs, file = paste0(result_file, "/posterior_objs.Rdata"))
}

cl <- makeCluster(15)
registerDoParallel(cl)
foreach (id = 1:nrow(params), .packages = c('BayesianTools', 'extraDistr', 'corpcor', 'DEoptim'), .export = c("params")) %dopar% wrapper_model(id)
