#Librairies
library(BayesianTools)
library(foreach)
library(doParallel)


# Retrieve model characteristics
params <- expand.grid(traits_biotic = c("none","noTr", "Height","SLA", "poros", "2tr_HS","2tr_HP","2tr_PS","3tr"),
                      abio = c(TRUE, FALSE), stringsAsFactors = F)



wrapper_model <- function(id){
  #Load functions
  source("lib/interactionMatrix.R")
  source("lib/traitspace_and_banquo.R")
  
  abio = params[id, "abio"]
  traits_biotic = params[id, "traits_biotic"]
  
  result_file = paste0("results_noPCA/abio_", abio,"/",traits_biotic)
  if (!file.exists(result_file)) dir.create(result_file, recursive = T)

  # Data
  # Beginning preparation ##############
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
  out.DEopt1 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = 1000))#1000
  out.DEopt2 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = 1000))#1000
  out.DEopt3 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = 1000))#1000
  out.DEopt4 = DEoptim(fn =  negLikelihod, lower = bounds[,1], upper = bounds[,2], control = DEoptim.control(itermax = 1000))#1000
  
  
  # Setup of the mcmc
  prior <- createPrior(density = density, sampler = sampler, lower = bounds[,1], upper = bounds[,2])
  bayesianSetup <- createBayesianSetup(likelihoodAb, prior, names = list_params)
  
  # # settings for the sampler
  source("main/priors.R", local = T)
  newZ <- rbind(out.DEopt1$member$bestmemit, out.DEopt2$member$bestmemit, out.DEopt3$member$bestmemit, out.DEopt4$member$bestmemit)
  newZ <- newZ[sample(1:nrow(newZ)),]
  settings <- list(iterations = 3*4e5, nrChains = 4, Z= newZ) #50000 per chains in the one chain triplet

  out <- runMCMC(bayesianSetup = bayesianSetup, settings = settings, sampler = 'DEzs') # Run the mcmc

  # # # Compute posterior and basic statistics
  thres <- 3e5
  dis <- getSample(out, parametersOnly = F, start = thres, thin = 50)

  dic.val <- DIC(out, start = thres)
  conv.val <- gelmanDiagnostics(out, start = thres)


  dis <- getSample(out, start = 3e5, parametersOnly = F)
  list.objs <- c("dis", "conv.val", "dic.val", "tr", "abio", "traits_biotic", "list_params")
  
  save(list = list.objs, file = paste0(result_file, "/posterior_objs.Rdata"))
}

cl <- makeCluster(15)
registerDoParallel(cl)
foreach (id = 1:nrow(params), .packages = c('BayesianTools', 'extraDistr', 'corpcor', 'ade4', 'DEoptim'), .export = c("params")) %dopar% wrapper_model(id)
