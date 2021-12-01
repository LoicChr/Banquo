likelihoodCov <- function(obs, pred, phi,  offset = 0.0005, cover_class){
  # This likelihood function actually calculates the likelihood of the mean covers according to a beta distribution.
  
  obs[obs == 0] <- offset
  mean_cover <- as.numeric(pred)
  mean_cover[mean_cover < offset] <- offset
  alpha = mean_cover * phi #Conversion in alpha/beta shape parameters
  beta = (1-mean_cover) * phi
  j <- match(obs, cover_class$cover_class)
  
  prob_cat <- log(pbeta(cover_class$maxCov[j], alpha, beta) - pbeta(cover_class$minCov[j], alpha, beta))
  logLik <- matrix(prob_cat, nrow = nrow(obs), ncol = ncol(obs), byrow = F)
  
  
  return(sum(logLik))
}


likelihoodAb <-  function(pars, spxp = F){
  names(pars) <- list_params
  # This likelihood function uses the two model markers (abio (logical), traits_biotic (character)) to choose how calculate species abundance across plots.
  # The option spxp is used in a non-MCMC computation context to output the predicted species covers.
  
  # Transformation to carrying capacities
  if (abio){
    P_S_E_tr2 <- pars["a2"]*(P_S_E_tr/max(P_S_E_tr))^pars["c2"]
  }else{
    P_S_E_tr2 <- pars["a2"]*(P_S_E_tr/max(P_S_E_tr))^0
  }
  
  # Biotic interaction model
  if (traits_biotic == "none"){
    pred <- P_S_E_interactions <- P_S_E_tr2
  }  else if (traits_biotic == "NoTr"){
    P_S_E_interactions <- banquo(P_S_E_tr2, tr = NULL, intercept = pars["intercept"],
                                 mu = NULL, sigma = NULL, out_alpha = spxp)
    if (spxp){
      pred <- as.matrix(P_S_E_interactions[[1]])
    }else{
      pred <- as.matrix(P_S_E_interactions)
    }
  } else {
    P_S_E_interactions <- banquo(P_S_E = P_S_E_tr2, tr= tr, 
                                 intercept = pars["intercept"],
                                 mu = pars[grep("mu", list_params)],
                                 sigma = pars[grep("sigma", list_params)],
                                 rho = pars[grep("rho", list_params)],
                                 out_alpha = spxp, intra = 1)
    if (spxp){
      pred <- as.matrix(P_S_E_interactions[[1]])
    }else{
      pred <- as.matrix(P_S_E_interactions)
    }
    
  }
  #Computation of the likelihood
  loglik <- likelihoodCov(obs = obs.comm, pred, phi = pars["phi"], offset = 0.0005, cover_class = cover_class)
  if (spxp){
    out <- list(loglik, P_S_E_interactions)
  }else{
    out <- loglik
  }
  
  return(out)
}