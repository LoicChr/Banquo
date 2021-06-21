likelihoodCov <- function(obs, pred, phi,  offset = 0.005, cover_class){
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
  # This likelihood function is used to choose the option to calculate species abundance across sites

  # Transformation to carrying capacities
  if (abio){
    P_S_E_tr2 <- pars["a2"]*(P_S_E_tr/max(P_S_E_tr))^pars["c2"]
  }else{
    P_S_E_tr2 <- pars["a2"]*(P_S_E_tr/max(P_S_E_tr))^0
  }
  
  # Biotic interaction model
  if (traits_biotic == "none"){
    pred <- P_S_E_interactions <- P_S_E_tr2
  }
  if (traits_biotic == "NoTr"){
    P_S_E_interactions <- banquo(P_S_E_tr2, tr = NULL, avg.out = T, intercept = pars["intercept"],
                                 mu = NULL, sigma = NULL, intra = 1, std = 'FALSE', out_alpha = spxp)
    pred <- as.matrix(P_S_E_interactions)
    if (spxp){
      pred <- as.matrix(P_S_E_interactions[[1]])
    }else{
      pred <- as.matrix(P_S_E_interactions)
    }
  }
  if (traits_biotic %in% c("Height", "SLA")){
    P_S_E_interactions <- banquo(P_S_E_tr2, as.data.frame(trait), avg.out = T, intercept = pars["intercept"],
                                 mu = pars["mu"], sigma = pars["sigma"], intra = 1, std = std, out_alpha = spxp)
    pred <- as.matrix(P_S_E_interactions)
    if (spxp){
      pred <- as.matrix(P_S_E_interactions[[1]])
    }else{
      pred <- as.matrix(P_S_E_interactions)
    }
  }
  if (traits_biotic == "2tr_norm"){
    P_S_E_interactions <- banquo(P_S_E_tr2, tr= cbind(trait1, trait2), avg.out = T,
                                 intercept = pars["intercept"],
                                 mu = c(pars["mu1"], pars["mu2"]),
                                 sigma = c(pars["sigma1"],pars["sigma2"]),
                                 rho = ifelse(rho_pars, pars["rho"],0),
                                 intra =1, std = std, out_alpha = spxp)
    if (spxp){
      pred <- as.matrix(P_S_E_interactions[[1]])
    }else{
      pred <- as.matrix(P_S_E_interactions)
    }

  }
  #Computation of the likelihood
  loglik <- likelihoodCov(comm.red, pred, phi = pars["phi"], offset = 0.005, cover_class = cover_class)

  if (spxp){
   out <- list(loglik, P_S_E_interactions)
  }else{
    out <- loglik
  }

  return(out)
}