###############################################################################
#                                                                             #
#   Parametrization of the Banquo model                                       #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin        #
#   *contact                                                                  #
###############################################################################
### ML_interactions: computation of a single trait Banquo model
### ML_interactions_bi and ML_interactions_bi2: computation of a two trait Banquo model. 
###               These functions differs in their way to parametrize the peak position.
###               ML_interactions_bi computes it in cartesian coordinates.
###               ML_interactions_bi2 computes it in polar coordinates. 
###               This choice affects the bounds of the space in which the nloptr algorithm searchs 
###               for a solution. The first one looks into a rectangle determined by trait1 and trait2 range. The second in an ellipse.

### Parameter controlling the interaction matrix are transformed prior to the banquo computation. 
#### The internal function params_tr specifies the transformation

#### The three functions returns the output in the same format. It is a list with the following elements:
######## $pars.tr: the final parametrization of the interaction matrix (after parameter transformation)
######## $negll: the final value of the negative likelihood
######## $output_nlptor: the output as provided  by the nloptr function
######## $P_S_E_interactions: the predicted species abundances in each communities
######## $alphas: the interaction matrix as parametrized by pars.tr

ML_interactions <- function(P_S_E_tr, trait, obs.comm, lb = c(0, -1, 0.04),  ub = c(20, 1, 15),
                            opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1000, xtol_rel = 10e-4)){
  require(nloptr)
  require(topicmodels)
  ### Arguments
  # P_S_E_tr: the traitspace object
  # trait: the trait to use to calibrate Banquo
  # lb: lower bounds of the parameter intra, mu and sigma.
  # up: upper bounds of the parameter intra, mu and sigma.
  # opts: arguments for the function 'nloptr'.
  
  list_params <- c("intra", "mu", "sigma")
  rg.trdiff <- max(range(outer(scale(trait), scale(trait), "-")))
  
  # Transformation of parameters 
  params_tr <- function(pars){
    pars.tr <- pars
    
    pars.tr["intra"] <- exp(pars["intra"])
    pars.tr["mu"] <-rg.trdiff*pars["mu"]
    pars.tr["sigma"] <- exp(pars["sigma"])
    return(pars.tr)
  }

  # Likelihood function
  Neglikelihood <- function(pars){
    names(pars) <- list_params
    
    # Transformation of the parameters
    pars.tr <- params_tr(pars)
   
    P_S_E_interactions <- banquo(P_S_E_tr, as.data.frame(trait), avg.out = T, intercept = 1, 
                                 mu = pars.tr["mu"], sigma = pars.tr["sigma"], intra = pars.tr["intra"])
    
    # Likelihood of the output
    comm.red <- obs.comm[row.names(P_S_E_interactions), colnames(P_S_E_interactions)]
    comm.red <- as.matrix(sweep(comm.red, 1, rowSums(comm.red), "/"))
    P_S_E_interactions <- as.matrix(P_S_E_interactions)

    #Computation of the likelihood
    neglogLik <-  -sum(sapply(1:nrow(comm.red), function(j){ 
      distHellinger(t(as.matrix(comm.red[j,])), t(as.matrix(P_S_E_interactions[j,])))
    }))
    
    
    return(neglogLik)
  }
  pars.init = sapply(1:length(lb), function(i) runif(1, lb[i], ub[i]))
  
  out = nloptr(pars.init, Neglikelihood,  lb = lb, ub = ub, 
               opts = opts)
  pars = out$solution
  names(pars)  <- list_params
  pars.tr <- params_tr(pars)

  pred <- banquo(P_S_E_tr, as.data.frame(trait), avg.out = T, 
                               intercept = 1 , 
                               mu = pars.tr["mu"], sigma = pars.tr["sigma"], intra = pars.tr["intra"],out_alpha = TRUE)
  P_S_E_interactions <- pred[[1]]
  alphas <- pred[[2]]
  
  res <- list(pars.tr = pars.tr, negll = out$objective, output_nlptor = out, P_S_E_interactions = P_S_E_interactions, alphas = alphas)
  return(res)
}


ML_interactions_bi_cart <- function(P_S_E_tr, trait1, trait2, obs.comm, lb = c(0, -1, 0.04,-1, 0.04,-1),ub = c(20, 1, 15, 1, 15, 1), det_lim = NULL,
                               opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1000, xtol_rel = 10e-4)){
  require(nloptr)
  require(topicmodels)
  ### Arguments
  # P_S_E_tr: the traitspace object
  # trait: the trait to use to calibrate Banquo
  # lb: lower bounds of the parameter intra, mu1, sigma1, mu2, sigma2, rho.
  # up: upper bounds of the parameter intra, mu1, sigma1, mu2, sigma2, rho.
  # det_lim: optional lower bound for the quantity sigma1*sigma2*sqrt(1-rho^2)
  # opts: arguments for the function 'nloptr'.
  
  list_params <- c("intra", "mu1", "sigma1", "mu2", "sigma2", "rho")
 
  rg.trdiff1 <- max(range(outer(scale(trait1), scale(trait1), "-")))
  rg.trdiff2 <- max(range(outer(scale(trait2), scale(trait2), "-")))
  
  params_tr <- function(pars, det_lim = NULL, sigma2_lb, sigma2_ub){
    pars.tr <- pars
      
    pars.tr["intra"] <- exp(pars["intra"])
    pars.tr["mu1"] <-rg.trdiff1*pars["mu1"]
    pars.tr["mu2"] <-rg.trdiff2*pars["mu2"]
    
    pars.tr["sigma1"] <- exp(pars["sigma1"])
      
      if (!is.null(det_lim)){
        if (det_lim/pars.tr["sigma1"] > exp(sigma2_lb)){
          std <- (pars["sigma2"] - sigma2_lb)/(sigma2_ub- sigma2_lb)
          new_low <- det_lim/pars.tr["sigma1"]
          pars.tr["sigma2"] <- exp(std*(sigma2_ub - new_low) + new_low)
        }else{
          pars.tr["sigma2"] <- exp(pars["sigma2"])
        }
      }else{
        pars.tr["sigma2"] <- exp(pars["sigma2"])
      }
      ratio <- max(0, 1 - (det_lim/(pars.tr["sigma1"]*pars.tr["sigma2"]))^2)
      pars.tr["rho"] <- pars['rho']*max(0, sqrt(ratio))
      return(pars.tr)
  }
  
  Neglikelihood <- function(pars){
    names(pars) <- list_params
    
    # Transformation of the parameters
    pars.tr <- params_tr(pars, det_lim = det_lim, sigma2_lb = lb[5], sigma2_ub = ub[5])

    #Banquo run
    P_S_E_interactions <- banquo(P_S_E_tr, tr= cbind(trait1, trait2), avg.out = T, 
                                 intercept = 1, 
                                 mu = c(pars.tr["mu1"], pars.tr["mu2"]), 
                                 sigma = c(pars.tr["sigma1"],pars.tr["sigma2"]),
                                 rho = pars.tr["rho"],
                                 intra =pars.tr["intra"]) 
    
    # Likelihood of the output
    comm.red <- obs.comm[row.names(P_S_E_interactions), colnames(P_S_E_interactions)]
    comm.red <- as.matrix(sweep(comm.red, 1, rowSums(comm.red), "/"))
    P_S_E_interactions <- as.matrix(P_S_E_interactions)
    
    #Computation of the likelihood
    neglogLik <-  -sum(sapply(1:nrow(comm.red), function(j){ 
      distHellinger(t(as.matrix(comm.red[j,])), t(as.matrix(P_S_E_interactions[j,])))
    }))
    
    return(neglogLik)
  }
  pars.init = sapply(1:length(lb), function(i) runif(1, lb[i], ub[i]))
  
  out = nloptr(pars.init, Neglikelihood,  lb = lb, ub = ub, opts = opts)
  pars = out$solution
  names(pars)  <- names(pars.init) <- list_params
  pars.tr <- params_tr(pars, det_lim = det_lim,sigma2_lb = lb[5], sigma2_ub = ub[5])

  pred <- banquo(P_S_E_tr, cbind(trait1, trait2), avg.out = T, 
                               intercept = 1 , 
                               mu = c(pars.tr["mu1"], pars.tr["mu2"]), 
                               sigma = c(pars.tr["sigma1"],pars.tr["sigma2"]),
                               rho =  pars.tr["rho"], intra = pars.tr["intra"], out_alpha = TRUE)
  P_S_E_interactions <- pred[[1]]
  alphas <- pred[[2]]

  res <- list(pars.tr = pars.tr, negll = out$objective, out = out, P_S_E_interactions = P_S_E_interactions, alphas = alphas)
  return(res)
}

ML_interactions_bi_pol <- function(P_S_E_tr, trait1, trait2, obs.comm, lb = c(0, 0, 0.04,0, 0.04,-1),ub = c(20, 1, 15, pi*2, 15, 1), 
                               opts = list(algorithm = "NLOPT_LN_SBPLX", maxeval = 1000, xtol_rel = 10e-4), det_lim = NULL){
  require(nloptr)
  ### Arguments
  # P_S_E_tr: the traitspace object
  # trait: the trait to use to calibrate Banquo
  # lb: lower bounds of the parameter intra, r_dis, sigma1, theta, sigma2, rho.
  # up: upper bounds of the parameter intra, r_dis, sigma1, theta, sigma2, rho.
  # det_lim: optional lower bound for the quantity sigma1*sigma2*sqrt(1-rho^2)
  # opts: arguments for the function 'nloptr'.
  
  list_params <- c("intra", "r_dis", "sigma1", "theta", "sigma2", "rho")

  rg.trdiff1 <- max(range(outer(scale(trait1), scale(trait1), "-")))
  rg.trdiff2 <- max(range(outer(scale(trait2), scale(trait2), "-")))
  
  params_tr <- function(pars, det_lim = NULL, sigma2_lb, sigma2_ub){
    pars.tr <- pars
    names(pars.tr)[c(2,4)] <- c("mu1", "mu2")
    
    pars.tr["intra"] <- exp(pars["intra"])
    pars.tr["sigma1"] <- exp(pars["sigma1"])
    
    r_ellipse <- function(theta, a, b){
      a*b/(sqrt((b*cos(theta))^2 +(a*sin(theta))^2))
    }
    
    r_lim <- r_ellipse(pars["theta"], a = rg.trdiff1, b = rg.trdiff2)
    r_new <- sqrt(pars["r_dis"])*r_lim
    pars.tr["mu1"] <- r_new*cos(pars["theta"])
    pars.tr["mu2"] <- r_new*sin(pars["theta"])
    
    if (!is.null(det_lim)){
      if (det_lim/pars.tr["sigma1"] > exp(sigma2_lb)){
        std <- (pars["sigma2"] - sigma2_lb)/( sigma2_ub - sigma2_lb)
        new_low <- det_lim/pars.tr["sigma1"]
        pars.tr["sigma2"] <- exp(std*( sigma2_ub - new_low) + new_low)
      }else{
        pars.tr["sigma2"] <- exp(pars["sigma2"])
      }
    }else{
      pars.tr["sigma2"] <- exp(pars["sigma2"])
    }
    ratio <- max(0, 1 - (det_lim/(pars.tr["sigma1"]*pars.tr["sigma2"]))^2)
    pars.tr["rho"] <- pars['rho']*max(0, sqrt(ratio))
    return(pars.tr)
  }
  
  Neglikelihood <- function(pars){
    names(pars) <- list_params
    # Transformation of the parameters
    pars.tr <- params_tr(pars, det_lim = det_lim, sigma2_lb = lb[5], sigma2_ub = ub[5])
   
    #Banquo run
    P_S_E_interactions <- banquo(P_S_E_tr, tr= cbind(trait1, trait2), avg.out = T, 
                                 intercept = 1, 
                                 mu = c(pars.tr["mu1"], pars.tr["mu2"]), 
                                 sigma = c(pars.tr["sigma1"],pars.tr["sigma2"]),
                                 rho = pars.tr["rho"],
                                 intra = pars.tr["intra"]) 
    
    # Likelihood of the output
    comm.red <- obs.comm[row.names(P_S_E_interactions), colnames(P_S_E_interactions)]
    comm.red <- as.matrix(sweep(comm.red, 1, rowSums(comm.red), "/"))
    P_S_E_interactions <- as.matrix(P_S_E_interactions)
    
    #Computation of the likelihood
    neglogLik <-  -sum(sapply(1:nrow(comm.red), function(j){ 
      distHellinger(t(as.matrix(comm.red[j,])), t(as.matrix(P_S_E_interactions[j,])))
    }))
    
    return(neglogLik)
  }

  pars.init = sapply(1:length(lb), function(i) runif(1, lb[i], ub[i]))
  out = nloptr(pars.init, Neglikelihood,  lb = lb, ub = ub, opts = opts)
  pars = out$solution
  names(pars)  <-list_params
  pars.tr <- params_tr(pars, det_lim = det_lim, sigma2_lb = lb[5], sigma2_ub = ub[5])
  

  out <- banquo(P_S_E_tr, cbind(trait1, trait2), avg.out = T, 
                               intercept = 1 , 
                               mu = c(pars.tr["mu1"], pars.tr["mu2"]), 
                               sigma = c(pars.tr["sigma1"],pars.tr["sigma2"]),
                               rho = pars.tr["rho"], intra = pars.tr["intra"], 
                               out_alpha = TRUE)
  P_S_E_interactions <- pred[[1]]
  alphas <- pred[[2]]
  
  res <- list(pars.tr = pars.tr, negll = out$objective, out = out, P_S_E_interactions = P_S_E_interactions, alphas = alphas)
  return(res)
}
