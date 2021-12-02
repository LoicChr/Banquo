###############################################################################
#                                                                             #
#   Banquo and Traitspace functions                                           #
#   From the article: Predictions of biodiversity are improved by             #
#    integrating trait-based competition with abiotic filtering               #
#                                                                             #
###############################################################################


traitspace <- function(trait.model, env, PT.Sk, N = 100, avg.out = TRUE){
  ##### Arguments.
  # trait.model: linear model predicting the trait values as a function of environmental variables.
  # env: data.frame containing the environmental values of the sites under study.
  # PT.Sk: a list containing the distribution of species traits. Here Mclust objects.
  # N: the number of random draws used to the MonteCarlo integration.
  # avg.out: should the random draws be averaged out for the final output?
  
  require(mclust)
  require(mvtnorm)
  
  if (is.null(row.names(env))){
    row.names(env) <- paste("Site", 1:nrow(env), "-")
  }else{
    row.names(env) <- gsub("_", "-", row.names(env))
  }
  
  #Extraction of the trait-environment prediction
  pred_multi <- as.data.frame(predict(trait.model, newdata=env, interval="prediction", level=0.95))
  cov_mat <- cov(as.matrix(residuals(trait.model)))
  if (is.null(trait.model$weights)){
    cov_mat <- cov.wt(as.matrix(residuals(trait.model)))$cov #no weights
    
  }else{
    cov_mat <- cov.wt(as.matrix(residuals(trait.model)), wt = trait.model$weights)$cov #add weights
    
  }
  # Construction of the simulated data
  ## STEP 2a: Drawing samples from P(T/E)
  env.to.fit <- data.frame(env[gl(nrow(env), N),])
  dimnames(env.to.fit) <- list(paste(row.names(env)[gl(nrow(env), N)], rep(1:N, nrow(env)), sep = "_"), colnames(env))
  tr_mean_pred <- pred_multi[gl(nrow(env), N), ,drop=F]
  row.names(tr_mean_pred) <- paste(row.names(pred_multi)[gl(nrow(pred_multi), N)], rep(1:N, nrow(pred_multi)), sep = "_")
  
  if (all(dim(cov_mat) == c(1,1))) tr_mean_pred <- tr_mean_pred[, 1, drop = F]
  
  tr_sample <- as.data.frame(t(apply(tr_mean_pred, 1, function(x) rmvnorm(1, x, cov_mat))))
  if (all(dim(cov_mat) == c(1,1))) tr_sample <- as.data.frame(t(tr_sample))
  
  ## computing(P(T/E))
  P_T_E <-do.call(rbind, lapply(1:nrow(tr_sample), function(j) dmvnorm(as.numeric(tr_sample[j,]),as.numeric(tr_mean_pred[j,]), cov_mat)))
  row.names(P_T_E) <- row.names(tr_sample)
  
  ## Step 2b: Computing the likelihood P(T/Sk) using Mclust done earlier
  P_T_S <- lapply(PT.Sk, function(pdf){
    if (ncol(tr_sample) == 1){
      mclust::dens(pdf$modelName,tr_sample[,1],parameters=pdf$parameters)
    }else{
      mclust::dens(pdf$modelName,tr_sample,parameters=pdf$parameters)
    }
  })
  
  P_T_S <- do.call(cbind, P_T_S)
  colnames(P_T_S) <- names(PT.Sk)
  
  ## Step 2c: Computing posterior P(Sk/T,E) using Bayes theorem
  P_S_T_E <- exp(sweep(log(P_T_S), 1, log(rowSums(P_T_S)), "-"))
  
  ## Step 2d: Posterior P(Sk/T) by integrating out T's (with log)
  P_S_E_all <- exp(sweep(log(P_S_T_E), 1, log(P_T_E), "+"))
  
  row.names(P_S_E_all) <- row.names(P_T_E)
  
  sites <- sapply(strsplit(row.names(P_S_E_all), "_"), function(x) x[1])
  
  ### Step 2d.1: Monte Carlo integration across trait samples
  P_S_E_all.sp <- split(as.data.frame(P_S_E_all), sites)
  P_S_E_traitspace_unnorm <- do.call(rbind, lapply(P_S_E_all.sp, function(x) apply(x, 2, mean)))
  
  return(as.matrix(P_S_E_traitspace_unnorm))
}

banquo <- function(P_S_E, tr = NULL, intercept=0.5, mu=-NULL, sigma = NULL, rho = NULL, out_alpha = FALSE, intra = 1){
  require(corpcor)
  ### Banquo computes species abundances based on the traitspace object, the trait that controls competition 
  #### and the parameters
  ##### Arguments.
  # P_S_E: the output of the traitspace function
  # tr : the trait(s) used to compute the interaction matrix
  # intercept, mu, sigma, intra, rho. The parameters to generate the interaction matrix
  # out_alpha : should the interaction matrix be saved in the output?

  # Compute the pairwise interaction matrix
  if (is.null(mu) & is.null(sigma) & !is.null(intercept)){
    alphas <- matrix(intercept, nrow = ncol(P_S_E), ncol = ncol(P_S_E))
    diag(alphas) <- intra
  }else{
    if (length(mu) == ncol(tr) & length(sigma) == ncol(tr) & ncol(tr) == ncol(tr) & ncol(P_S_E) == nrow(tr)){
      tr = apply(tr, 2, scale)
      alphas <- interaction_matrix(tr = tr,
                                    intercept = intercept, mu = mu, sigma = sigma, rho = rho,
                                    intra = intra)
    }else{
      stop("tr, mu, sigma, rho, P_S_E do not match")
    }
  }
  
  #Compute the inverse of the pairwise interaction matrix
  alphas.inv <- corpcor::pseudoinverse(alphas)
  
  P_S_E_all_interactions <- P_S_E
  
  # Function to calculate species abundance at equilibrium as a function of species carrying capacities (Ks) 
  integration.biotic <- function(Ks){
    Ns <- alphas.inv %*% as.matrix(Ks)
    # at least one species has been competitively excluded
    excluded <- c()
    if(any(Ns < 0)){
      A <- alphas
      Ks_red <- Ks
      excluded <- c(excluded, which.min(Ns))
      Ks_red <- Ks[-excluded]
      A <- A[-excluded, -excluded]
      Ns <- Ks*0
      Ns[-excluded] <- corpcor::pseudoinverse(A) %*% as.matrix(Ks_red)
      
      while(any(Ns < 0)){
        excluded <- c(excluded, which.min(Ns))
        A <- alphas
        Ks_red <- Ks
        Ks_red <- Ks[-excluded]
        A <- as.matrix(A[-excluded, -excluded])
        Ns <- Ks*0
        Ns[-excluded] <- corpcor::pseudoinverse(A) %*% as.matrix(Ks_red)
      }
    }
    return(Ns)
  }
  
  # Determine abundances at equilibrium for each trait sample
  for(k in 1:nrow(P_S_E)){
    Ns <- integration.biotic(Ks = P_S_E[k,])
    P_S_E_all_interactions[k,] <- Ns
  }
  dimnames(P_S_E_all_interactions) <- dimnames(P_S_E)
  
  # Generating the function output
  if (out_alpha){
    out = list(P_S_E_all_interactions, alphas)
  } else{
    out = P_S_E_all_interactions
  }
  return(out)
}


