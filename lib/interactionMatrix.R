###############################################################################
#                                                                             #
#   Functions to generate the interaction matrix                               #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin       #
#   *contact                                                                  #
###############################################################################

# rho_gen calculates a correlation matrix based on angle parameters (thetas) that determines the Cholesky factorisation of the final correlation matrix.
# It was coded directly from Forrester et Zhang 2020 (Journal of Multivariate Analysis)

# Interaction_matrix calculates the pairwise competition matrix as a function of the intercept, mu, sigma and rhos parameters. Intra specifies the intraspecific interaction coefficients.
rho_gen <- function(thetas){
  N <- (1+sqrt(1+8*length(thetas)))/2
  if (any(thetas  < 0 | thetas > pi)) stop("thetas must be between 0 and pi")
  if (!all.equal(N, as.integer(N)))  stop("number of theta parameters is wrong")
  
  theta_mat <- matrix(0, ncol = N, nrow = N)
  theta_mat[lower.tri(theta_mat, diag = F)] <- thetas
  
  B_mat <- matrix(0, nrow = nrow(theta_mat), ncol = ncol(theta_mat))
  for (i in 1:nrow(B_mat)){
    for (j in 1:i){
      if (j == 1){
        B_mat[i,j] <- cos(theta_mat[i,j])
      } else if (j >= 2 & j <= (i -1)){
        B_mat[i,j] <- cos(theta_mat[i,j])*prod(sin(theta_mat[i,1:(j-1)]))
      } else{
        B_mat[i,j] <- prod(sin(theta_mat[i,1:(j-1)]))
      }
      
    }
  }
  rho_mat <- B_mat %*% t(B_mat)
  return(rho_mat )
}

interaction_matrix <- function(tr, intercept, mu, sigma, rho = NULL, intra = 1){
  if (is.null(rho)){
    rho_mat <- diag(ncol(tr))
  }else{
    rho_mat <- rho_gen(acos(rho))
  }
  if (ncol(tr) == 1){
    inv_sigma_mat <- as.matrix(1/sigma)
  }else{
    inv_sigma_mat <- solve(diag(sigma) %*% rho_mat %*% diag(sigma))  
  }

  aij = matrix(NA, nrow(tr), nrow(tr))
  for (i in 1:nrow(tr)){
    for (j in 1:nrow(tr)){
      dt <- as.matrix(tr[i,]- tr[j,])
      aij[i,j] = exp(-0.5* t(dt -mu) %*% inv_sigma_mat %*% (dt -mu))
    }
  }
  aij = intercept*aij
  diag(aij) <- intra
  
  return(aij)
}