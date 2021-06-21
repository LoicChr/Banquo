###############################################################################
#                                                                             #
#   Function to generate the interaction matrix                               #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier, Daniel B. Stouffer, Daniel C. Laughlin        #
#   *contact                                                                  #
###############################################################################

interaction.matrix <- function(trait, intercept=1, mu=0, sigma=0.3, intra=1, std = T){
  # Arguments
  ## trait: he single trait across species.
  ## intercept: the multiplying coefficients on interaction coefficients.
  ## mu, sigma: the peak point and the width of the bell-shaped curve between trait difference and pairwise interaction coefficients
  ## intra: the additional coefficient for intraspecific coefficients.
  
  A <- matrix(0, length(trait),length(trait))
  for(ii in 1:nrow(A)){
    for(jj in 1:nrow(A)){
      A[ii,jj] <- -0.5*((trait[ii] - trait[jj] - mu)/sigma)^2
    }
  }
  if (std){
    A <- intercept * exp(A)/(sqrt(2*pi)*sigma)
  }else{
    A <- intercept * exp(A)
  }
 
  diag(A) <- intra 
  return(A)
}

interaction.matrix_bi <- function(trait1, trait2, intercept=1, mu1=0, sigma1=0.3, mu2 = 0, sigma2 = 0.3, rho = 0, intra= 1, std = T){
  # Arguments
  ## trait1, trait2: the two traits across species.
  ## intercept: the multiplying coefficients on interaction coefficients.
  ## mu1, mu2, sigma1, sigma2, rho: the peak point, the width of the bell-shaped curve (in two dimensions) and the interaction coefficients among the two dimensions between trait difference and pairwise interaction coefficients
  ## intra: the additional coefficient for intraspecific coefficients.
  
   if (length(trait1) != length(trait2)) stop("trait1 and trait2 do not match")
  A <- matrix(0, length(trait1), length(trait1))
  for(ii in 1:nrow(A)){
    for(jj in 1:nrow(A)){
      A[ii,jj] <- - (1/(2*(1-rho^2)))*(((trait1[ii] - trait1[jj] - mu1)/sigma1)^2 
                  + ((trait2[ii] - trait2[jj] - mu2)/sigma2)^2 
                  - 2*rho*((trait1[ii] - trait1[jj] - mu1)*(trait2[ii] - trait2[jj] - mu2)/(sigma2*sigma1))) 
    }
  }
  # take the exponential
  if (std){
    A <- intercept * exp(A)/(2*pi*sigma1*sigma2*sqrt(1-rho^2))
  }else{
    A <- intercept * exp(A)
  }
  diag(A) <- intra
  return(A)
}
