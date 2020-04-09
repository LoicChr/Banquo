load("ozab_data.Rdata")

# Pick a community
pred <- P_S_E_tra[50,]
obs <- obs[50,]

logit <- function(x) log(x/(1-x))
expit <- function(x) 1/(1+exp(-x))
# Param
phi <- 0.6
a_ab <- 0.5
b_ab <- 0.1
sigma <- 2
delta <- 0.8
a_pa <- 0.5
b_pa <- 0.1

ozab <- function(obs, pred){
  # Probability of observing the presence of species
  
  # If phi depends on predicted species cover
  # phi <- expit(a_pa*pred[x]+b_pa)
  prob_pres <- dbinom(obs > 0, 1, phi)
  
  # If species is present
  x <- obs > 0
  mean_cover <- expit(a_ab*pred[x]+b_ab) 
        #Note I thought how having a error term determined by sigma here. 
        ##But I don't know how to integrate it in that equation.
  
  alpha = mean_cover * delta #Conversion in alpha/beta shape parameters
  beta = (1-mean_cover) * delta
  j <- match(obs[x], cover_class$cover_class)
  
  prob_cat <- pbeta(cover_class$maxCov[j], alpha, beta) - pbeta(cover_class$minCov[j], alpha, beta)
  prob_ab <- prob_pres
  prob_ab[x] <- prob_ab[x]*prob_cat
  
  logLik <- sum(log(prob_ab))
  return(logLik)
}
