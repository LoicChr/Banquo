#Data prep
obs <- read.csv("data/kettlehole_dataset/community.cover.15sp.csv", header=T, row.names = 1)
load("data/tra2003/traitspace.Rdata")
cover_class <- sort(unique(unlist(obs[obs > 0])))
min_cut <- 0.5*(cover_class + c(0,cover_class[-length(cover_class)]))/100
max_cut <- 0.5*(cover_class + c(cover_class[-1], 100))/100
max_cut[which.max(max_cut)] <- 1
min_cut[which.min(min_cut)] <- 0

obs <- obs[row.names(P_S_E_tra), colnames(P_S_E_tra)]
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
  j <- match(obs[x], cover_class)
  
  prob_cat <- pbeta(max_cut[j], alpha, beta) - pbeta(min_cut[j], alpha, beta)
  prob_ab <- prob_pres
  prob_ab[x] <- prob_ab[x]*prob_cat
  
  logLik <- sum(log(prob_ab))
  return(logLik)
}
