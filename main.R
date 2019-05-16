###############################################################################
#                                                                             #
#   Main script to perform the analysis                                       #
#   From the article: Integrating traits, competition and abiotic filtering   #
#                 improves biodiversity predictions                           #
#   Authors: Loïc Chalmandrier*, Daniel B. Stouffer, Daniel C. Laughlin       #
#   *contact                                                                  #
###############################################################################

rm(list=ls())

# Load data
load("data.Rdata")
#KH.data: 
# comm: species cover in each site.
# env: number of flooded days (sub) in each site
# KH.data: data.frame with the characteristics of each sampled plant individual.
### Columns: the site name (QuadID), the species names (species), the number of flood days (sub), the root aerenchyma (poros), sla, root dry mass content (rdmc) and the height. 

# libraries
library(mvtnorm)
library(mclust)
library(topicmodels)
library(nloptr)
library(plyr)
library(ade4)
print("libraries loaded")

# Useful functions
source("interactionMatrix.R")
source("traitspace_and_banquo.R")
source("nlptr_comput.R")

## Preparation of data
t.avg <- ddply(KH.data, c("species"), summarise,
               maxht  = quantile(height,0.975), sla = mean(sla, na.rm = T))
t.avg$maxht <- log(t.avg$maxht)
t.avg$slaRes  <- residuals(lm(sla ~ maxht, data = t.avg)) 

dudi.tr <- dudi.pca(t.avg[,c("maxht", "sla")], nf = 2, scannf = F)
t.avg$Axis1 <- dudi.tr$li[,1]
t.avg$Axis2 <- dudi.tr$li[,2]

# Traitspace
data <- KH.data
data$poros <- log(data$poros)
multi.mod <- lm(cbind(poros,rdmc) ~ sub, data = data)

pdf_species <- lapply(levels(data$species), function(x){
  Mclust(na.omit(data[data$species==x, c("poros","rdmc")]), warn=TRUE)
})
names(pdf_species) <- levels(data$species)
P_S_E_tra <- traitspace(multi.mod, env, pdf_species, N = 20, avg.out =T)

### Height model
banquo_out_H <- ML_interactions(P_S_E_tra, trait = t.avg$maxht, comm, lb = c(-5, -0.9, -2), ub = c(-0.5, 0.9, 0.75), 
                                opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval = 3e6, xtol_rel = 10e-10, population = 10))

### SLA model
banquo_out_S <- ML_interactions(P_S_E_tra, trait = t.avg$sla, comm, lb = c(-5, -0.9, -2), ub = c(-0.5, 0.9, 0.75), 
                                opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval = 3e6, xtol_rel = 10e-10, population = 1000))

### HeightSLA models
banquo_out_HS <- ML_interactions_bi(P_S_E_tra, trait1 = t.avg$maxht,trait2 = t.avg$slaRes, comm, det_lim = exp(-1.5),
                                 lb = c(-5, -0.9,-1.5,-0.9, -1.5,-0.9), ub = c(-1, 0.9, 2, 0.9,2,0.9),
                                 opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval = 3e6, xtol_rel = 10e-10, population = 2000))

banquo_out_HS2 <- ML_interactions_bi2(P_S_E_tr=P_S_E_tra, trait1 = t.avg$Axis1, trait2 = t.avg$Axis2, obs.comm = comm, det_lim = exp(-1.5),
                                    lb = c(-5, 0,-1.5,0, -1.5,-0.9), ub = c(-1, 0.8, 2, 2*pi,2,0.9),
                                    opts = list(algorithm = "NLOPT_GN_CRS2_LM", maxeval = 3e6, xtol_rel = 10e-10, population = 2000))

