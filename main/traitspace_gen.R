##############
#
#     New traitspace objects - linear models
#
##############

library(mclust)
library(plyr)
library(ade4)

source("main/data_prep.R")
source("lib/traitspace_and_banquo.R")

traits = c('sla', "height", "poros")
x1 <- data$enclos == 0 & data$elev < 0.6
x2 <- env_p$enclos == 0 & env_p$elev < 0.6

# Threshold
tr.eff <- sapply(split(KH.data[,c("sla", "height", "poros")], f = as.factor(KH.data$species)), function(x) min(colSums(!is.na(x))))

# Species trait distribution
pdf_species <- lapply(levels(data$species), function(x){
  Mclust(na.omit(data[data$species==x  , traits]),warn=TRUE, G = 1)
})
names(pdf_species) <- levels(data$species)

# Reduction of the dataset
env_p <- env_p[x2, ]
obs.comm2 <- obs.comm[row.names(env_p),names(tr.eff)[tr.eff >=20] ]
obs.comm2 <- obs.comm2[rowSums(obs.comm2)/rowSums(obs.comm[row.names(env_p),]) >=0.8,]
obs.comm2 <- obs.comm2[,colSums(obs.comm2) > 0]
data <- data[data$species %in% colnames(obs.comm2) & data$QuadID %in% row.names(obs.comm2),]
data$species <- droplevels(data$species)

pdf_species <- pdf_species[levels(data$species)]

## Trait model
multi.mod <- lm(cbind("poros", "sla", "height") ~ DAYSUB, data = na.omit(data), weights = Cover)
P_S_E_tr <- traitspace(multi.mod, env_p, pdf_species, N = 500, avg.out = T)
obs.comm <- as.matrix(obs.comm[row.names(P_S_E_tr), colnames(P_S_E_tr)])/100
P_S_E_tr <- as.matrix(P_S_E_tr)
save(list = c("P_S_E_tr", 'env_p', "obs.comm", 'pdf_species',"multi.mod"), file = "data/traitspace_object.Rdata")
