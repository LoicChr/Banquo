##############
#
#     New traitspace objects - data preparation
#
##############
library(tidyverse)
library(plyr)
library(ade4)
library(mclust)

obs.comm <- read.csv("data/community.cover.23sp.csv", header=T, row.names = 1)
KH.data <- read.csv("data/kettleholedata.23sp.csv")

env_p <- unique(KH.data[,c("DAYSUB","enclos", "elev", "QuadID")])
row.names(env_p) <- env_p[,"QuadID"]
env_p <- env_p[row.names(obs.comm),]

obs.comm2 <- as.data.frame(cbind(QuadID=row.names(obs.comm), obs.comm))
obs.comm2[,-1] <- apply(obs.comm2[,-1], 2, as.numeric)
covers <- pivot_longer(obs.comm2, cols = colnames(obs.comm2)[-1], names_to = c("species"), values_to = "Cover")
KH.data <- left_join(KH.data, covers) 
rm(list = c("obs.comm2"))

# Construct traitspace object
KH.data$height <- log(KH.data$height)
KH.data$sla <- log(KH.data$sla)
KH.data$poros <- qlogis(KH.data$poros/100)

# Species trait distribution
pdf_species <- lapply(unique(KH.data$species), function(x){
  Mclust(na.omit(KH.data[KH.data$species==x  , c("poros", "sla", "height")]),warn=TRUE, G = 1)
})
names(pdf_species) <- unique(KH.data$species)

# Check the coverage
# Threshold
Ntr.per.species <- sapply(split(KH.data[,c("sla", "height", "poros")], f = as.factor(KH.data$species)), function(x) min(colSums(!is.na(x))))
sel.species <- names(Ntr.per.species)[Ntr.per.species >=20] # kept species

# Reduction of the dataset to quadrats outside the experimental treatments and that experience flooding
x2 <- env_p$enclos == 0 & env_p$elev < 0.6
env_p <- env_p[x2, ]
obs.comm2 <- obs.comm[row.names(env_p), ]
sel.quadrats <- row.names(obs.comm2)[rowSums(obs.comm2[,sel.species])/rowSums(obs.comm2) > 0.8]

obs.comm <- obs.comm[sel.quadrats, sel.species]
obs.comm <- obs.comm[,colSums(obs.comm) > 0]
obs.comm <- as.matrix(obs.comm/100)
data <- filter(KH.data, species %in% colnames(obs.comm) & QuadID %in% row.names(obs.comm))
env_p <- env_p[row.names(obs.comm),]
pdf_species <- pdf_species[colnames(obs.comm)]

# Average traits
t.avg <- ddply(data, c("species"), summarise,
               maxht  = quantile(height,0.975), sla = median(sla, na.rm = T))

dudi.tr <- dudi.pca(t.avg[,c("maxht", "sla")], nf = 2, scannf = F)
t.avg$Axis1 <- dudi.tr$li[,1]
t.avg$Axis2 <- dudi.tr$li[,2]
save(list = c("t.avg", "pdf_species", "env_p", "obs.comm", "data"), file = "data/data_ready.Rdata")
