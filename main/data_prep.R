##############
#
#     New traitspace objects - data preparation
#
##############
library(tidyverse)
library(plyr)
library(ade4)
obs.comm <- read.csv("data/community.cover.23sp.csv", header=T, row.names = 1)
KH.data <- read.csv("data/kettleholedata.23sp.csv")

env_p <- unique(KH.data[,c("DAYSUB","enclos", "elev", "QuadID")])
row.names(env_p) <- env_p[,"QuadID"]
env_p <- env_p[row.names(obs.comm),]

KH.data$poros[KH.data$poros == 23.49123724] <- NA
KH.data$sla[KH.data$sla == 15.88267954] <- NA
KH.data$rdmc[KH.data$rdmc == 0.260117717] <- NA
KH.data$srl[KH.data$srl == 75.37941805] <- NA

KH.data$species <- as.character(KH.data$species)

obs.comm2 <- as.data.frame(cbind(QuadID=row.names(obs.comm), obs.comm))
obs.comm2[,-1] <- apply(obs.comm2[,-1], 2, as.numeric)
covers <- pivot_longer(obs.comm2, cols = colnames(obs.comm2)[-1], names_to = c("species"), values_to = "Cover")
KH.data <- left_join(KH.data, covers) 
rm(list = c("obs.comm2"))

# Construct traitspace object
data <- na.omit(KH.data[,c("poros", "height","sla", "DAYSUB","QuadID", "elev","enclos", "species", "Cover")])
data$height <- log(data$height)
data$sla <- log(data$sla)
data$poros <- qlogis(data$poros/100)


t.avg <- ddply(KH.data, c("species"), summarise,
               ht = mean(height),maxht  = quantile(height,0.975), sla = mean(sla, na.rm = T), poros = mean(poros, na.rm = T), rdmc = mean(rdmc, na.rm=T))
t.avg$maxht <- log(t.avg$maxht)
t.avg$sla <- log(t.avg$sla)


dudi.tr <- dudi.pca(t.avg[,c("maxht", "sla")], nf = 2, scannf = F)
t.avg$Axis1 <- dudi.tr$li[,1]
t.avg$Axis2 <- dudi.tr$li[,2]
save(t.avg, file = "data/tra2009/Travg.rdata")
