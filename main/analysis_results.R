 library(tidyverse)
 library(BayesianTools)
# library(ggplot2)
 library(pROC)
 library(ggpubr)
 library(mgcv)
# library(rcartocolor)
# library(colorBlindness)

source("lib/abgFunctions.R")
chains.to.analyze <- list.files("results", full.names = T, recursive = T, pattern = "chain")

species <- c(AGRcap = "Agrostis capillaris", AMPflu = "Amphibromus fluitans", ANTodo = "Anthoxanthum odoratum", CARgau = "Carex gaudichaudiana", ELAacu = "Eleocharis acuta", ELApus =  "Eleocharis pusilla", EPIang = "Epilobium angustum", GALper = "Galium perpusillum", JUNart = "Juncus articulatus", LACstr = "Lachnagrostis striata", LAClya = "Lachnagrostis lyallii", LILrut = "Lilaeopsis ruthiana", LOBper = "Lobelia perpusilla", PARcan = "Parahebe canescens", PILoff = "Pilosella officinarum", PILpil = "Pilosella piloselloides")


###
out_results<- list()
for (ik in 1:length(chains.to.analyze)){
  chain = chains.to.analyze[ik]
  b=load(chain)
  source("./lib/traitspace_and_banquo.R")
  ### Isolate posterior
  nstep <- nrow(getSample(out, start = 0, parametersOnly = F, thin = 0))/(4*3)
  thres <- nstep - 2e5
  dis <- getSample(out, parametersOnly = F, start = thres, thin = 50)
  
  pars <- apply(dis, 2, median)
  
  ### General performance stats
  dic.val <- DIC(out, start = thres)$DIC
  conv.val <- gelmanDiagnostics(out, start = thres)$mpsrf
  
  ### Trait differences
  if (grepl("Height|SLA", chain)){
    diff <- outer(trait, trait, '-')
  } else if (grepl("2tr", chain)){
    diff <- list( outer( t.avg$maxht,  t.avg$maxht, '-'),  outer( t.avg$sla,  t.avg$sla, '-'))
  }else{
    diff = NULL
  }

  ### Compute spxp matrix, likelihood matrix at the median posterior
  pred.med = likelihoodAb(pars, spxp = T)
  if (traits_biotic == "none"){
    spxp.pred =  pred.med[[2]]
    int = diag(ncol(spxp.pred))
  }else{
    spxp.pred =  pred.med[[2]][[1]]
    int <- pred.med[[2]][[2]]
  }
  ### Distribution of diversities
  div_distri <- do.call(rbind, lapply(sample(1:nrow(dis), 500, replace = T), function(i){
    if (traits_biotic == "none"){
      spxp.pred = likelihoodAb(dis[i,], spxp = T)[[2]]
    }else{
      spxp.pred = likelihoodAb(dis[i,], spxp = T)[[2]][[1]]
    }
  
    spxp.simul <- matrix(NA, nrow = nrow(spxp.pred), ncol = ncol(spxp.pred))
    for (i in 1:nrow(spxp.pred)){
      for (j in 1:ncol(spxp.pred)){
        alpha = spxp.pred[i,j] * dis[i,"phi"]
        beta = (1 - spxp.pred[i,j]) * dis[i,"phi"]
        spxp.simul[i,j] <-  rbeta(1, alpha, beta)
      }
    }
    c(unlist(abgDecompQ(spxp.pred, q= 2)[2:3]), unlist(abgDecompQ(spxp.simul, q= 2)[2:3]))
    
  }))
  colnames(div_distri) <- c("Beta", "mAlpha", "Beta_simul", "mAlpha_simul")
  div_distri <- as.data.frame(div_distri)
  roc1 <- roc(as.numeric(obs.comm > 0),as.numeric(  spxp.pred ))
  
  out_results[[ik]] <- list(DIC =  dic.val , conv = conv.val, diff_tr = diff, alpha = int, alpha_beta = div_distri, ll = median(dis[,"Llikelihood"]), spxp = spxp.pred, roc = roc1, Nobs = prod(dim(obs.comm)), ll_025 = quantile(dis[,"Llikelihood"], 0.025), ll_975 = quantile(dis[,"Llikelihood"], 0.975))
  rm(out)
  print(ik)
}
names(out_results) = chains.to.analyze

dat <- data.frame(names(out_results), conv = map_dbl(out_results, 'conv'), DIC = map_dbl(out_results, "DIC") )
lls <- map_dbl(out_results, 'll')
Nobs <- map_dbl(out_results, 'Nobs')
ll_H0 <- lls[which(grepl("abio_FALSE/none", names(out_results)))]
dat$pseudoR2 <- (1-exp(2/Nobs*(ll_H0-lls)))/(1-exp(2/Nobs*ll_H0))
dat$pseudoR2_025 <- (1-exp(2/Nobs*(ll_H0-map_dbl(out_results, 'll_025'))))/(1-exp(2/Nobs*ll_H0))
dat$pseudoR2_975 <- (1-exp(2/Nobs*(ll_H0-map_dbl(out_results, 'll_975'))))/(1-exp(2/Nobs*ll_H0))
dat$AUC  <- sapply(map(out_results, 'roc'), auc)

dat2 <- dat
dat2 <- separate(as.data.frame(dat2), col = 1, into = c("folder", "abio", "bio", "chain"),sep= "/")
row.names(dat2) <- NULL
dat2 <- select(dat2, -folder, -chain)

dat2$conv <- round(dat2$conv, 3)
dat2$DIC <- round(dat2$DIC, 1)
dat2[, c("pseudoR2_025","pseudoR2_975")] <- replace(dat2[, c("pseudoR2_025","pseudoR2_975")]  , dat2[, c("pseudoR2_025","pseudoR2_975")] < 0, 0)
dat2[, c("pseudoR2_025","pseudoR2_975")] <- round(dat2[, c("pseudoR2_025","pseudoR2_975")] , 3)
dat2$pseudo_char <- paste0('[', round(dat2$pseudoR2_025,3), ', ', round(dat2$pseudoR2_975,3), "]")
dat2$AUC <- round(dat2$AUC, 3)
dat2 <- dat2[order(dat2$pseudoR2, decreasing = T),]

dat2 <- select(dat2, -pseudoR2, -pseudoR2_025, -pseudoR2_975)
dat2$Assembly_model <- 'Null model'
dat2$Assembly_model[dat2$abio == "abio_FALSE" & dat2$bio != "none"] <- "Biotic model"
dat2$Assembly_model[dat2$abio == "abio_TRUE" & dat2$bio != "none"] <- "Abiotic & Biotic model"
dat2$Assembly_model[dat2$abio == "abio_TRUE" & dat2$bio == "none"] <- "Abiotic model"
dat2$bio <- factor(dat2$bio, levels = c("SLA","2tr_rho","2tr","Height","none","NoTr"),  labels = c("SLA","Height + SLA (with int.)","Height + SLA (no int.)","Height","No interactions","No traits"))

dat2 <- dat2[,c("Assembly_model", "bio", "conv", "DIC", "pseudo_char", "AUC")]
write.csv2(dat2, file = "results/tab.csv", row.names = F, quote = F)

###### Figure 1 ####
load(chains.to.analyze[1])
data <- multi.mod$model
colnames(data)
data <- cbind(data[, "cbind(poros, sla, height)" ], data[,c("DAYSUB", "(weights)")])
colnames(data)[5] <- "Cover"
traits <- c("sla", "height", "poros")
mod <-lapply(traits, function(tr){
  reg <- lm(as.formula(paste0(tr, "~ DAYSUB")), data = data, weights = Cover)
  pred = predict(reg, interval = "confidence")
  colnames(pred) <- paste(tr, colnames(pred), sep = "_")
  pred
})

colnames(data)[1:3] <- paste(colnames(data)[1:3], "resp",sep = "_")
data_ext <- cbind(data, do.call(cbind, mod))
data_ext[,grepl("poros",colnames(data_ext))] <- 100*apply(data_ext[,grepl("poros",colnames(data_ext))],2, plogis)
data_ext[,grepl("sla|height",colnames(data_ext))] <- exp(data_ext[,grepl("sla|height",colnames(data_ext))])
data_long <- pivot_longer(data_ext, cols = colnames(data_ext)[grepl('sla|height|poros',colnames(data_ext))], names_to = c("trait",".value"),names_sep = "_")

lbs = setNames(c("'SLA ('*mm^2*'/mg)'", "Height~(cm)", "'Root porosity (%)'"), c("sla", "height", "poros"))

Fig1 <- ggplot(data_long, aes(x = DAYSUB, y = resp)) + 
  theme_bw() +
  geom_point(aes(size = Cover), shape = 1, color = "black") + 
  geom_line( aes(y = fit), size = 1, color = "#CC6677") +
  theme(axis.text = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(size = "Cover (%)", x = "Days submerged per year", y = "Trait value") +
  geom_ribbon( aes(ymin = lwr, ymax = upr, fill = "#CC6677"), alpha = 0.5) +
  guides(fill= FALSE) +
  facet_wrap(~trait, scales = "free",labeller=as_labeller(lbs, label_parsed))
ggsave("results/Fig1.jpg", Fig1, width = 8, height = 3)
######

###### Figure 2 ###############
sel.models <- which(grepl("abio_TRUE", names(out_results)) & !grepl("none", names(out_results)))
max_int <- max(sapply(map(out_results[sel.models],'alpha'), function(x){
  diag(x) <- NA
  max(x, na.rm =T)
}))
interaction_plots <- lapply(sel.models,function(i){
  
  if (!grepl("2tr", names(out_results)[i])){
    if (!grepl("NoTr",  names(out_results)[i])){
      diff = out_results[[i]]$diff
      title_name <- trait_name <- ifelse(grepl("SLA",  names(out_results)[i]), "SLA", "Height")
    }else{
      diff = out_results[[sample(which(grepl("Height", names(out_results))),1)]]$diff
      trait_name <- "Height"
      title_name <- "No traits"
    }
    mat = out_results[[i]]$alpha
    diag(mat) <- NA
    dat.int = tibble(diff = as.numeric(diff), int = as.numeric(mat))
    gp <- ggplot(dat.int, aes(x = diff, y = int, fill = int))+ scale_fill_viridis_c(option = "viridis", limits = c(0, max_int))   +
      theme(panel.grid.major = element_blank(), plot.margin = unit(rep(1.5,4), "lines"), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      geom_point(shape = 21, size = 2.5) + ggtitle(title_name) + ylim(c(0,max_int)) +
      labs(y='Pairwise interaction', x = bquote(paste(.(trait_name),' log-ratio (log(',t[i],'/',t[j],')')),fill = "Interaction\nstrength") + 
      geom_vline(xintercept = 0, linetype = 2, color = "black")
    
  }else {
    title_name <- ifelse(grepl("rho_FALSE",  names(out_results)[i]), "Height & SLA (no int.)", "Height & SLA (with int.)")
    
    diff1 = out_results[[i]]$diff[[1]]
    diff2 = out_results[[i]]$diff[[2]]
    mat = out_results[[i]]$alpha
    diag(mat) <- NA
    dat.int = tibble(diff1 = as.numeric(diff1), diff2 = as.numeric(diff2), int = as.numeric(mat)) %>% drop_na()
    
    gp <- ggplot(dat.int, aes(x = diff1, y = diff2, fill = int))  +
      theme(panel.grid.major = element_blank(), plot.margin = unit(rep(1.5,4), "lines"), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
      geom_point(shape = 21, size = 2.5) +  geom_vline(xintercept = 0, linetype = 2, color = "black") + geom_hline(yintercept = 0, linetype = 2, color = "black")+
      labs(y=expression(paste('SLA log-ratio (log(',t[i],'/',t[j],')')), x = expression(paste('Height log-ratio (log(',t[i],'/',t[j],')')), fill = "Interaction\nstrength") + ggtitle(title_name ) + scale_fill_viridis_c(option = "viridis", limits = c(0, max_int), direction = 1) +
      annotate("point", x = as.numeric(filter(dat.int, int == max(int))[1,"diff1"]), y = as.numeric(filter(dat.int, int == max(int))[1,"diff2"]), fill = "red", shape=23, size = 4.5)
  }
  return(gp)
})
names(interaction_plots) <- map_chr(strsplit(names(out_results)[sel.models],"/"),3)
Fig2 <- ggarrange(plotlist = interaction_plots[c("NoTr", "Height", "SLA", "2tr", "2tr_rho")], nrow = 3, ncol = 2, common.legend = T, legend = 'right', labels = LETTERS[1:length(interaction_plots)]) 
ggsave("results/Fig2.jpg", Fig2, width = 9.26, height = 10)
#######

####### Figure 3 ################
b=load(chains.to.analyze[1])
id.abio <- which(grepl("abio_TRUE", names(out_results)) & grepl("none", names(out_results)))
id.best <- which(grepl("abio_TRUE", names(out_results)) & grepl("SLA", names(out_results)))
species <- species[colnames(obs.comm)]

format_data <- function(out_results, id.abio, id.best){
  dat_abio <- as.data.frame.table(out_results[[id.abio]]$spxp)
  dat_best <- as.data.frame.table(out_results[[id.best]]$spxp)
  dat_obs <-  as.data.frame.table(obs.comm)
  
  dat.dummy <- merge(dat_abio, dat_best, by = c("Var1",   "Var2"), suffixes = c(".abio", ".best"))
  dat.tot <- merge(dat.dummy, dat_obs, by = c("Var1",   "Var2"))

  colnames(dat.tot)[c(1,2,5)] <- c("QuadID", "Species",'Freq.obs')
  colnames(dat.tot) <- gsub("Freq", "Cover", colnames(dat.tot))
  dat.tot2 <- merge(dat.tot, data.frame(QuadID = row.names(env_p), DAYSUB=env_p[,"DAYSUB"]), by = "QuadID")
  dat.tot3 <- pivot_longer(dat.tot2, cols = starts_with('Cover'), names_to = 'Model', values_to = 'Cover') %>% 
    mutate(Species = forcats::fct_recode(Species,!!! setNames(names(species), species))) %>%
    mutate(Model = factor(Model, levels = c("Cover.obs", "Cover.abio","Cover.best"), labels = c("Observed", "Abiotic model","Abiotic + Biotic model")))
  
  rm(list = c('dat_abio', 'dat_best', 'dat_obs', 'dat.tot', 'dat.dummy', 'dat.tot2'))
  return(dat.tot3)
}

dat <- out_results %>% format_data(id.abio = id.abio, id.best = id.best) 
Fig3 <- ggplot(dat, aes(x = DAYSUB, y = Cover*100)) +
  geom_point(data = filter(dat,Model== 'Observed'), fill = 'black') + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
  geom_smooth(mapping = aes(color = Model), se = F) + 
  scale_color_manual(values=c("#000000", "#CC6677", "#117733")) +
  facet_wrap(~Species, scales = 'free') + scale_y_continuous(limits =c(0,80),breaks = c(0,5,20,50,80)) +  coord_trans(y = "sqrt") +
  labs(x = "Days submerged per year", y = "Species cover (%)") + guides(linetype = FALSE) 
ggsave("results/Fig3.jpg", Fig3, width = 8, height = 8)
################

#### Supplementary Figure 4 ############
load(chains.to.analyze[1])

mat_pred <- sapply(1:ncol(obs.comm), function(i){
  resp <- obs.comm[,i]
  resp[resp == 0] <- min(obs.comm[obs.comm > 0])
  reg <- gam(resp ~ s(env_p$DAYSUB), family = "betar")
  y = fitted(reg)
})
source("lib/abgFunctions.R")
ab.obs.med = abgDecompQ(mat_pred, q = 2)
ab.obs = abgDecompQ(obs.comm, q = 2)

format_data2 <- function(x, pos){
  x %>% map(pos) %>% bind_rows(.id = "Chain") %>% 
    separate(Chain, sep = "/", remove = F, into = c("folder","abio", "bio_traits", "chain_number")) %>%
    mutate(Assembly_model = fct_cross(abio, bio_traits, sep='_'),
          Assembly_model = fct_collapse(Assembly_model, 
                                           "Null model" = "abio_FALSE_none",
                                           "Abiotic model" ="abio_TRUE_none",
                                           'Biotic models' = str_subset(levels(Assembly_model),"abio_FALSE_(?!none)"),
                                           'Abiotic & Biotic models' = str_subset(levels(Assembly_model),"abio_TRUE_(?!none)")),
            Assembly_model = factor(Assembly_model, levels = c("Null model", "Abiotic model",'Biotic models', 'Abiotic & Biotic models')),
          bio_traits = factor(bio_traits, levels = c('none', 'NoTr','Height', 'SLA', '2tr', '2tr_rho')),
           bio_traits = fct_recode(bio_traits, "No interactions" = "none", 
                                   "No traits"="NoTr",
                                    "Height & SLA (no int.)" = "2tr",
                                     "Height & SLA (with int.)" = "2tr_rho"),
           bio_traits = factor(bio_traits, c("No interactions","No traits","Height", "SLA", "Height & SLA (no int.)","Height & SLA (with int.)"))) %>%
          select(-folder, -abio, -chain_number, -Chain)
    
}

dat2 <- out_results %>% format_data2('alpha_beta') %>% pivot_longer(cols = c("Beta", 'mAlpha', 'Beta_simul', 'mAlpha_simul'), names_to = "diversity_name", values_to = "diversity")
dat2$diversity_name <- factor(dat2$diversity_name, levels = c("Beta", "Beta_simul", "mAlpha_simul", "mAlpha"),
                              labels = c(expression(paste(beta,'-diversity (fixed)')),
                                         expression(paste(beta,'-diversity (stochastic)')),
                                         expression(paste(alpha,'-diversity (stochastic)')),
                                         expression(paste(alpha,'-diversity (fixed)'))))
ref_tab = as_tibble(data.frame(diversity_name = levels(dat2$diversity_name), ref_line = c(ab.obs.med[["Beta"]], ab.obs[['Beta']],ab.obs[['mAlpha']],ab.obs.med[["mAlpha"]])))

FigS4 <- ggplot(dat2, aes_string(y= 'diversity',x = 'Assembly_model', fill = 'bio_traits')) + 
  theme_bw() +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(data= ref_tab,aes(yintercept = ref_line), linetype="dashed", size = 1) +
    facet_wrap( ~diversity_name, drop = T, scale = 'free_y', labeller = label_parsed) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text = element_text(colour = "black"),axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(fill = "Competition traits") +scale_fill_manual(values=c("#FFFFFF",carto_pal(n=7, name="Safe")[-c(5,7)]))
ggsave("results/FigS4.jpg", FigS4, width = 8, height = 6.5)

#############

## Supplementary Figure 5 ######
p1 <- ggplot(dat, aes(x = Cover*100, col = Model)) + scale_color_manual(values = cols[-3]) +
  theme(legend.position = "top", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + labs(x = 'Species cover (%)') + 
  scale_x_continuous(trans = 'sqrt', breaks = c(0.5,2, 5,10,20,40,80), labels = as.character(c(0.5,2, 5,10,20,40,80)))+
  geom_density(size = 1, adjust = 2.5 )

dat3 <- dat%>% select(-DAYSUB) %>% pivot_wider(values_from = "Cover", names_from = "Model") %>% 
  mutate(Observed = factor(Observed*100), `Abiotic model` = 100*`Abiotic model`, `Abiotic + Biotic model` = 100*`Abiotic + Biotic model`)

p2 <- ggplot(dat3, aes(y = `Abiotic + Biotic model`, x = Observed)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + 
  labs(x = 'Observed cover class (%)', y = "Percent cover (Abiotic & biotic model)") + ylim(c(0, max(dat.tot4$`Abiotic model`, dat.tot4$`Abiotic + Biotic model`))) +
  geom_boxplot(varwidth = FALSE, fill = cols[4], alpha = 0.5)

p3 <- ggplot(dat3, aes(y = `Abiotic model`, x = Observed)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + 
  labs(x = 'Observed cover class (%)', y = "Percent cover (Abiotic model)") + ylim(c(0, max(dat.tot4$`Abiotic model`, dat.tot4$`Abiotic + Biotic model`))) +
  geom_boxplot(varwidth = FALSE, fill = cols[2], alpha = 0.5) 
p4 <- ggplot(dat.tot4, aes(x = `Abiotic model`, y = `Abiotic + Biotic model`)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black"))  + 
  geom_point(shape = 19, alpha = 0.3) + labs(x = 'Percent cover (Abiotic model)', y = "Percent cover (Abiotic & biotic model)")
FigS5 <- ggarrange(plotlist = list(p3, NULL, p2, NULL, NULL, NULL, p4,NULL, p1), heights = c(1,0.1,1), widths = c(1,0.05,1), nrow = 3, ncol = 3, labels = c('A',NA,'B', rep(NA,3), 'C', NA, 'D'), common.legend = T) 
ggsave("results/FigS5.jpg", FigS5, width = 9, height = 8)
###########


## Supplementary Figure 6 ######
rocs <- map(out_results, 'roc')
rocs.dat <- bind_rows(lapply(rocs, function(x) data.frame(spec= 1-x$specificities, sens = x$sensitivities)), .id = "Chain") %>% 
  separate(Chain, sep = "/", remove = F, into = c("folder","abio", "bio_traits", "chain_number")) %>%
  mutate(Assembly_model = fct_cross(abio, bio_traits, sep='_'),
         Assembly_model = fct_collapse(Assembly_model, 
                                       "Null model" = "abio_FALSE_none",
                                       "Abiotic model" ="abio_TRUE_none",
                                       'Biotic models' = str_subset(levels(Assembly_model),"abio_FALSE_(?!none)"),
                                       'Abiotic & Biotic models' = str_subset(levels(Assembly_model),"abio_TRUE_(?!none)")),
         Assembly_model = factor(Assembly_model, levels = c("Null model", "Abiotic model",'Biotic models', 'Abiotic & Biotic models')))

cols = c("#888888","#CC6677","#332288","#117733") 
FigS6 <- ggplot(rocs.dat, aes(x = spec, y = sens, group = Chain, color = Assembly_model)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
  geom_line() + labs(x = "False positive rate", y = 'True positive rate') + 
  scale_colour_manual(name = "grp",values = cols) 
ggsave("results/FigS6.jpg", FigS6, width = 6, height = 6)


# Posterior distribution figures
library(extraDistr)
names_params <- c("a","b", expression(phi), 'log(C)', expression(mu), expression(paste('log(',sigma,')')), expression(mu[1]), expression(paste('log(',sigma[1],')')), expression(mu[2]), expression(paste('log(',sigma[2],')')), expression(rho))
all.params <- c("a2", "c2", "phi", "intercept","mu","sigma","mu1","sigma1", "mu2", "sigma2", "rho")
names(all.params) <- names_params
names_models <- gsub("results/|/chains.Rdata", "", chains.to.analyze)
names_models <- factor(names_models , levels = c('abio_FALSE/none','abio_TRUE/none','abio_FALSE/NoTr','abio_FALSE/Height','abio_FALSE/SLA', 'abio_FALSE/2tr', 'abio_FALSE/2tr_rho',
                                                            'abio_TRUE/NoTr','abio_TRUE/Height','abio_TRUE/SLA', 'abio_TRUE/2tr', 'abio_TRUE/2tr_rho'), labels = paste("Model", 1:12, sep = "_"))
for (ip in 1:length(chains.to.analyze)){
  chain = chains.to.analyze[ip]
  b=load(chain)
  ### Isolate posterior
  nstep <- nrow(getSample(out, start = 0, parametersOnly = F, thin = 0))/(4*3)
  thres <- nstep - 2e5
  post.dis <- as.data.frame(getSample(out, parametersOnly = T, start = thres, thin = 1000))
  colnames(post.dis) <- list_params
  post.dis$distri <- 'Posterior'
  prior.dis <- data.frame(sampler(nrow(post.dis)*100))
  colnames(prior.dis) <- list_params
  for (i in 1:ncol(prior.dis)){
    prior.dis[prior.dis[,i] < bounds[i,1] |     prior.dis[,i] > bounds[i,2],i] <- NA
  }
  prior.dis[prior.dis$phi > 10,"phi"] <- NA
  prior.dis$distri <- 'Prior'
  bounds$param = factor(row.names(bounds), levels=all.params[all.params %in% list_params])
  bounds<-bounds %>% mutate(param = fct_recode(param, !!!all.params[all.params %in% all_of(list_params)]))
  distris <- rbind(prior.dis, post.dis)
  distris[,grep('intercept|sigma', colnames(distris))] <- log(distris[,grep('intercept|sigma', colnames(distris))])
  bounds[grep('intercept|sigma', row.names(bounds)),c('lower', 'upper')] <- log(  bounds[grep('intercept|sigma', row.names(bounds)),c('lower', 'upper')] )
  distris_long <- distris %>% pivot_longer(cols = list_params, names_to = "param") %>% 
    mutate(param = fct_recode(param, !!!all.params[all.params %in% all_of(list_params)])) %>% 
    mutate(param = factor(param, levels = names(all.params)[all.params %in% list_params]))
  gp1 = ggplot(distris_long, aes(x = value, col = distri)) + geom_density() + 
    geom_vline(data = bounds, aes(xintercept = lower))+ theme(legend.position="top") +
    geom_vline(data = bounds, aes(xintercept = upper)) +
    facet_wrap(~param, scales = "free", labeller = label_parsed) +labs(x = "Parameter value", y = "Density", col = "Distribution")
  if (nrow(bounds) <=3){
    height  = 4.2
  } else if(nrow(bounds) > 3 & nrow(bounds) <= 6){
    height = 5.3
  }else{
    height = 8.42
  }
  
  ggsave(gp1, file = paste0("results/post_models/post_",as.character(names_models[ip]) ,".jpg"),width = 8.16,height = height,units = 'in')
}

