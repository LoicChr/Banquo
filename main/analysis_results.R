################################################################
#                                                              #
#                   Figures and tables                         #
#                                                              #
################################################################

# Libraries
library(tidyverse)
library(BayesianTools)
library(pROC)
library(ggpubr)
library(mgcv)
library(purrr)

# Load useful functions
source("lib/abgFunctions.R")
 source("lib/traitspace_and_banquo.R")
 source("lib/interactionMatrix.R")
 source("lib/Likelihood.R")
 
# Load useful data
load("data/data_ready.Rdata")
load("data/traitspace_objs.Rdata")
cover_class <- read.table("data/cover_class.txt")
chains.to.analyze <- c(list.files("results", full.names = T, recursive = T, pattern = "posterior_objs.Rdata"))

# Complete species names
species <- c(AGRcap = "Agrostis capillaris", AMPflu = "Amphibromus fluitans", ANTodo = "Anthoxanthum odoratum", CARgau = "Carex gaudichaudiana", ELAacu = "Eleocharis acuta", ELApus =  "Eleocharis pusilla", EPIang = "Epilobium angustum", GALper = "Galium perpusillum", JUNart = "Juncus articulatus", LACstr = "Lachnagrostis striata", LAClya = "Lachnagrostis lyallii", LILrut = "Lilaeopsis ruthiana", LOBper = "Lobelia perpusilla", PARcan = "Parahebe canescens", PILoff = "Pilosella officinarum", PILpil = "Pilosella piloselloides")


### Extract or compute the useful features from the posterior distributions
out_results<- list()
for (ik in 1:length(chains.to.analyze)){
  chain = chains.to.analyze[ik]
  b=load(chain)

  pars <- apply(dis, 2, median)
  
  ### Compute the predicted site by species matrix, likelihood matrix at the median posterior
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
  
  out_results[[ik]] <- list(DIC =  dic.val$Dbar +dic.val$pV, conv = conv.val$mpsrf, alpha = int, alpha_beta = div_distri, ll = median(dis[,"Llikelihood"]), spxp = spxp.pred, roc = roc1, Nobs = prod(dim(obs.comm)), ll_025 = quantile(dis[,"Llikelihood"], 0.025), ll_975 = quantile(dis[,"Llikelihood"], 0.975))
  rm(out)
  rm(list = b)
  print(ik)
}
names(out_results) = chains.to.analyze

dat <- data.frame(names(out_results), conv = map_dbl(out_results, 'conv'), DIC = map_dbl(out_results, "DIC") )
lls <- map_dbl(out_results, 'll')
Nobs <- map_dbl(out_results, 'Nobs')
ll_H0 <- mean(lls[which(grepl("abio_FALSE/none", names(out_results)))])

# Extraction of the performance statistics
dat$pseudoR2 <- (1-exp(2/Nobs*(ll_H0-lls)))/(1-exp(2/Nobs*ll_H0))
dat$pseudoR2_025 <- (1-exp(2/Nobs*(ll_H0-map_dbl(out_results, 'll_025'))))/(1-exp(2/Nobs*ll_H0))
dat$pseudoR2_975 <- (1-exp(2/Nobs*(ll_H0-map_dbl(out_results, 'll_975'))))/(1-exp(2/Nobs*ll_H0))
dat$AUC  <- sapply( purrr::map(out_results, 'roc'), auc)

##### Generation of Table 1 ##########
dat2 <- dat
dat2 <- separate(as.data.frame(dat2), col = 1, into = c("folder", "abio", "bio", "chain"),sep= "/")
row.names(dat2) <- NULL

  
dat2$conv <- round(dat2$conv, 3)
dat2$DIC <- round(dat2$DIC, 1)
dat2[, c("pseudoR2_025","pseudoR2_975")] <- replace(dat2[, c("pseudoR2_025","pseudoR2_975")]  , dat2[, c("pseudoR2_025","pseudoR2_975")] < 0, 0)
dat2$pseudo_char <- paste0('[', format(round(dat2$pseudoR2_025,3), nsmall = 3), ', ', format(round(dat2$pseudoR2_975,3), nsmall = 3), "]")
dat2$pseudoR2 <- format(round(dat2$pseudoR2,3), nsmall = 3)
dat2$AUC <- format(round(dat2$AUC, 3), nsmall = 3)

dat2 <- select(dat2, -pseudoR2_025, -pseudoR2_975)
dat2$Assembly_model <- 'Null model'
dat2$Assembly_model[dat2$abio == "abio_FALSE" & dat2$bio != "none"] <- "Biotic model"
dat2$Assembly_model[dat2$abio == "abio_TRUE" & dat2$bio != "none"] <- "Abiotic & Biotic model"
dat2$Assembly_model[dat2$abio == "abio_TRUE" & dat2$bio == "none"] <- "Abiotic model"
dat2$Assembly_model <- factor(dat2$Assembly_model, levels= c('Null model',"Abiotic model","Biotic model", "Abiotic & Biotic model"), ordered = T)
dat2$bio <- factor(dat2$bio, levels = c("none","NoTr", "Height", "poros","SLA","2tr_HS","2tr_HP","2tr_PS","3tr"),  labels = c("No interactions","No traits", "Height","Porosity","SLA","Height + Porosity", "Height + SLA","Porosity + SLA" ,"Height + SLA + Porosity"), ordered = T)
x1 <- order(dat2$Assembly_model, dat2$bio)
dat2 <- dat2[x1,] 
dat <- dat[x1, ]
out_results <- out_results[x1]
dat2$Model <- dat$Model <- 1:nrow(dat2)
dat2 <- dat2[order(dat2$pseudoR2, decreasing = T),]

dat2 <- dat2[,c("Model","Assembly_model", "bio", "conv", "DIC", "pseudoR2","pseudo_char", "AUC")]
write.csv2(dat2, file = "results/tab.csv", row.names = F, quote = F)
################

###### Generation of Figure 2 ####
load("data/traitspace_objs.Rdata")
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

Fig2 <- ggplot(data_long, aes(x = DAYSUB, y = resp)) + 
  theme_bw() +
  geom_point(aes(size = Cover), shape = 1, color = "black") + 
  geom_line( aes(y = fit), size = 1, color = "#CC6677") +
  theme(axis.text = element_text(colour = "black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  labs(size = "Cover (%)", x = "Days submerged per year", y = "Trait value") +
  geom_ribbon( aes(ymin = lwr, ymax = upr, fill = "#CC6677"), alpha = 0.5) +
  guides(fill= 'none') +
  facet_wrap(~trait, scales = "free",labeller=as_labeller(lbs, label_parsed), ncol = 1, nrow = 3)
ggsave("results/Fig2.eps", Fig2, width = 82, height = 180, units = "mm", scale =1.5, bg = "white")

######

###### Generation of Figure 3 ###############
sel.models <- which(grepl("abio_TRUE", names(out_results)) & grepl("3tr|NoTr|Height|SLA|poros", names(out_results))) #Models to be plotted in Figure 3
max_int <- max(sapply(purrr::map(out_results[sel.models],'alpha'), function(x){
  diag(x) <- NA
  max(x, na.rm =T)
}))
tr <- t.avg[,c("maxht", "poros", "sla")]
diff_tr <- lapply(1:ncol(tr), function(i) outer(tr[,i], tr[,i], '-'))
names(diff_tr) <- c("maxht", "poros", "sla")

interaction_plots <- lapply(sel.models,function(i){
  
  if (!grepl("3tr", names(out_results)[i])){
    if (grepl("SLA",  names(out_results)[i])){
      diff = diff_tr$sla
      title_name <- paste0("SLA (Model ",i,")")
      xlabel <- bquote(paste('SLA log-ratio (log(',t[i],'/',t[j],'))'))
    }else if(grepl("Height",  names(out_results)[i])){
      diff = diff_tr$maxht
      title_name <- paste0("Height (Model ",i,")")
      xlabel <- bquote(paste('Height log-ratio (log(',t[i],'/',t[j],'))'))
    } else if (grepl("poros",  names(out_results)[i])){
      diff = diff_tr$poros
      title_name <- paste0("Root porosity (Model ",i,")")
      xlabel <- bquote(paste('Porosity difference (',t[i],'-',t[j],')'))
    }else{
      diff = diff_tr$maxht
      trait_name <- "Height"
      title_name <- paste0("No traits (Model ",i,")")
      xlabel <- bquote(paste('Height log-ratio (log(',t[i],'/',t[j],')'))
    }
    mat = out_results[[i]]$alpha
    diag(mat) <- NA
    dat.int = tibble(diff = as.numeric(diff), int = as.numeric(mat))
    gp <- ggplot(dat.int, aes(x = diff, y = int, fill = int))+ scale_fill_viridis_c(option = "viridis", limits = c(0, max_int))   +
      theme(panel.grid.major = element_blank(), plot.margin = unit(rep(1.5,4), "lines"), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+
      geom_point(shape = 21, size = 2.5) + ggtitle(title_name) + ylim(c(0,max_int)) +
      labs(y='Pairwise interaction', x = xlabel,fill = "Interaction\nstrength") + 
      geom_vline(xintercept = 0, linetype = 2, color = "black")
      return(gp)
  }else {
    title_name <- paste0("Height, SLA & Root porosity (Model ",i,")")
    
    mat = out_results[[i]]$alpha
    diag(mat) <- NA
    dat.int = tibble(diff1 = as.numeric(diff_tr$maxht), 
                     diff2 = as.numeric(diff_tr$poros), 
                     diff3 = as.numeric(diff_tr$sla), int = as.numeric(mat)) %>% drop_na()
    
    gp1 <- ggplot(dat.int, aes(x = diff3, y = diff1, fill = int))  +
      theme(panel.grid.major = element_blank(), plot.margin = unit(rep(1.5,4), "lines"), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
      geom_point(shape = 21, size = 2.5) +  geom_vline(xintercept = 0, linetype = 2, color = "black") + geom_hline(yintercept = 0, linetype = 2, color = "black")+
      labs(x=expression(paste('SLA log-ratio (log(',t[i],'/',t[j],'))')), y = expression(paste('Height log-ratio (log(',t[i],'/',t[j],'))')), fill = "Interaction\nstrength") + ggtitle(title_name ) + scale_fill_viridis_c(option = "viridis", limits = c(0, max_int), direction = 1) +
      annotate("point", x = as.numeric(filter(dat.int, int == max(int))[1,"diff3"]), y = as.numeric(filter(dat.int, int == max(int))[1,"diff1"]), fill = "red", shape=23, size = 4.5)
    gp2 <- ggplot(dat.int, aes(x = diff3, y = diff2, fill = int))  +
      theme(panel.grid.major = element_blank(), plot.margin = unit(rep(1.5,4), "lines"), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"))+ 
      geom_point(shape = 21, size = 2.5) +  geom_vline(xintercept = 0, linetype = 2, color = "black") + geom_hline(yintercept = 0, linetype = 2, color = "black")+
      labs(x=expression(paste('SLA log-ratio (log(',t[i],'/',t[j],'))')), y = bquote(paste('Porosity difference (',t[i],'-',t[j],')')), fill = "Interaction\nstrength") + ggtitle(title_name ) + scale_fill_viridis_c(option = "viridis", limits = c(0, max_int), direction = 1) +
      annotate("point", x = as.numeric(filter(dat.int, int == max(int))[1,"diff3"]), y = as.numeric(filter(dat.int, int == max(int))[1,"diff2"]), fill = "red", shape=23, size = 4.5)  }
  return(list(gp1, gp2))
})
names(interaction_plots) <- map_chr(strsplit(names(out_results)[sel.models],"/"),3)

interaction_plots$'3tr_1' <- interaction_plots$`3tr`[[1]]
interaction_plots$'3tr_2' <- interaction_plots$`3tr`[[2]]
Fig3 <- ggarrange(plotlist = interaction_plots[c("NoTr", "Height", "poros","SLA", "3tr_1", "3tr_2")], nrow = 3, ncol = 2, common.legend = T, legend = 'right', labels = LETTERS[1:length(interaction_plots)]) 

ggsave("results/Fig3.eps", Fig3, width = 110, height = 130, units = "mm", scale = 2, bg = "white")
#######

####### Generation of Figure 4 ################

id.abio <- which(grepl("abio_TRUE", names(out_results)) & grepl("none", names(out_results)))
id.best <- which(grepl("abio_TRUE", names(out_results)) & grepl("3tr", names(out_results)))
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

dat.spDis <- out_results %>% format_data(id.abio = id.abio, id.best = id.best) 
Fig4 <- ggplot(dat.spDis, aes(x = DAYSUB, y = Cover*100)) +
  geom_point(data = filter(dat.spDis,Model== 'Observed'), fill = 'black') + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
  geom_smooth(mapping = aes(color = Model), se = F) + 
  scale_color_manual(values=c("#000000", "#CC6677", "#117733")) +
  facet_wrap(~Species, scales = 'free') + scale_y_continuous(limits =c(0,80),breaks = c(0,5,20,50,80)) +  coord_trans(y = "sqrt") +
  labs(x = "Days submerged per year", y = "Species cover (%)") + guides(linetype = "none") 
ggsave("results/Fig4.eps", Fig4, width = 173, height = 173, units = "mm", scale = 1.2, bg = "white")
################

#### Generation of Supplementary Figure 4 ############

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
  out <- x %>% purrr::map(pos) %>% bind_rows(.id = "Chain") %>% 
    separate(Chain, sep = "/", remove = F, into = c("folder","abio", "bio_traits", "chain_number")) %>%
    mutate(Assembly_model = fct_cross(abio, bio_traits, sep='_'))
  out <- out %>% mutate(Assembly_model = fct_collapse(Assembly_model, 
                                           "Null model" = "abio_FALSE_none",
                                           "Abiotic model" ="abio_TRUE_none",
                                           'Biotic models' = str_subset(levels(Assembly_model),"abio_FALSE_(?!none)"),
                                           'Abiotic & Biotic models' = str_subset(levels(Assembly_model),"abio_TRUE_(?!none)")),
                        bio_traits = factor(bio_traits, levels = c('none', 'NoTr','Height', 'SLA',"poros", '2tr_rho', "2tr_HP_rho", "2tr_SP_rho", "3tr_rho")))
  out <- out %>% mutate(Assembly_model = factor(Assembly_model, levels = c("Null model", "Abiotic model",'Biotic models', 'Abiotic & Biotic models')),
                        bio_traits = fct_recode(bio_traits, "No interactions" = "none", 
                                   "No traits"="NoTr",
                                   'Porosity' = "poros",
                                    "Height & SLA" = "2tr_rho",
                                     "Height & Porosity" = "2tr_HP_rho",
                                   "Porosity & SLA" = "2tr_SP_rho",
                                   "Height, Porosity & SLA" = "3tr_rho")) %>%
          select(-folder, -abio, -chain_number, -Chain)
    return(out)
}

dat.div <- out_results %>% format_data2('alpha_beta') %>% pivot_longer(cols = c("Beta", 'mAlpha', 'Beta_simul', 'mAlpha_simul'), names_to = "diversity_name", values_to = "diversity")
dat.div$diversity_name <- factor(dat.div$diversity_name, levels = c("Beta", "Beta_simul", "mAlpha_simul", "mAlpha"),
                              labels = c(expression(paste(beta,'-diversity (fixed)')),
                                         expression(paste(beta,'-diversity (stochastic)')),
                                         expression(paste(alpha,'-diversity (stochastic)')),
                                         expression(paste(alpha,'-diversity (fixed)'))))
ref_tab = as_tibble(data.frame(diversity_name = levels(dat.div$diversity_name), ref_line = c(ab.obs.med[["Beta"]], ab.obs[['Beta']],ab.obs[['mAlpha']],ab.obs.med[["mAlpha"]])))

FigS4 <- ggplot(dat.div, aes_string(y= 'diversity',x = 'Assembly_model', fill = 'bio_traits')) + 
  theme_bw() +
  geom_boxplot(outlier.shape = NA) + 
  geom_hline(data= ref_tab,aes(yintercept = ref_line), linetype="dashed", size = 1) +
    facet_wrap( ~diversity_name, drop = T, scale = 'free_y', labeller = label_parsed) +
  theme(axis.title.y=element_blank(), axis.title.x=element_blank(), axis.text = element_text(colour = "black"),axis.text.x = element_text(angle = 35, hjust = 1)) +
  labs(fill = "Competition traits") +scale_fill_manual(values=c("#FFFFFF",carto_pal(n=9, name="Safe")))
ggsave("results/FigS4.jpg", FigS4, width = 8, height = 6.5, bg = "white")

#############

##  Generation of Supplementary Figure 5 ######
p1 <- ggplot(dat.spDis, aes(x = Cover*100, col = Model)) + scale_color_manual(values = c("#000000", "#CC6677", "#117733")) +
  theme(legend.position = "top", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + labs(x = 'Species cover (%)') + 
  scale_x_continuous(trans = 'sqrt', breaks = c(0.5,2, 5,10,20,40,80), labels = as.character(c(0.5,2, 5,10,20,40,80)))+
  geom_density(size = 1, adjust = 2.5 )

dat3 <- dat.spDis%>% select(-DAYSUB) %>% pivot_wider(values_from = "Cover", names_from = "Model") %>% 
  mutate(Observed = factor(Observed*100), `Abiotic model` = 100*`Abiotic model`, `Abiotic + Biotic model` = 100*`Abiotic + Biotic model`)

p2 <- ggplot(dat3, aes(y = `Abiotic + Biotic model`, x = Observed)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + 
  labs(x = 'Observed cover class (%)', y = "Percent cover (Abiotic & biotic model)") + ylim(c(0, max(dat3$`Abiotic model`, dat3$`Abiotic + Biotic model`))) +
  geom_boxplot(varwidth = FALSE, fill =  "#117733", alpha = 0.5)

p3 <- ggplot(dat3, aes(y = `Abiotic model`, x = Observed)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) + 
  labs(x = 'Observed cover class (%)', y = "Percent cover (Abiotic model)") + ylim(c(0, max(dat3$`Abiotic model`, dat3$`Abiotic + Biotic model`))) +
  geom_boxplot(varwidth = FALSE, fill = "#CC6677", alpha = 0.5) 
p4 <- ggplot(dat3, aes(x = `Abiotic model`, y = `Abiotic + Biotic model`)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black"))  + 
  geom_point(shape = 19, alpha = 0.3) + labs(x = 'Percent cover (Abiotic model)', y = "Percent cover (Abiotic & biotic model)")
FigS5 <- ggarrange(plotlist = list(p3, p2, p4, p1), heights = c(1,1), widths = c(1,1), nrow = 2, ncol = 2, labels = c('A','B',  'C',  'D'), common.legend = T) 
ggsave("results/FigS5.jpg", FigS5, width = 9, height = 8, bg = "white")
###########

## Generation of Supplementary Figure 6 : ROC curves ######
rocs <- purrr::map(out_results, 'roc')
rocs.dat <- bind_rows(lapply(rocs, function(x) data.frame(spec= 1-x$specificities, sens = x$sensitivities)), .id = "Chain") %>% 
  separate(Chain, sep = "/", remove = F, into = c("folder","abio", "bio_traits", "chain_number")) %>%
  mutate(Assembly_model = fct_cross(abio, bio_traits, sep='_')) %>% 
  mutate(Assembly_model = fct_collapse(Assembly_model, 
                                       "Null model" = "abio_FALSE_none",
                                       "Abiotic model" ="abio_TRUE_none",
                                       'Biotic models' = str_subset(levels(Assembly_model),"abio_FALSE_(?!none)"),
                                       'Abiotic & Biotic models' = str_subset(levels(Assembly_model),"abio_TRUE_(?!none)"))) %>%
  mutate(Assembly_model = factor(Assembly_model, levels = c("Null model", "Abiotic model",'Biotic models', 'Abiotic & Biotic models')))

cols = c("#888888","#CC6677","#332288","#117733") 
FigS6 <- ggplot(rocs.dat, aes(x = spec, y = sens, group = Chain, color = Assembly_model)) + 
  theme(legend.position = "bottom", legend.title=element_blank(), panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.text = element_text(colour = "black"), axis.line = element_line(colour = "black")) +
  geom_line() + labs(x = "False positive rate", y = 'True positive rate') + 
  scale_colour_manual(name = "grp",values = cols) 
ggsave("results/FigS6.jpg", FigS6, width = 6, height = 6, bg = "white")


## Generation of Supplementary Figure for the posterior distributions of each assembly models######
library(extraDistr)
names_params <- c("a","b", expression(phi), 'C', 
                  expression(mu[1]), expression(sigma[1]), 
                  expression(mu[2]), expression(sigma[2]), 
                  expression(mu[3]), expression(sigma[3]), 
                  expression(rho[1]), expression(rho[2]), expression(rho[3]))
all.params <- c("a2", "c2", "phi", "intercept","mu1","sigma1","mu2","sigma2", "mu3", "sigma3", "rho1", "rho2", "rho3")
names(all.params) <- names_params
names_models <- dat$Model

for (ip in 1:length(chains.to.analyze)){
  chain = dat$names.out_results.[ip]
  b=load(chain)
  source("main/priors.R")
  ### Isolate posterior
  post.dis <- as.data.frame(dis[,colnames(dis) %in% all.params])
  colnames(post.dis) <- list_params
  post.dis$distri <- 'Posterior'
  prior.dis <- data.frame(sampler(nrow(post.dis)))
  colnames(prior.dis) <- list_params
  for (i in 1:ncol(prior.dis)){
    prior.dis[prior.dis[,i] < bounds[i,1] |     prior.dis[,i] > bounds[i,2],i] <- NA
  }
  prior.dis[prior.dis$phi > 10,"phi"] <- NA
  prior.dis$distri <- 'Prior'
  bounds$param = factor(row.names(bounds), levels=all.params[all.params %in% list_params])
  bounds<-bounds %>% mutate(param = fct_recode(param, !!!all.params[all.params %in% all_of(list_params)]))
  distris <- rbind(prior.dis, post.dis)
   distris_long <- distris %>% pivot_longer(cols = list_params, names_to = "param") %>% 
    mutate(param = fct_recode(param, !!!all.params[all.params %in% all_of(list_params)])) %>% 
    mutate(param = factor(param, levels = names(all.params)[all.params %in% list_params]))
  levels(distris_long$param) <- replace(levels(distris_long$param), levels(distris_long$param) == "rho[1]", expression(rho[21]))
  levels(distris_long$param) <- replace(levels(distris_long$param), levels(distris_long$param) == "rho[2]", expression(rho[31]))
  levels(distris_long$param) <- replace(levels(distris_long$param), levels(distris_long$param) == "rho[3]", expression(rho[32]))
  levels(bounds$param) <- replace(levels(bounds$param), levels(bounds$param) == "rho[1]", expression(rho[21]))
  levels(bounds$param) <- replace(levels(bounds$param), levels(bounds$param) == "rho[2]", expression(rho[31]))
  levels(bounds$param) <- replace(levels(bounds$param), levels(bounds$param) == "rho[3]", expression(rho[32]))
  
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
  
  ggsave(gp1, file = paste0("results/post_models/post_Model_",as.character(names_models[ip]) ,".jpg"),width = 8.16,height = height,units = 'in')
  print(ip)
}

