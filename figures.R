
rm(list = ls())

# Load libraries
library(ade4)
library(plyr)
library(shape)

# Load data
load("data.Rdata")
load("results.Rdata")

t.avg <- ddply(KH.data, c("species"), summarise,
               maxht  = quantile(height,0.975), sla = mean(sla, na.rm = T))
t.avg$maxht <- log(t.avg$maxht)

dudi.tr <- dudi.pca(t.avg[,c("maxht", "sla")], nf = 2, scannf = F)
t.avg$Axis1 <- dudi.tr$li[,1]
t.avg$Axis2 <- dudi.tr$li[,2]

### Figure 2 : Visualisation of community structure along the environmental gradient ############

traitspace_graph <- function(ComMat, env, xlab = "", f = 12,cols = NULL, ylim= NULL, sqrt = F, cex.legend = 1,cex.lab = 1.1,cex.axis = 1,cex.point = 1,...){
  require(mgcv)
  if (is.null(cols)) cols <- rainbow(ncol(ComMat))
  # layout(matrix(c(1,1,2), nrow = 1))
  if (is.null(ylim))  ylim=c(0.001,max(ComMat))
  if (sqrt) ComMat <- sqrt(ComMat)
  plot(env[,1], ComMat[,1], type='n', xlim=range(env[,1]),ylim = ylim,
       xlab='', 
       col="white", font.lab=1, font=1, pch = 19, cex.lab = cex.lab, cex.axis = cex.axis, cex = cex.point,...)
  mtext(side = 1,xlab, line = 3, cex = cex.axis)
  for( i in which(colSums(ComMat) > 0)){
    gradient = env[,1]
    mod1.gam <- gam(ComMat[,i] ~ s(gradient,k=f,bs="tp"))
    
    xvec =  seq(min(gradient), max(gradient), length.out = 1000)
    yvec <- predict(mod1.gam ,list(gradient =xvec))
    yvec[yvec <0] <-0
    lines(xvec, yvec, col=cols[i], lwd=4, lty=1) 
  }
  for( i in which(colSums(ComMat) > 0)){
    points(env[,1],ComMat[,i], col=cols[i], pch=16, cex = cex.point)
  }
}

jpeg("Figure2.jpeg", height = 7.5, width = 10, res = 600, units = "in", quality = 1)
par(oma = c(1,1,1,1), cex.lab = 2, cex.axis = 2, mar = c(4,4.5,5,2), cex = 1.8, mfrow = c(2,2))
layout(matrix(c(1,1,1,2,2,2,6,6,6, 3,3,3,4,4,4,5,5,5), nrow = 2, ncol = 9, byrow = T))
f.param = 5
cex.lab = 1.3
cex.main = 1.4
cex.axis = 1
cex.point = 1
cols <- rainbow(ncol(P_S_E_tr))

traitspace_graph(comm,env = env, xlab = 'Days submerged', f= f.param, ylim = c(0,1), sqrt =T, cols = cols, main = "Observed", cex.lab = cex.lab, ylab = expression('Species relative abundance'^'1/2'), cex.main = cex.main, cex.axis = cex.axis, cex.point = cex.point)
mtext("(A)", 2, adj=6, las=1, padj=-10, line = -8, cex =1.3, font = 2)

traitspace_graph(P_S_E_tr,env = env, xlab = 'Days submerged', f= f.param, ylim = c(0,1), sqrt = T, cols = cols, main = "Abiotic filtering (TM)", cex.lab = cex.lab, ylab = expression('Species relative abundance'^'1/2'), cex.main = cex.main, cex.axis = cex.axis, cex.point = cex.point)
mtext("(B)", 2, adj=6, las=1, padj=-10, line = -8, cex =1.3, font = 2)

traitspace_graph(banquo_out_S$P_S_E_interactions,env = env, xlab = 'Days submerged', f= f.param, ylim = c(0,1), sqrt = T, cols = cols, main = "Banquo (SLA)", cex.lab = cex.lab, ylab = expression('Species relative abundance'^'1/2'), cex.main = cex.main)
mtext("(C)", 2, adj=6, las=1, padj=-10, line = -10, cex = 1.3, font = 2)

traitspace_graph(banquo_out_H$P_S_E_interactions,env = env, xlab = 'Days submerged', f= f.param, ylim = c(0,1), sqrt = T, cols = cols, main = "Banquo (Height)", cex.lab = cex.lab, ylab = expression('Species relative abundance'^'1/2'), cex.main = cex.main)
mtext("(D)", 2, adj=6, las=1, padj=-10, line = -10, cex = 1.3, font = 2)

traitspace_graph(banquo_out_HS2$P_S_E_interactions,env = env, xlab = 'Days submerged', f= f.param, ylim = c(0,1), sqrt = T, cols = cols, main = "Abiotic & biotic filtering\n(BM-HS)", 
                   cex.lab = cex.lab, ylab = expression('Species relative abundance'^'1/2'), cex.main = cex.main, cex.axis = cex.axis, cex.point = cex.point)
mtext("(E)", 2, adj=6, las=1, padj=-10, line = -8, cex =1.3, font = 2)
par(xpd = T)
plot(c(0,2,3,5), type = "n", bty ="none", xaxt = "n",yaxt = "n", ylab ="", xlab = "", mar = c(0,0,0,0))
legend(1.1,5.5, legend = species[,1], col = cols, pch = 19, cex = 1.2, bty = "n" )
par(xpd = F)
dev.off()

### Figure 3 + extension : Visualisation of the interaction matrices ############
height_diff <- outer(scale(log(t.avg$maxht)), scale(log(t.avg$maxht)), "-")
sla_diff <- outer(scale(t.avg$sla), scale(t.avg$sla), "-")
Aij_h <- banquo_out_H$alphas
Aij_s <- banquo_out_S$alphas
Aij_hs <- banquo_out_HS2$alphas

jpeg("Figure3.jpeg", width = 17, height = 6.61, units = "cm", res= 400)
  cols <- colorRampPalette(c("white","skyblue","blue","purple", "red"))
  par(cex = 0.4, mar = c(3,3.5,2,0.7), cex.main = 1.15, cex.axis = 1, oma = c(1,1,1,1))
  layout(matrix(c(1,1,1,2,2,2,3,3,3,4), nrow = 1, byrow = T))
  aij_vec.tot <- c(as.numeric(Aij_s), as.numeric(Aij_h), as.numeric(Aij_hs))
  cex.lab <- 0.7
  cex.point = 1.3
  
  aii_h = Aij_h[1,1]
  Aij_h2 <- Aij_h
  diag(Aij_h2) <- NA
  plot(height_diff, Aij_h2,  bg = cols(14)[cut(Aij_h2, 14)],
       pch = 21, cex =cex.point, lwd = 1, col = "black",
       xlab = "", ylab = '', main = "BM-H (Height)", cex.lab = cex.lab )
  points(0,aii_h, bg = "yellow", cex = cex.point*1.5, pch = 23)
  abline( v= 0, lty = 2)
  mtext(side = 1, text = "log-ratio of height", cex = cex.lab, line = 2.5)
  mtext(side = 2, text = expression('a'['ij']), cex = cex.lab*1.2, line = 2)
  mtext("(A)", 2, adj=6, las=1, padj=-9.5, line = -5, cex = 0.8, font = 2)
  
  aii_s = Aij_s[1,1]
  Aij_s2 <- Aij_s
  diag(Aij_s2) <- NA
  plot(sla_diff, Aij_s2,  bg = cols(14)[cut(as.numeric(Aij_s2), 14)],
       pch = 21, cex =cex.point, lwd = 1,
       xlab = "", ylab = '', main = "BM-S (SLA)", cex.lab = cex.lab )
  points(0, aii_s, bg = "yellow", cex = cex.point*1.5, pch = 23)
  mtext(side = 1, text = "SLA difference", cex = cex.lab, line = 2.5)
  mtext(side = 2, text = expression('a'['ij']), cex = cex.lab*1.2, line = 2)
  abline( v= 0, lty = 2)
  mtext("(B)", 2, adj=6, las=1, padj=-9.5, line = -5, cex = 0.8, font = 2)
  
  aii_hs = round(Aij_hs[1,1],3)
  Aij_hs2 <- Aij_hs
  diag(Aij_hs2) <- NA
  plot(as.numeric(height_diff), as.numeric(sla_diff), bg = cols(14)[cut(Aij_hs2, 14)], pch = 21, cex = cex.point,xlab = "", ylab = "", cex.lab = cex.lab, main = "BM-HS (SLA + Height)")
  abline( v= 0, h = 0, lty = 2)
  text(max(height_diff), 0.95* max(sla_diff), bquote(paste('a'['ii'], '=', .(aii_hs))), pos = 2, cex = 1.2)
  mtext("(C)", 2, adj=6, las=1, padj=-9.5, line = -5, cex = 0.8, font = 2)
  mtext(side = 1, text = "log-ratio of height", cex = cex.lab, line = 2.5)
  mtext(side = 2, text = "SLA difference", cex = cex.lab, line = 2.5)
  plot(0,1, type ="n", axes = F, xlab ="", ylab = "")
  colorlegend(col =cols(100), zlim = c(0,1), zval = c(0,0.5,1), digit = 2, posx = c(0.30,0.45),posy = c(0.05,0.8),main = "", cex = 1., main.cex =1.4, font = 1)
  text(0.45, 1.5, labels = expression('a'['ij']), cex = 1.2, font = 2)
  
dev.off()
##############