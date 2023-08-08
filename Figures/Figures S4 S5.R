################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figures S4 and S5.


######################### 1. SETTINGS ##########################################


## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.
if(!require(scales)){
  install.packages("scales")
  library(scales)
}

if(!require(purrr)){
  install.packages("purrr")
}
library(purrr)

if(!require(grid)){
  install.packages("grid")
}
library(grid)

if(!require(dplyr)){
  install.packages("dplyr")
}
library(dplyr)

## Source the functions.
source("Code/genfun.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")



######################### 2. CALCULATIONS ######################################
msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}

lisUn<-mat1
if (T){
  
  #lisUexv2<-lapply(lisUn, function(x){ genfun(t(x$large), x$res_abun)})
  lisUexv2<-lapply(lisUn, function(x){ genfun(t(x$preference*1000000), rep(1,nrow(x$preference)))})
  lisUexv2<-bind_rows(lisUexv2)
  
  lisUexv2_10<-lapply(lisUn, function(x){ genfun(t(x$small1), x$res_abun)})
  lisUexv2_10<-bind_rows(lisUexv2_10)
  
  lisUexv2_50<-lapply(lisUn, function(x){ genfun(t(x$small2), x$res_abun)})
  lisUexv2_50<-bind_rows(lisUexv2_50)
  
  lisUexv2_100<-lapply(lisUn, function(x){ genfun(t(x$small3), x$res_abun)})
  lisUexv2_100<-bind_rows(lisUexv2_100)
  
  lisUexv2_500<-lapply(lisUn, function(x){ genfun(t(x$small4), x$res_abun)})
  lisUexv2_500<-bind_rows(lisUexv2_500)
  
}

l5<-nrow(lisUexv2)/3
l10<-nrow(lisUexv2)/3*2
l50<-nrow(lisUexv2)

summ_Bs<-data.frame(S10=lisUexv2_10$Bs-lisUexv2$Bs, 
                    S50=lisUexv2_50$Bs-lisUexv2$Bs, 
                    S100=lisUexv2_100$Bs-lisUexv2$Bs,
                    S500=lisUexv2_500$Bs-lisUexv2$Bs, 
                    res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_Bprime<-data.frame(S10=lisUexv2_10$`B'`-lisUexv2$`B'`, 
                        S50=lisUexv2_50$`B'`-lisUexv2$`B'`, 
                        S100=lisUexv2_100$`B'`-lisUexv2$`B'`,
                        S500=lisUexv2_500$`B'`-lisUexv2$`B'`, 
                        res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_W<-data.frame(S10=lisUexv2_10$W-lisUexv2$W, 
                   S50=lisUexv2_50$W-lisUexv2$W, 
                   S100=lisUexv2_100$W-lisUexv2$W,
                   S500=lisUexv2_500$W-lisUexv2$W, 
                   res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_PS<-data.frame(S10=lisUexv2_10$PS-lisUexv2$PS, 
                    S50=lisUexv2_50$PS-lisUexv2$PS, 
                    S100=lisUexv2_100$PS-lisUexv2$PS,
                    S500=lisUexv2_500$PS-lisUexv2$PS, 
                    res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_FT<-data.frame(S10=lisUexv2_10$FT-lisUexv2$FT, 
                   S50=lisUexv2_50$FT-lisUexv2$FT, 
                   S100=lisUexv2_100$FT-lisUexv2$FT,
                   S500=lisUexv2_500$FT-lisUexv2$FT, 
                   res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_d<-data.frame(S10=lisUexv2_10$`1-d'`-lisUexv2$`1-d'`, 
                   S50=lisUexv2_50$`1-d'`-lisUexv2$`1-d'`, 
                   S100=lisUexv2_100$`1-d'`-lisUexv2$`1-d'`,
                   S500=lisUexv2_500$`1-d'`-lisUexv2$`1-d'`, 
                   res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

summ_gen<-data.frame(S10=lisUexv2_10$gen-lisUexv2$gen, 
                    S50=lisUexv2_50$gen-lisUexv2$gen, 
                    S100=lisUexv2_100$gen-lisUexv2$gen,
                    S500=lisUexv2_500$gen-lisUexv2$gen, 
                    res=c(rep(5, l5),rep(15, l5), rep(55, l5)))


######################### 3. PLOTTING ######################################

x1_xit<-jitter(rep(1, l5), factor=6)
x2_xit<-jitter(rep(2, l5), factor=6)
x3_xit<-jitter(rep(3, l5), factor=6)

x_xit<-c(x1_xit,x2_xit,x3_xit)
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
spar_col<-as.raster(matrix(colfunc(2000), nrow=1))
spar_col<-rep(rep(spar_col, each=1), 3)

############## S4
png(filename="Figures/Exported/FigureS4.png", width=6000, height=5000, res=600)
if (T){
  
  par(mar= c(1,3,2,0),las=1, xpd=F)
  layout(matrix(c(25,20,20,20,20,19,
                  17, 1, 2, 3, 4,19,
                  17,21,21,21,21,19, 
                  17, 5, 6, 7, 8,19,
                  17,22,22,22,22,19,
                  17, 9,10,11,12,19,
                  17,23,23,23,23,19,
                  17,13,14,15,16,19,
                  24,18,18,18,18,19), ncol=6, byrow=T),
         widths=c(6,45/2,45/2,45/2,45/2,11), heights=c(6,40,6,40,6,40,6,40,10))
  layout.show(25)
  ## Bs
  
  boxplot(S10~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), main="10 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_Bs$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bs[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bs)[which(summ_Bs$res==5),1])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==15),1])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="50 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bs$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bs[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bs)[which(summ_Bs$res==5),2])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==15),2])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="100 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bs$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bs[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bs)[which(summ_Bs$res==5),3])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==15),3])/l5,2),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="500 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bs$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_Bs, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bs[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bs)[which(summ_Bs$res==5),1])/l5,4),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==15),1])/l5,4),
         round(sum(is.na(summ_Bs)[which(summ_Bs$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ## Bprime
  
  
  boxplot(S10~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_Bprime$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bprime[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==5),1])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==15),1])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bprime$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bprime[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==5),2])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==15),2])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bprime$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bprime[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==5),3])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==15),3])/l5,2),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_Bprime$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_Bprime, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_Bprime[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==5),1])/l5,4),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==15),1])/l5,4),
         round(sum(is.na(summ_Bprime)[which(summ_Bprime$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ## W
  
  boxplot(S10~res, data=summ_W, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_W$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_W, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_W[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_W)[which(summ_W$res==5),1])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==15),1])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_W$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_W[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_W)[which(summ_W$res==5),2])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==15),2])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_W$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_W[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_W)[which(summ_W$res==5),3])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==15),3])/l5,2),
         round(sum(is.na(summ_W)[which(summ_W$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_W$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_W, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_W[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_W)[which(summ_W$res==5),1])/l5,4),
         round(sum(is.na(summ_W)[which(summ_W$res==15),1])/l5,4),
         round(sum(is.na(summ_W)[which(summ_W$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ### PS
  
  
  boxplot(S10~res, data=summ_PS, ylim=c(-1,1),  pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_PS$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_PS, ylim=c(-1,1),  pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_PS[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_PS)[which(summ_PS$res==5),1])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==15),1])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_PS, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_PS$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_PS, ylim=c(-1,1), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_PS[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_PS)[which(summ_PS$res==5),2])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==15),2])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_PS, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_PS$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_PS, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_PS[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_PS)[which(summ_PS$res==5),3])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==15),3])/l5,2),
         round(sum(is.na(summ_PS)[which(summ_PS$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_PS, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_PS$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_PS, ylim=c(-1,1), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_PS[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_PS)[which(summ_PS$res==5),1])/l5,4),
         round(sum(is.na(summ_PS)[which(summ_PS$res==15),1])/l5,4),
         round(sum(is.na(summ_PS)[which(summ_PS$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ### Y Label
  
  par(mar= c(0,0.5,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.5, labels="Difference between the observed and expected degree of generalization", cex=1, srt=90)
  
  ### X Label
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels="Number of potential resources", cex=1)
  
  
  ###### Legend
  
  par(mar= c(4,2,4,0),las=1)
  plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialization\nparameter",
       cex.main=1, ylim=c(0,1))
  legend_image <- as.raster(rev(colfunc(20)))
  grid.raster(legend_image, width=0.03, height = 0.8, x = unit(0.94, "npc"))
  par(xpd=T)
  abline(v=-0.2, lty=3)
  text(y=c(-0.02,1.02), x =0.7 , labels = c(0.1, "60"), cex=1)
  par(xpd=F)
 
  
  ### Bs main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels=expression(italic("Bs'")), cex=2)
  
  ### Bprime main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.3, labels=expression(italic("B''")), cex=2)
  
  ### W main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.3, labels=expression(italic("W'")), cex=2)
  
  ### PS main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.3, labels=expression(italic("PS'")), cex=2)
}
dev.off()


############## S5 ##########################################################################
png(filename="Figures/Exported/FigureS5.png", width=6000, height=4200, res=600)
if (T){
  
  par(mar= c(1,3,2,0),las=1, xpd=F)
  layout(matrix(c(20,16,16,16,16,15,
                  13, 1, 2, 3, 4,15,
                  13,17,17,17,17,15, 
                  13, 5, 6, 7, 8,15,
                  13,18,18,18,18,15,
                  13, 9,10,11,12,15,
                  19,14,14,14,14,15), ncol=6, byrow=T),
         widths=c(6,45/2,45/2,45/2,45/2,11), heights=c(6,40,6,40,6,40,10))
  
  ## FT
  
  boxplot(S10~res, data=summ_FT, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), main="10 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_FT$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_FT, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_FT[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_FT)[which(summ_FT$res==5),1])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==15),1])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="50 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_FT$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_FT[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_FT)[which(summ_FT$res==5),2])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==15),2])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="100 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_FT$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_FT[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_FT)[which(summ_FT$res==5),3])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==15),3])/l5,2),
         round(sum(is.na(summ_FT)[which(summ_FT$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="500 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_FT$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_FT, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_FT[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_FT)[which(summ_FT$res==5),1])/l5,4),
         round(sum(is.na(summ_FT)[which(summ_FT$res==15),1])/l5,4),
         round(sum(is.na(summ_FT)[which(summ_FT$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ## d
  
  
  boxplot(S10~res, data=summ_d, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_d$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_d, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_d[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_d)[which(summ_d$res==5),1])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==15),1])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_d$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_d[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_d)[which(summ_d$res==5),2])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==15),2])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_d$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_d[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_d)[which(summ_d$res==5),3])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==15),3])/l5,2),
         round(sum(is.na(summ_d)[which(summ_d$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_d$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_d, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_d[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_d)[which(summ_d$res==5),1])/l5,4),
         round(sum(is.na(summ_d)[which(summ_d$res==15),1])/l5,4),
         round(sum(is.na(summ_d)[which(summ_d$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ## gen
  
  boxplot(S10~res, data=summ_gen, ylim=c(-1,1), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_gen$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_gen, ylim=c(-1,1), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_gen[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_gen)[which(summ_gen$res==5),1])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==15),1])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_gen, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_gen$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_gen, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_gen[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_gen)[which(summ_gen$res==5),2])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==15),2])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_gen, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_gen$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_gen, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_gen[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_gen)[which(summ_gen$res==5),3])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==15),3])/l5,2),
         round(sum(is.na(summ_gen)[which(summ_gen$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_gen, ylim=c(-1,1), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_gen$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_gen, ylim=c(-1,1),  yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_gen[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_gen)[which(summ_gen$res==5),1])/l5,4),
         round(sum(is.na(summ_gen)[which(summ_gen$res==15),1])/l5,4),
         round(sum(is.na(summ_gen)[which(summ_gen$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ### Y Label
  
  par(mar= c(0,0.5,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.5, labels="Difference between the observed and expected degree of generalization", cex=1, srt=90)
  
  ### X Label
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels="Number of potential resources", cex=1)
  
  
  ###### Legend
  
  par(mar= c(4,2,4,0),las=1)
  plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialization\nparameter",
       cex.main=1, ylim=c(0,1))
  legend_image <- as.raster(rev(colfunc(20)))
  grid.raster(legend_image, width=0.03, height = 0.8, x = unit(0.94, "npc"))
  par(xpd=T)
  abline(v=-0.2, lty=3)
  text(y=c(-0.02,1.02), x =0.7 , labels = c(0.1, "60"), cex=1)
  par(xpd=F)
  
  
  ### FT main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels=expression(italic("FT'")), cex=2)
  
  ### d main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels=expression(italic("1-d'")), cex=2)
  
  ### gen main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels=expression(italic("gen'")), cex=2)
  

}
dev.off()
