################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure 3.


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

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/wcfun.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")
mat2<-readRDS("Data/sim_data2.rds")


######################### 2. CALCULATIONS ######################################
msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}

biasfun<-function(exp,obs){
  bias<-(mean(obs)-exp)
  return(bias)
} #Bias



lisUn<-mat1
if (T){

  lisUexv2<-lapply(lisUn, function(x){ alpha_PDI(t(x$preference), rep(1,nrow(x$preference)), corrected=F)})
  lisUexv2<-unlist(lisUexv2)
  
  lisUexv2_10<-lapply(lisUn, function(x){ alpha_PDI(t(x$small1), x$res_abun)$corrected_aPDI})
  lisUexv2_10<-unlist(lisUexv2_10)
  
  lisUexv2_50<-lapply(lisUn, function(x){ alpha_PDI(t(x$small2), x$res_abun)$corrected_aPDI})
  lisUexv2_50<-unlist(lisUexv2_50)
  
  lisUexv2_100<-lapply(lisUn, function(x){ alpha_PDI(t(x$small3), x$res_abun)$corrected_aPDI})
  lisUexv2_100<-unlist(lisUexv2_100)
  
  lisUexv2_500<-lapply(lisUn, function(x){ alpha_PDI(t(x$small4), x$res_abun)$corrected_aPDI})
  lisUexv2_500<-unlist(lisUexv2_500)

}


l5<-length(lisUexv2)/3
l10<-length(lisUexv2)/3*2
l50<-length(lisUexv2)




summ_aPDI<-data.frame(S10=lisUexv2_10-lisUexv2, 
                      S50=lisUexv2_50-lisUexv2, 
                      S100=lisUexv2_100-lisUexv2,
                      S500=lisUexv2_500-lisUexv2, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))


if (T){
  
  lisUexv2_raw<-lapply(lisUn, function(x){ alpha_PDI(t(x$preference), rep(1,nrow(x$preference)), corrected=F)})
  lisUexv2_raw<-unlist(lisUexv2_raw)
  
  lisUexv2_10_raw<-lapply(lisUn, function(x){ alpha_PDI(t(x$small1), x$res_abun, corrected=F)})
  lisUexv2_10_raw<-unlist(lisUexv2_10_raw)
  
  lisUexv2_50_raw<-lapply(lisUn, function(x){ alpha_PDI(t(x$small2), x$res_abun, corrected=F)})
  lisUexv2_50_raw<-unlist(lisUexv2_50_raw)
  
  lisUexv2_100_raw<-lapply(lisUn, function(x){ alpha_PDI(t(x$small3), x$res_abun, corrected=F)})
  lisUexv2_100_raw<-unlist(lisUexv2_100_raw)
  
  lisUexv2_500_raw<-lapply(lisUn, function(x){ alpha_PDI(t(x$small4), x$res_abun, corrected=F)})
  lisUexv2_500_raw<-unlist(lisUexv2_500_raw)
  
}

summ_aPDI_raw<-data.frame(S10=lisUexv2_10_raw-lisUexv2_raw, 
                      S50=lisUexv2_50_raw-lisUexv2_raw, 
                      S100=lisUexv2_100_raw-lisUexv2_raw,
                      S500=lisUexv2_500_raw-lisUexv2_raw, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))


lisUn2<-mat2
posslm1 <- possibly(.f = wcfun, otherwise = NA)
if (T){
  
  lisUwc<-map2(lisUn, lisUn2, ~ posslm1(t(.y$preference), rep(1,nrow(.x$preference))))
  lisUwc<-sapply(lisUwc, function(x){ x[[1]]})
  
  lisUwc_10<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small1), .x$res_abun))
  lisUwc_10<-sapply(lisUwc_10, function(x){ x[[1]]})
  
  lisUwc_50<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small2), .x$res_abun))
  lisUwc_50<-sapply(lisUwc_50, function(x){ x[[1]]})
  
  lisUwc_100<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small3), .x$res_abun))
  lisUwc_100<-sapply(lisUwc_100, function(x){ x[[1]]})
  
  lisUwc_500<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small4), .x$res_abun))
  lisUwc_500<-sapply(lisUwc_500, function(x){ x[[1]]})
  
}


summ_wc<-data.frame(S10=lisUwc_10-lisUwc, 
                      S50=lisUwc_50-lisUwc, 
                      S100=lisUwc_100-lisUwc,
                      S500=lisUwc_500-lisUwc, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))







######################### 3. PLOTTING ######################################

png(filename="Figures/Exported/Figure3.png", width=6000, height=4000, res=600)

if (T){
  x1_xit<-jitter(rep(1, l5), factor=6)
  x2_xit<-jitter(rep(2, l5), factor=6)
  x3_xit<-jitter(rep(3, l5), factor=6)
  
  x_xit<-c(x1_xit,x2_xit,x3_xit)
  colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
  spar_col<-as.raster(matrix(colfunc(2000), nrow=1))
  spar_col<-rep(rep(spar_col, each=1), 3)
  
  par(mar= c(1,3,2,0),las=1, xpd=F)
  layout(matrix(c(20,16,16,16,16,15,
                  13, 1, 2, 3, 4,15,
                  13,17,17,17,17,15, 
                  13, 5, 6, 7, 8,15,
                  13,18,18,18,18,15,
                  13, 9,10,11,12,15,
                  19,14,14,14,14,15), ncol=6, byrow=T),
         widths=c(6,45/2,45/2,45/2,45/2,11), heights=c(6,40,6,40,6,40,10))
  #layout.show(20)
  
  
  ## aPDI'
  boxplot(S10~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), main="10 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI[,1], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="50 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI[,2], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),2])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),2])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="100 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI[,3], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),3])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),3])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), main="500 recorded interactions\n")
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI[,4], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),1])/l5,4),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),1])/l5,4),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  ## aPDI
  par(mar= c(1.5,3,1.5,0),las=1, xpd=F)
  boxplot(S10~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI_raw$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI_raw[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI_raw[,1], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI_raw$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI_raw[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI_raw[,2], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==5),2])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==15),2])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI_raw$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI_raw[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI_raw[,3], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==5),3])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==15),3])/l5,2),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI_raw$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_aPDI_raw, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_aPDI_raw[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_aPDI_raw[,4], ncol=3, byrow=F), 2, function(x){biasfun(0,x)}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==5),1])/l5,4),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==15),1])/l5,4),
         round(sum(is.na(summ_aPDI_raw)[which(summ_aPDI_raw$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  ## Wc
  par(mar= c(2,3,1,0),las=1, xpd=F)
  boxplot(S10~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  abline(h=0, lty=2)
  points(summ_wc$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_wc[,1], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_wc[,1], ncol=3, byrow=F), 2, function(x){biasfun(0,na.omit(x))}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),1])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),1])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),1])/l5,3))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  
  boxplot(S50~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_wc$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_wc[,2], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_wc[,2], ncol=3, byrow=F), 2, function(x){biasfun(0,na.omit(x))}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),2])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),2])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),2])/l5,3))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_wc$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_wc[,3], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_wc[,3], ncol=3, byrow=F), 2, function(x){biasfun(0,na.omit(x))}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),3])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),3])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),3])/l5,3))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S500~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_wc$S500~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S500~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  
  msqe<-round(apply(matrix(summ_wc[,4], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  
  biasc<-round(apply(matrix(summ_wc[,4], ncol=3, byrow=F), 2, function(x){biasfun(0,na.omit(x))}),3)
  biasc<-sprintf(biasc, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("Bias = "*.(biasc[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("Bias = "*.(biasc[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("Bias = "*.(biasc[3])), cex=0.7)
  
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),4])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),4])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),4])/l5,3))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.8, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.8, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.8, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  
  ### Y Label
  
  par(mar= c(0,0.5,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.5, labels="Difference between the observed and expected degree of generalisation", cex=1, srt=90)
  
  ### X Label
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.7, labels="Number of potential resources", cex=1)
  
  
  ###### Legend
  
  par(mar= c(4,2,4,0),las=1)
  plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialisation\nparameter",
       cex.main=1, ylim=c(0,1))
  legend_image <- as.raster(rev(colfunc(20)))
  grid.raster(legend_image, width=0.03, height = 0.8, x = unit(0.94, "npc"))
  par(xpd=T)
  abline(v=-0.2, lty=3)
  text(y=c(-0.02,1.02), x =0.7 , labels = c(0.1, "60"), cex=1)
  par(xpd=F)
  par(new=T, mar=c(42.5,4.8,6,0.5))
  
  x <- seq(-4, 4, length=1000)
  y1 <- dnorm(x,0,0.1)/dnorm(c(0),0,0.1)
  y1.5 <- dnorm(x,0,0.5)/dnorm(c(0),0,0.5)
  y2 <- dnorm(x,0,1)/dnorm(c(0),0,1)
  y2.5 <- dnorm(x,0,5)/dnorm(c(0),0,5)
  y3 <- dnorm(x,0,20)/dnorm(c(0),0,20)
  
  plot(x,y1, type = "l", lwd = 1, xlab = "", ylab = "", ylim=c(0,max(y1)),
       xaxt="n", yaxt="n",  cex=1.1, bty="l", col=colfunc(5)[5])
  
  par(new=T,mar=c(33.5,4.8,15,0.5))
  plot(x,y1.5, type = "l", lwd = 1, xlab = "", ylab = "", ylim=c(0,max(y1.5)),
       xaxt="n", yaxt="n", cex=1.1, bty="l", col=colfunc(5)[4])
  
  par(new=T,mar=c(24.25,4.8,24.25,0.5))
  plot(x,y2, type = "l", lwd = 1.5, xlab = "", ylab = "", ylim=c(0,max(y2)),
       xaxt="n", yaxt="n", cex=1.1, bty="l", col=colfunc(5)[3])
  
  par(new=T,mar=c(15,4.8,33.5,0.5))
  plot(x,y2.5, type = "l", lwd = 1, xlab = "", ylab = "", ylim=c(0,max(y2.5)),
       xaxt="n", yaxt="n",  cex=1.1, bty="l", col=colfunc(5)[2])
  
  par(new=T,mar=c(6,4.8,42.5,0.5))
  plot(x,y3, type = "l", lwd = 1, xlab = "", ylab = "", ylim=c(0,max(y3)),
       xaxt="n", yaxt="n",  cex=1.1, bty="l", col=colfunc(5)[[1]])
  
  par(new=F)
  
  ### aPDI' main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.5, labels=expression(alpha*"PDI'"), cex=2)
  
  ### aPDI main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.3, labels=expression(alpha*"PDI"), cex=2)
  
  ### Wc main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.3, labels=expression(italic("Wc'")), cex=2)
  
}

dev.off()

