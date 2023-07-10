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

## Source the functions.
source("Code/alpha_PDI.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")
mat2<-readRDS("Data/sim_data2.rds")


######################### 2. CALCULATIONS ######################################
msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}


lisUn<-mat1
if (T){

  lisUexv2<-lapply(lisUn, function(x){ alpha_PDI(t(x$large), x$res_abun)$corrected_aPDI})
  lisUexv2<-unlist(lisUexv2)
  
  lisUexv2_10<-lapply(lisUn, function(x){ alpha_PDI(t(x$small10), x$res_abun)$corrected_aPDI})
  lisUexv2_10<-unlist(lisUexv2_10)
  
  lisUexv2_50<-lapply(lisUn, function(x){ alpha_PDI(t(x$small50), x$res_abun)$corrected_aPDI})
  lisUexv2_50<-unlist(lisUexv2_50)
  
  lisUexv2_80<-lapply(lisUn, function(x){ alpha_PDI(t(x$small80), x$res_abun)$corrected_aPDI})
  lisUexv2_80<-unlist(lisUexv2_80)
  
  lisUexv2_100<-lapply(lisUn, function(x){ alpha_PDI(t(x$small100), x$res_abun)$corrected_aPDI})
  lisUexv2_100<-unlist(lisUexv2_100)

}


l5<-length(lisUexv2)/3
l10<-length(lisUexv2)/3*2
l50<-length(lisUexv2)




summ_aPDI<-data.frame(S10=lisUexv2_10-lisUexv2, 
                      S50=lisUexv2_50-lisUexv2, 
                      S80=lisUexv2_80-lisUexv2,
                      S100=lisUexv2_100-lisUexv2, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))


lisUn2<-mat2
posslm1 <- possibly(.f = wcfun, otherwise = NA)
if (T){
  
  lisUwc<-map2(lisUn, lisUn2, ~ posslm1(t(.y$large), .x$res_abun))
  lisUwc<-sapply(lisUwc, function(x){ x[[1]]})
  
  lisUwc_10<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small10), .x$res_abun))
  lisUwc_10<-sapply(lisUwc_10, function(x){ x[[1]]})
  
  lisUwc_50<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small50), .x$res_abun))
  lisUwc_50<-sapply(lisUwc_50, function(x){ x[[1]]})
  
  lisUwc_80<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small80), .x$res_abun))
  lisUwc_80<-sapply(lisUwc_80, function(x){ x[[1]]})
  
  lisUwc_100<-map2(lisUn, lisUn2, ~ posslm1(t(.y$small100), .x$res_abun))
  lisUwc_100<-sapply(lisUwc_100, function(x){ x[[1]]})
  
}


summ_wc<-data.frame(S10=lisUwc_10-lisUwc, 
                      S50=lisUwc_50-lisUwc, 
                      S80=lisUwc_80-lisUwc,
                      S100=lisUwc_100-lisUwc, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))




x1_xit<-jitter(rep(1, l5), factor=6)
x2_xit<-jitter(rep(2, l5), factor=6)
x3_xit<-jitter(rep(3, l5), factor=6)

x_xit<-c(x1_xit,x2_xit,x3_xit)
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
spar_col<-as.raster(matrix(colfunc(2000), nrow=1))
spar_col<-rep(rep(spar_col, each=Nbee), 3)

png(filename="Figures/Exported/Figure3.1.png", width=6000, height=3000, res=600)
if (T){
  par(mar= c(2,3,1,0),las=1, xpd=F)
  layout(matrix(c(15,12,12,12,12,11,
                   9, 1, 2, 3, 4,11,
                   9,13,13,13,13,11, 
                   9, 5, 6, 7, 8,11,
                  14,10,10,10,10,11), ncol=6, byrow=T),
         widths=c(6,45/2,45/2,45/2,45/2,11), heights=c(6,40,6,40,10))
  #layout.show(15)
  
  
  ## aPDI
  boxplot(S10~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S50~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S50~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI[,2], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),2])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),2])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S80~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S80~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S80~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI[,3], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),3])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),3])/l5,2),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),3])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S100~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_aPDI$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_aPDI, ylim=c(-1,1), xaxt="n", yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI[,4], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==5),1])/l5,4),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==15),1])/l5,4),
         round(sum(is.na(summ_aPDI)[which(summ_aPDI$res==55),1])/l5,4))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  ## Wc
  boxplot(S10~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  abline(h=0, lty=2)
  points(summ_wc$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S10~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_wc[,1], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),1])/l5,2),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),1])/l5,2),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  
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
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),2])/l5,2),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),2])/l5,2),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),2])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  boxplot(S80~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_wc$S80~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S80~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_wc[,3], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),1])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),1])/l5,3),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),1])/l5,3))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
 
  
  boxplot(S100~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0))
  axis(1, at=3, labels="55")
  axis(2, at=seq(-1,1, length=5), labels=F)
  abline(h=0, lty=2)
  points(summ_wc$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  par(new=T)
  boxplot(S100~res, data=summ_wc, ylim=c(-1,1), xlim=c(0.5,3.5), yaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4))
  par(new=F)
  msqe<-round(apply(matrix(summ_wc[,4], ncol=3, byrow=F), 2, function (x) {msqerrorfun(na.omit(x))}),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_wc)[which(summ_wc$res==5),1])/l5,4),
         round(sum(is.na(summ_wc)[which(summ_wc$res==15),1])/l5,4),
         round(sum(is.na(summ_wc)[which(summ_wc$res==55),1])/l5,4))
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
  text(x=0.535, y=0.7, labels="Number of potential resources", cex=1)
  
  
  ###### Legend
  
  par(mar= c(4,2,4,0),las=1)
  plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialization\nparameter",
       cex.main=1, ylim=c(0,1))
  legend_image <- as.raster(rev(colfunc(20)))
  grid.raster(legend_image, width=0.03, height = 0.8, x = unit(0.94, "npc"))
  par(xpd=T)
  abline(v=-0.2, lty=3)
  text(y=c(-0.03,1.03), x =0.7 , labels = c(0.01, "3000"), cex=1)
  par(xpd=F)
  
  ### aPDI main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.4, labels=expression(alpha*italic(PDI)*"'"), cex=2)
  
  ### Wc main
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.4, labels=expression(italic(Wc)), cex=2)
  
}
dev.off()

