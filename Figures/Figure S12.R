################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure S12.


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

## Load the data
lisUn<-readRDS("Data/sim_data.rds")

######################### 2. CALCULATIONS ######################################

if (T){
  
  lisUexv2_m1<-lapply(lisUn, function(x){ alpha_PDI(t(x$preference), rep(1,nrow(x$preference)), corrected=F)})
  lisUexv2_m1<-unlist(lisUexv2_m1)
  
  lisUexv2_curr_m1<-lapply(lisUn, function(x){ alpha_PDI(t(x$current), x$res_abun, corrected=F)})
  lisUexv2_curr_m1<-unlist(lisUexv2_curr_m1)
  
}

l5<-length(lisUexv2_m1)/3
l10<-length(lisUexv2_m1)/3*2
l50<-length(lisUexv2_m1)


summ_aPDI_m1<-data.frame(dif=lisUexv2_curr_m1-lisUexv2_m1, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

if (T){
  
  lisUexv2_m5<-lapply(lisUn, function(x){ alpha_PDI(t(x$preference), rep(1,nrow(x$preference)), corrected=F, m=5)})
  lisUexv2_m5<-unlist(lisUexv2_m5)
  
  lisUexv2_curr_m5<-lapply(lisUn, function(x){ alpha_PDI(t(x$current), x$res_abun, corrected=F, m=5)})
  lisUexv2_curr_m5<-unlist(lisUexv2_curr_m5)
}

summ_aPDI_m5<-data.frame(diff=lisUexv2_curr_m5-lisUexv2_m5,
                         res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

if (T){
  
  lisUexv2_m10<-lapply(lisUn, function(x){ alpha_PDI(t(x$preference), rep(1,nrow(x$preference)), corrected=F, m=10)})
  lisUexv2_m10<-unlist(lisUexv2_m10)
  
  lisUexv2_curr_m10<-lapply(lisUn, function(x){ alpha_PDI(t(x$current), x$res_abun, corrected=F, m=10)})
  lisUexv2_curr_m10<-unlist(lisUexv2_curr_m10)
  
}

summ_aPDI_m10<-data.frame(dif=lisUexv2_curr_m10-lisUexv2_m10, res=c(rep(5, l5),rep(15, l5), rep(55, l5)))

msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}


######################### 3. PLOTTING ##########################################


x1_xit<-jitter(rep(1, l5), factor=6)
x2_xit<-jitter(rep(2, l5), factor=6)
x3_xit<-jitter(rep(3, l5), factor=6)


x_xit<-c(x1_xit,x2_xit,x3_xit)
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
spar_col<-as.raster(matrix(colfunc(2000), nrow=1))
spar_col<-rep(rep(spar_col, each=1), 3)


png(filename="Figures/Exported/FigureS12.png", width=5000, height=1300, res=600)

if (T) {
  layout(matrix(c(1,2,3), ncol=3))
  
  par(las=1, mar=c(4,5,2,1))
  boxplot(summ_aPDI_m1$dif~res, ylim=c(-1,1), pch=8, xaxt="n", main=expression(italic("m")*" = 1"), col=alpha("gray",0), 
          ylab=expression("   Difference between the observed\nand expected degree of generalization"), xlab="")
  points(summ_aPDI_m1$dif~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  axis(1, at=c(1,2,3), labels=c(5,15,55))
  par(new=T)
  boxplot(summ_aPDI_m1$dif~res, data=summ_aPDI_m1, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4), ylab="", xlab="")
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI_m1[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  
  boxplot(summ_aPDI_m5$dif~res, ylim=c(-1,1), pch=8, xaxt="n", main=expression(italic("m")*" = 5"), col=alpha("gray",0), ylab="", xlab="Number of potential resources")
  points(summ_aPDI_m5$dif~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  axis(1, at=c(1,2,3), labels=c(5,15,55))
  par(new=T)
  boxplot(summ_aPDI_m5$dif~res, data=summ_aPDI_m5, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4), ylab="", xlab="")
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI_m5[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI_m5)[which(summ_aPDI_m5$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI_m5)[which(summ_aPDI_m5$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI_m5)[which(summ_aPDI_m5$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
  
  
  boxplot(summ_aPDI_m1$dif~res, ylim=c(-1,1), pch=8, xaxt="n", main=expression(italic("m")*" = 10"), col=alpha("gray",0), ylab="", xlab="")
  points(summ_aPDI_m1$dif~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  axis(1, at=c(1,2,3), labels=c(5,15,55))
  par(new=T)
  boxplot(summ_aPDI_m1$dif~res, data=summ_aPDI_m1, ylim=c(-1,1), xaxt="n", pch=8, col=alpha("gray",0), border=alpha("black",0.4), ylab="", xlab="")
  par(new=F)
  msqe<-round(apply(matrix(summ_aPDI_m1[,1], ncol=3, byrow=F), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  Nas<-c(round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==5),1])/l5,2),
         round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==15),1])/l5,2),
         round(sum(is.na(summ_aPDI_m1)[which(summ_aPDI_m1$res==55),1])/l5,2))
  Nas<-sprintf(Nas, fmt='%#.3f')
  text(x=1, y=0.9, labels=bquote("NA = "*.(Nas[1])), cex=0.7)
  text(x=2, y=0.9, labels=bquote("NA = "*.(Nas[2])), cex=0.7)
  text(x=3, y=0.9, labels=bquote("NA = "*.(Nas[3])), cex=0.7)
  
}

dev.off()
