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

## Source the functions.
source("Code/alpha_PDI.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")
mat2<-readRDS("Data/sim_data2.rds")


######################### 2. CALCULATIONS ######################################
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



summ_uneven<-data.frame(S10=lisUexv2_10-lisUexv2, 
                      S50=lisUexv2_50-lisUexv2, 
                      S80=lisUexv2_80-lisUexv2,
                      S100=lisUexv2_100-lisUexv2, res=c(rep(5, l5),rep(10, l5), rep(50, l5)))

x1_xit<-jitter(rep(1, l5), factor=6)
x2_xit<-jitter(rep(2, l5), factor=6)
x3_xit<-jitter(rep(3, l5), factor=6)

x_xit<-c(x1_xit,x2_xit,x3_xit)
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
spar_col<-as.raster(matrix(colfunc(2000), nrow=1))
spar_col<-rep(rep(spar_col, each=Nbee), 3)

layout(matrix(c(1,2,3,4,
                5,6,7,8), ncol=4, byrow=T))
layout.show(8)



boxplot(S10~res, data=summ_uneven, ylim=c(-1,1))
points(summ_uneven$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S50~res, data=summ_uneven, ylim=c(-1,1))
points(summ_uneven$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S80~res, data=summ_uneven, ylim=c(-1,1))
points(summ_uneven$S80~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S100~res, data=summ_uneven, ylim=c(-1,1))
points(summ_uneven$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)


### Bprime


summ_uneven_Bprime<-data.frame(S10=lisUexvs2_10_o$`B'`-lisUexv2_o$`B'`, 
                               S50=lisUexvs2_50_o$`B'`-lisUexv2_o$`B'`, 
                               S80=lisUexvs2_80_o$`B'`-lisUexv2_o$`B'`,
                               S100=lisUexvs2_100_o$`B'`-lisUexv2_o$`B'`, res=c(rep(5, l5),rep(10, l5), rep(50, l5)))

boxplot(S10~res, data=summ_uneven_Bprime, ylim=c(-1,1))
points(summ_uneven_Bprime$S10~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S50~res, data=summ_uneven_Bprime, ylim=c(-1,1))
points(summ_uneven_Bprime$S50~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S80~res, data=summ_uneven_Bprime, ylim=c(-1,1))
points(summ_uneven_Bprime$S80~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)

boxplot(S100~res, data=summ_uneven_Bprime, ylim=c(-1,1))
points(summ_uneven_Bprime$S100~x_xit, bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
