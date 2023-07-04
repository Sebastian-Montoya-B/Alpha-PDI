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
source("Code/QNM.R")

## Generate the list of vectors with even (lisEv) and uneven (lisUn) resource 
## abundance distributions, using the quantitative niche model of 
## Fründ et al. (2016)

Nbee <- 1 # number of consumers. If > 1 it will generate matrices.
nsim<-1
MaUn<-NULL
MaEv<-NULL
#spen<-c(seq(1/10000000000, 50, length=200)) # Specialization parameter
#spen<-sfsmisc::lseq(0.0000001, 60, length=200)
spen<-seq(0.01, 50, length=200)
length(spen)
lisUn<-NULL
lisEv<-NULL
spelisUn<-NULL
spelisEv<-NULL
counter<-1
lisnam<-NULL
for (Nplant in c(5, 10, 50)){ # Number of potential resources
  
  for (spe in spen){
    
    for (i in 1:nsim){
      MaUn[[i]]<-gen_uneven2(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),maxsamp=1000, make="random" )
      #MaEv[[i]]<-gen_even2(Nbee,Nplant, spe,minsamp=c(10,50,80,100),maxsamp=1000, samp=T, make="random")
      
    }
    lisnam[[counter]]<-Nplant
    lisUn[[counter]]<-MaUn
    #lisEv[[counter]]<-MaEv
    counter<-counter+1
  }
  
}

## For each consumer in lisEv and lisUn there are eight vectors: 
##   (1) the resource abundance distribution ($res_abun)
##   (2) the true preferences ($preference)
##   (3) the current pattern of resource use ($current)
##   (4) the observed pattern of resource use in a case with 10^6 observations ($large)
##   (5) the observed pattern of resource use in a case with 10 observations ($small10)
##   (6) the observed pattern of resource use in a case with 50 observations ($small50)
##   (7) the observed pattern of resource use in a case with 80 observations ($small80)
##   (8) the observed pattern of resource use in a case with 100 observations ($small100)

######################### 2. CALCULATIONS ######################################
lisEv[1]
if (T){
  lisUexv2<-lapply(lisUn, function(x){
    lapply(x, function(x){
      alpha_PDI(t(x$large), x$res_abun)$corrected_aPDI})})
  lisUexv2<-unlist(lisUexv2)
  
  lisUexvs2_10<-lapply(lisUn, function(x){
    lapply(x, function(x){
      alpha_PDI(t(x$small10), x$res_abun)$corrected_aPDI})})
  lisUexvs2_10<-unlist(lisUexvs2_10)
  
  lisUexvs2_50<-lapply(lisUn, function(x){
    lapply(x, function(x){
      alpha_PDI(t(x$small50), x$res_abun)$corrected_aPDI})})
  lisUexvs2_50<-unlist(lisUexvs2_50)
  
  lisUexvs2_80<-lapply(lisUn, function(x){
    lapply(x, function(x){
      alpha_PDI(t(x$small80), x$res_abun)$corrected_aPDI})})
  lisUexvs2_80<-unlist(lisUexvs2_80)
  
  lisUexvs2_100<-lapply(lisUn, function(x){
    lapply(x, function(x){
      alpha_PDI(t(x$small100), x$res_abun)$corrected_aPDI})})
  lisUexvs2_100<-unlist(lisUexvs2_100)
}


l5<-length(lisUexv2)/3
l10<-length(lisUexv2)/3*2
l50<-length(lisUexv2)



summ_uneven<-data.frame(S10=lisUexvs2_10-lisUexv2, 
                      S50=lisUexvs2_50-lisUexv2, 
                      S80=lisUexvs2_80-lisUexv2,
                      S100=lisUexvs2_100-lisUexv2, res=c(rep(5, l5),rep(10, l5), rep(50, l5)))

x1_xit<-jitter(rep(1, l5), factor=6)
x2_xit<-jitter(rep(2, l5), factor=6)
x3_xit<-jitter(rep(3, l5), factor=6)

x_xit<-c(x1_xit,x2_xit,x3_xit)
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
spar_col<-as.raster(matrix(colfunc(length(spen)), nrow=1))
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
