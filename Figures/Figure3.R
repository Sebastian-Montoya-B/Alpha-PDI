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

## Load the list of vectors generated using even resource abundance distributions.
lisEv <- readRDS("Data/samp_vectors_even.RDS") 

## Load the list of vectors generated using uneven resource abundance distributions.
lisUn <- readRDS("Data/samp_vectors_uneven.RDS") 

## For each consumer in lisEv and lisUn there are five vectors: 
##   (1) the resource abundance distribution ($res_abun)
##   (2) the true preferences ($preference)
##   (3) the current pattern of resource use ($current)
##   (4) the observed pattern of resource use in a case with 10^6 observations ($large)
##   (5) the observed pattern of resource use in a case with 100 observations ($small)


######################### 2. CALCULATIONS ######################################


lisEexv<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$large), rep(1,NROW(x$large)))$corrected_aPDI})})
lisEexv<-unlist(lisEexv)

lisEexvs<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small), rep(1,NROW(x$small)))$corrected_aPDI})})
lisEexvs<-unlist(lisEexvs)

lisUexv2<-lapply(lisUn, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$large), x$res_abun)$corrected_aPDI})})
lisUexv2<-unlist(lisUexv2)

lisUexvs2<-lapply(lisUn, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small), x$res_abun)$corrected_aPDI})})
lisUexvs2<-unlist(lisUexvs2)

l5<-length(lisEexv)/3
l10<-length(lisEexv)/3*2
l50<-length(lisEexv)

## Calculate mean bias

## Sampling bias when resource abundance distribution is even
mbias<-mean(100*(lisEexv-lisEexvs)/(max(lisEexv)-min(lisEexv)))#
sdbias<-sd(100*abs(lisEexv-lisEexvs)/(max(lisEexv)-min(lisEexv)))#
cor.test(lisEexv,lisEexvs)#0.99 

## Sampling bias when resource abundance distribution is uneven
mbias<-mean(100*(lisUexv2-lisUexvs2)/(max(lisUexv2)-min(lisUexv2)))#
sdbias<-sd(100*abs(lisUexv2-lisUexvs2)/(max(lisUexv2)-min(lisUexv2)))#
cor.test(lisUexv2,lisUexvs2)#0.98  


######################### 3. PLOTTING ##########################################


colw<-c("#00ceff", "#078ab5","#004c6d")

svg(filename="Figures/Figure3.svg", width=8, height=7)
x11()

par(mar= c(4,2,2,1),las=1)
layout(matrix(c(1,2,3,4), ncol=2, byrow=T))
layout.show(4)
plot(x=lisEexv[1:l5],y=lisEexvs[1:l5], col=colw[3], pch=21,
     bg=alpha(colw[3],0.5), cex=0.7, ylab="Observed value",
     xlab="Expected value", xlim=c(0,1), ylim=c(0,1))
points(x=lisEexv[(l5+1):l10],y=lisEexvs[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisEexv[(l10+1):l50],y=lisEexvs[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)

plot(x=lisUexv2[1:l5],y=lisUexvs2[1:l5], col=colw[3], pch=21, 
     bg=alpha(colw[3],0.5), cex=0.7, ylab="Observed value", 
     xlab="Expected value", xlim=c(0,1), ylim=c(0,1))
points(x=lisUexv2[(l5+1):l10],y=lisUexvs2[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisUexv2[(l10+1):l50],y=lisUexvs2[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)

dev.off()