################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure 2.


######################### 1. SETTINGS ##########################################


## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.
if(!require(zCompositions)){
  install.packages("zCompositions")
  library(zCompositions)
}

if(!require(bipartite)){
  install.packages("bipartite")
  library(bipartite)
}

if(!require(scales)){
  install.packages("scales")
  library(scales)
}

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/wcfun.R")
source("Code/QNM.R")

## Load a list of matrices and vectors generated using the quantitative niche
## model of Fründ et al. (2016)

Nbee <- 5
nsim<-1
MaUn<-NULL
MaEv<-NULL
spen<-c(seq(0.1, 60, by=0.1)) #Specialization parameter
length(spen)
lisUn<-NULL
lisEv<-NULL
spelisUn<-NULL
spelisEv<-NULL
counter<-1
lisnam<-NULL
for (Nplant in c(5, 10, 50)){
  
  for (spe in spen){
    
    for (i in 1:nsim){
      MaUn[[i]]<-gen_uneven(Nbee,Nplant, spe, samp=F, make="random" )
      
    }
    lisnam[[counter]]<-Nplant
    lisUn[[counter]]<-MaUn
    counter<-counter+1
  }

}

mat1<-lisUn

## For each consumer in mat1 there are four vectors/matrices: 
##   (1) the consumer abundance distribution ($con_abun)
##   (2) the resource abundance distribution ($res_abun)
##   (3) the true preferences ($preference) 
##   (4) the current resource use pattern ($current).


######################### 2. CALCULATIONS ######################################


## Calculating αPDI and other indices
exv<-NULL
obv<-NULL
lisExv<-NULL
lisObv<-NULL
for (i in 1:length(mat1)){
  for (j in 1:length(mat1[[1]])){
    
    ap<-alpha_PDI(t((mat1[[i]][[j]]$preference)), rep(1,NROW(mat1[[i]][[j]]$preference)),corrected=F)
    at<-alpha_PDI(t((mat1[[i]][[j]]$current)), mat1[[i]][[j]]$res_abun, corrected=F)
    gp<-genfun(t((mat1[[i]][[j]]$preference*100)), rep(1,NROW(mat1[[i]][[j]]$preference)))
    gt<-genfun(t(mat1[[i]][[j]]$current*100), mat1[[i]][[j]]$res_abun)
    
    
    exv[[j]]<-cbind(gp, aPDI=ap)
    obv[[j]]<-cbind(gt, aPDI=at)

  }
  lisExv[[i]]<-exv
  lisObv[[i]]<-obv
}


lisExv<-dplyr::bind_rows(lisExv)
lisObv<-dplyr::bind_rows(lisObv)
length(lisObv)

## Calculating Wc Pierotti et al. (2017)
wcp<-NULL
wct<-NULL
lisObvwc<-NULL
lisExvwc<-NULL
exvwc<-NULL
obvwc<-NULL
for (i in 1:length(mat1)){
  for (j in 1:length(mat1[[1]])){
    wcp<-wcfun(t(mat1[[i]][[j]]$preference), rep(1,NROW(mat1[[i]][[j]]$preference)))
    wct<-wcfun(t(mat1[[i]][[j]]$current), mat1[[i]][[j]]$res_abun)
    
    exvwc[[j]]<-wcp
    obvwc[[j]]<-wct
  }
  lisExvwc[[i]]<-exvwc
  lisObvwc[[i]]<-obvwc
}

lisExvwc<-unlist(lisExvwc)
lisObvwc<-unlist(lisObvwc)
## Depending on your data, some zeros may not be replaced using zCompositions.
## Error in if (any(X2[i, z] > colmins[z])) { : 
## missing value where TRUE/FALSE needed
## See error in https://stats.stackexchange.com/questions/477663/error-with-the-geometric-bayesian-multiplicative-replacement-of-count-zeros-with


######################### 3. PLOTTING ##########################################


colw<-c("#00ceff", "#078ab5","#004c6d")

svg(filename="Figures/Exported/Figure2.svg", width=8, height=9)
par(mar= c(4,2,2,1),las=1)
layout(matrix(seq(1,9), ncol=3, byrow=T))
layout.show(9)
l5<-length(lisObv$aPDI)/3
l10<-length(lisObv$aPDI)/3*2
l50<-length(lisObv$aPDI)

###aPDI
plot(y=lisObv$aPDI[1:l5],x=lisExv$aPDI[1:l5], xlim=c(0,1), ylim=c(0,1), pch=21, col=colw[3], bg=alpha(colw[3],0.4), cex=0.7, 
     ylab="Estimated value", xlab="", xaxt="n")
axis(1, labels=F)
points(y=lisObv$aPDI[(l5+1):l10],x=lisExv$aPDI[(l5+1):l10], xlim=c(0,1), ylim=c(0,1), col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$aPDI[(l10+1):l50],x=lisExv$aPDI[(l10+1):l50],xlim=c(0,1), ylim=c(0,1), col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.09,y=0.98, label=expression(alpha*italic(PDI)), cex=1.3)

###Bs
plot(y=lisObv$Bs[1:l5],x=lisExv$Bs[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7, 
     ylab="Estimated value", xlab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n")
axis(2, labels=F)
axis(1, labels=F)
points(y=lisObv$Bs[(l5+1):l10],x=lisExv$Bs[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$Bs[(l10+1):l50],x=lisExv$Bs[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(B[s])), cex=1.3)

###B
plot(y=lisObv$`B'`[1:l5],x=lisExv$`B'`[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7, 
     ylab="Estimated value", xlab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n")
axis(2, labels=F)
axis(1, labels=F)
points(y=lisObv$`B'`[l5+1:l10],x=lisExv$`B'`[l5+1:l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$`B'`[l10+1:l50],x=lisExv$`B'`[l10+1:l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic("B'")), cex=1.3)

###W
plot(y=lisObv$W[1:l5],x=lisExv$W[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7,
     ylab="Estimated value", xlab="", xlim=c(0,1), ylim=c(0,1), xaxt="n")
axis(1, labels=F)
points(y=lisObv$W[(l5+1):l10],x=lisExv$W[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$W[(l10+1):l50],x=lisExv$W[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(W)), cex=1.3)

###PS
plot(y=lisObv$PS[1:l5],x=lisExv$PS[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7,
     ylab="Estimated value", xlab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n")
axis(2, labels=F)
axis(1, labels=F)
points(y=lisObv$PS[(l5+1):l10],x=lisExv$PS[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$PS[(l10+1):l50],x=lisExv$PS[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(PS)), cex=1.3)

###FT
plot(y=lisObv$FT[1:l5],x=lisExv$FT[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7,
     ylab="Estimated value", xlab="", xlim=c(0,1), ylim=c(0,1), yaxt="n", xaxt="n")
axis(2, labels=F)
axis(1, labels=F)
points(y=lisObv$FT[(l5+1):l10],x=lisExv$FT[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$FT[(l10+1):l50],x=lisExv$FT[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(FT)), cex=1.3)

###1-d
plot(y=lisObv$`1-d'`[1:l5],x=lisExv$`1-d'`[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7,
     ylab="Estimated value", xlab="Expected value", xlim=c(0,1), ylim=c(0,1))
points(y=lisObv$`1-d'`[(l5+1):l10],x=lisExv$`1-d'`[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$`1-d'`[(l10+1):l50],x=lisExv$`1-d'`[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.09,y=0.98, label= expression(1 - italic("d'")), cex=1.3)

###gen
plot(y=lisObv$gen[1:l5],x=lisExv$gen[1:l5], col=colw[3], pch=21, bg=alpha(colw[3],0.4), cex=0.7,
     ylab="Estimated value", xlab="Expected value", xlim=c(0,1), ylim=c(0,1), yaxt="n")
axis(2, labels=F)
points(y=lisObv$gen[(l5+1):l10],x=lisExv$gen[(l5+1):l10],col=colw[2],pch=23, bg=alpha(colw[2],0.4),cex=0.7)
points(y=lisObv$gen[(l10+1):l50],x=lisExv$gen[(l10+1):l50],col=colw[1], pch=22,bg=alpha(colw[1],0.4),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(gen)), cex=1.3)

###Wc
plot(y=lisObvwc[1:l5],x=lisExvwc[1:l5], xlim=c(0,1), ylim=c(0,1), pch=21, col=colw[3], bg=alpha(colw[3],0.5), cex=0.7, 
     ylab="Estimated value", xlab="Expected value", yaxt="n")
axis(2, labels=F)
points(y=lisObvwc[(l5+1):l10],x=lisExvwc[(l5+1):l10], xlim=c(0,1), ylim=c(0,1), col=colw[2],pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(y=lisObvwc[(l10+1):l50],x=lisExvwc[(l10+1):l50],xlim=c(0,1), ylim=c(0,1), col=colw[1], pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
text(x=0.07,y=0.98, label=expression(italic(Wc)), cex=1.3)

dev.off()
