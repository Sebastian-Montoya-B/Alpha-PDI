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
nsim<-2
MaUn<-NULL
MaEv<-NULL
spen<-c(seq(0.1, 60, by=0.1)) # Specialization parameter
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
      MaUn[[i]]<-gen_uneven2(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),maxsamp=1000000, make="random" )
      MaEv[[i]]<-gen_even2(Nbee,Nplant, spe,minsamp=c(10,50,80,100),maxsamp=1000000, samp=T, make="random")
      
    }
    lisnam[[counter]]<-Nplant
    lisUn[[counter]]<-MaUn
    lisEv[[counter]]<-MaEv
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


lisEexv<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$large), rep(1,NROW(x$large)))$corrected_aPDI})})
lisEexv<-unlist(lisEexv)

lisEexvs10<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small10), rep(1,NROW(x$small10)))$corrected_aPDI})})
lisEexvs10<-unlist(lisEexvs10)

lisEexvs50<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small50), rep(1,NROW(x$small50)))$corrected_aPDI})})
lisEexvs50<-unlist(lisEexvs50)

lisEexvs80<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small80), rep(1,NROW(x$small80)))$corrected_aPDI})})
lisEexvs80<-unlist(lisEexvs80)

lisEexvs100<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$small100), rep(1,NROW(x$small100)))$corrected_aPDI})})
lisEexvs100<-unlist(lisEexvs100)

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

l5<-length(lisEexv)/3
l10<-length(lisEexv)/3*2
l50<-length(lisEexv)




######################### 3. PLOTTING ##########################################


colw<-c("#00ceff", "#078ab5","#004c6d")
png(filename="Figures/Exported/Figure3.png", width=4000, height=4600, res=600)
#svg(filename="Figures/Exported/Figure3.svg", width=8, height=7)


layout(matrix(c(17,17,1,1,
                17,17,2,3,
                4,5,9,10,
                4,6,11,12,
                4,7,13,14,
                4,8,15,16), ncol=4, byrow=T),
       widths=c(13,7,40,40), heights=c(5,5,rep(90/4,4)))
layout.show(17)

#1
par(mar= c(0,0,0,0),las=1)
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels="Resource abundance distribution", cex=1.6)

#2
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.575, y=0.55, labels="Even",  cex=1.6)

#3
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.575, y=0.55, labels="Uneven", cex=1.6)

#4
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels="Sampling intensity\n(Number of observations per consumer)",srt=90, cex=1.6)


#5
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.7, labels="10",  cex=1.6)

#6
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.7, labels="50",  cex=1.6)
#7
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.7, labels="80",  cex=1.6)
#8
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.7, labels="100",  cex=1.6)

#9
par(mar= c(3.4,4,0,1),las=1,cex.lab=1, mgp=c(2.3,0.5,0))

plot(x=lisEexv[1:l5],y=lisEexvs10[1:l5], col=colw[3], pch=21,
     bg=alpha(colw[3],0.5), cex=0.7, 
     ylab="Observed generality",
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisEexv[(l5+1):l10],y=lisEexvs10[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisEexv[(l10+1):l50],y=lisEexvs10[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisEexv,lisEexvs10)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#10
plot(x=lisUexv2[1:l5],y=lisUexvs2_10[1:l5], col=colw[3], pch=21, 
     bg=alpha(colw[3],0.5), cex=0.7, ylab="", 
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("","","","",""))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisUexv2[(l5+1):l10],y=lisUexvs2_10[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisUexv2[(l10+1):l50],y=lisUexvs2_10[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisUexv2,lisUexvs2_10)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#11

plot(x=lisEexv[1:l5],y=lisEexvs50[1:l5], col=colw[3], pch=21,
     bg=alpha(colw[3],0.5), cex=0.7, 
     ylab="Observed generality",
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisEexv[(l5+1):l10],y=lisEexvs50[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisEexv[(l10+1):l50],y=lisEexvs50[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)

r<-round(cor.test(lisEexv,lisEexvs50)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#12

plot(x=lisUexv2[1:l5],y=lisUexvs2_50[1:l5], col=colw[3], pch=21, 
     bg=alpha(colw[3],0.5), cex=0.7, ylab="", 
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("","","","",""))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisUexv2[(l5+1):l10],y=lisUexvs2_50[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisUexv2[(l10+1):l50],y=lisUexvs2_50[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisUexv2,lisUexvs2_50)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#13

plot(x=lisEexv[1:l5],y=lisEexvs80[1:l5], col=colw[3], pch=21,
     bg=alpha(colw[3],0.5), cex=0.7, 
     ylab="Observed generality",
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisEexv[(l5+1):l10],y=lisEexvs80[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisEexv[(l10+1):l50],y=lisEexvs80[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisEexv,lisEexvs80)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#14

plot(x=lisUexv2[1:l5],y=lisUexvs2_80[1:l5], col=colw[3], pch=21, 
     bg=alpha(colw[3],0.5), cex=0.7, ylab="", 
     xlab="", xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("","","","",""))
axis(1, at=seq(0,1, length=5), labels=c("","","","",""))
points(x=lisUexv2[(l5+1):l10],y=lisUexvs2_80[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisUexv2[(l10+1):l50],y=lisUexvs2_80[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisUexv2,lisUexvs2_80)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#15

plot(x=lisEexv[1:l5],y=lisEexvs100[1:l5], col=colw[3], pch=21,
     bg=alpha(colw[3],0.5), cex=0.7, 
     ylab="Observed generality",
     xlab=expression("Expected generality"),
     xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
axis(1, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
points(x=lisEexv[(l5+1):l10],y=lisEexvs100[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisEexv[(l10+1):l50],y=lisEexvs100[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisEexv,lisEexvs100)$estimate, digits=2)
if (r==1){
  r<-0.99
}
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#16

plot(x=lisUexv2[1:l5],y=lisUexvs2_100[1:l5], col=colw[3], pch=21, 
     bg=alpha(colw[3],0.5), cex=0.7, ylab="", 
     xlab=expression("Expected generality"), 
     xlim=c(0,1), ylim=c(0,1),xaxp=c(0,1,2), yaxp=c(0,1,2),
     yaxt="n", xaxt="n")
axis(2, at=seq(0,1, length=5), labels=c("","","","",""))
axis(1, at=seq(0,1, length=5), labels=c("0.0","","0.5","","1.0"))
points(x=lisUexv2[(l5+1):l10],y=lisUexvs2_100[(l5+1):l10],col=colw[2],
       pch=23, bg=alpha(colw[2],0.5),cex=0.7)
points(x=lisUexv2[(l10+1):l50],y=lisUexvs2_100[(l10+1):l50],col=colw[1],
       pch=22,bg=alpha(colw[1],0.5),cex=0.7)
abline(coef = c(0,1), lwd=1.5)
r<-round(cor.test(lisUexv2,lisUexvs2_100)$estimate, digits=2)
r2<-round(r^2, digits=2)
text(x=0.15, y=0.95, labels=bquote(rho == .(r)), cex=1.2)
text(x=0.15, y=0.8, labels=bquote("r"^2 == .(r2)), cex=1.2)

#17
par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
legend(x=0, y=1, legend=c("5", "10", "50"), pch=c(21,23,22),
       pt.bg=c(alpha(colw[3],0.5),alpha(colw[2],0.5),alpha(colw[1],0.5)),
       col=c(colw[3],colw[2],colw[1]), title="Potential resources",
       horiz=T, pt.cex=2, cex=1.3, bty="n")

dev.off()
