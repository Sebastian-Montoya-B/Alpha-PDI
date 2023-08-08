################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure S6.


######################### 1. SETTINGS ##########################################

######################### 1. SETTINGS ##########################################


## Clean the environment.
rm(list= ls())

## Source the functions.
source("Code/alpha_PDI.R")


######################### 2. CALCULATIONS ######################################
p<-sort(c(0.00001*10^(0:10),sfsmisc::lseq(1.2,50,8), seq(0.05,0.9,length=8)))


r<-seq(0,10)
Pr<-NULL
for (i in 1:length(r)){
  Pr[[i]]<-(1-r[[i]]/length(r))^p
}


Pr
pr<-do.call(rbind, Pr)
pr

lin<-alpha_PDI(t(pr), rep(1,nrow(pr)), corrected=F)
lin

lin2<-alpha_PDI(t(pr), rep(1,nrow(pr)), corrected=F, m=5)
lin2

lin3<-alpha_PDI(t(pr), rep(1,nrow(pr)), corrected=F, m=10)
lin3







######################### 3. PLOTTING ######################################

colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
#licol<-colfunc(5)

toplo<-data.frame(p,lin)
toplo
mycol<-colfunc(nrow(toplo))


png(filename="Figures/Exported/FigureS6.png", width=5000, height=1300, res=600)

par(mar=c(3,3,2,1))
layout(matrix(c(5,1,2,3,
                6,4,4,4), byrow=T, ncol=4), heights=c(90,10), widths=c(5,30,30,30))


#aPDI
par(las=1)
plot(p, lin, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic(m)*" = 1"),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))
#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=toplo$lin[i], x1=toplo$p[i+1], y1=toplo$lin[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#aPDI l=5
par(las=1)
plot(p, lin2, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic(m)*" = 5"),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))
#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2[i], x1=toplo$p[i+1], y1=lin2[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)




#aPDI l=10
par(las=1)
plot(p, lin3, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic(m)*" = 10"),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))
#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin3[i], x1=toplo$p[i+1], y1=lin3[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)

### X Label

par(mar= c(0,0,0,0),las=1)
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.53, y=0.7, labels=expression("Specialization parameter "*"("*rho*")"), cex=1.3) 

### Y Label

par(mar= c(0,0.5,0,0),las=1)
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.53, y=0.53, labels=expression(alpha*italic(PDI)), cex=1.3, srt=90)

dev.off()