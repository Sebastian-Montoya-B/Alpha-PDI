################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure S2.


######################### 1. SETTINGS ##########################################


## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.
if(!require(sfsmisc)){
  install.packages("sfsmisc")
  library(sfsmisc)
}

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/wcfun.R")

######################### 2. CALCULATIONS ######################################

p<-c(0.00001,0.3,1,3,100000)
x<-seq(0,1,length=20)


y<-NULL
for (i in 1:length(x)){
  y[[i]]<-(1-x[i])^p
}
y

mres<-do.call(rbind, y)

mres<-cbind(mres, x)
mres

p2<-1
x2<-seq(0,1,length=10)

y2<-NULL
for (i in 1:length(x2)){
  y2[[i]]<-(1-x2[i])^p2
}
y2
mres2<-do.call(rbind, y2)
######################### 3. PLOTTING ##########################################
png(filename="Figures/Exported/Figure S2l21.png", width=4000, height=1590, res=600)
par(las=1)
layout(matrix(c(1,2,3), ncol=3))
bar<-barplot(rep(1,10), xaxt="n", yaxt="n",  ylim=c(0,1.2), 
     main="Extreme generalist", xlab="", ylab="Preferences", col="gray")
lines(x= bar, y=rep(1,10), col="#df4f00", lwd=3)
axis(1, at=bar, labels=c(1,"",3,"",5,"",7,"",9,""))
axis(2, at=seq(0,1, length=5), labels=c("0.0", "",0.5,"","1.0"))
bar<-barplot(t(mres2), xaxt="n",yaxt="n", ylim=c(0,1.2), 
     main="Intermediate",xlab="Resource rank", ylab="", col="gray")
lines(x=bar, y=mres2, col="#f1f1f1", lwd=3)
axis(1, at=bar, labels=c(1,"",3,"",5,"",7,"",9,""))
axis(2, at=seq(0,1, length=5), labels=c("0.0", "",0.5,"","1.0"))
bar<-barplot(c(1,rep(0,9)),xaxt="n",yaxt="n",   ylim=c(0,1.2),
     main="Extreme specialist", xlab="", ylab="", col="gray")
lines(x=bar, y=c(1,rep(0,9)), col="#00918d", lwd=3)
axis(1, at=bar, labels=c(1,"",3,"",5,"",7,"",9,""))
axis(2, at=seq(0,1, length=5), labels=c("0.0", "",0.5,"","1.0"))
dev.off()


##################################################################
png(filename="Figures/Exported/Figure S2l1.png", width=4000, height=4000, res=600)

#colfunc <- colorRampPalette(c("#df4f00", "lightgray","#00918d"))
colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
licol<-colfunc(5)
par(las=1)
plot(mres[,6] , mres[,1], type="l", 
     xaxt="n", yaxt="n", 
     ylab="Preferences",
     xlab="Resource rank", cex.lab=1.2,
     col="#df4f00", ylim=c(0,1), xlim=c(0,1))
for (i in 1: (ncol(mres)-1)){
  lines(mres[,6] , mres[,i], col=licol[i], lwd=5)
}

axis(1, at=seq(0,1,length=20), labels=c(1,"",3,"",5,"",7,"",9,"",11,"",13,"",15,"",17,"",19,""))
axis(2, at=seq(0,1, length=5), labels=c("0.0", "",0.5,"","1.0"))
text(x=0.55, y=0.49, labels=1)
text(x=0.75, y=0.7, labels=0.3)
text(x=0.999, y=1, labels=1e-05)
text(x=0.37, y=0.3, labels=3)
text(x=0.12, y=0.025, labels=1e+05)

dev.off()
#########################################
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
lin2<-genfun(t(pr)*1000000, rep(1,nrow(pr)))
lin2

lin3<-wcfun(t(pr), rep(1,nrow(pr)))

toplo<-data.frame(p,lin)
toplo
mycol<-colfunc(nrow(toplo))

png(filename="Figures/Exported/Figure S2l.png", width=5000, height=3600, res=600)
par(mar=c(3,3,2,1))
layout(matrix(c(11,1,2,3,
                11,4,5,6,
                11,7,8,9,
                12,10,10,10), byrow=T, ncol=4), heights=c(30,30,30,5), widths=c(5,30,30,30))
#aPDI
par(las=1)
plot(p, lin, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(alpha*italic(PDI)),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=toplo$lin[i], x1=toplo$p[i+1], y1=toplo$lin[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#Bs
par(las=1)
plot(p, lin2$Bs, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("Bs'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$Bs[i], x1=toplo$p[i+1], y1=lin2$Bs[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#Bprime
par(las=1)
plot(p, lin2$`B'`, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("B''")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$`B'`[i], x1=toplo$p[i+1], y1=lin2$`B'`[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#W
par(las=1)
plot(p, lin2$W, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("W'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$W[i], x1=toplo$p[i+1], y1=lin2$W[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#PS
par(las=1)
plot(p, lin2$PS, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("PS'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$PS[i], x1=toplo$p[i+1], y1=lin2$PS[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)
text(x=12e-05, y= 0.92, label="Generalist" ,col="#df4f00", cex=1)
text(x=8000, y=0.08, "Specialist", col="#00918d", cex=1)


#FT
par(las=1)
plot(p, lin2$FT, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("FT'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$FT[i], x1=toplo$p[i+1], y1=lin2$FT[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#d
par(las=1)
plot(p, lin2$`1-d'`, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("1-d'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$`1-d'`[i], x1=toplo$p[i+1], y1=lin2$`1-d'`[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#gen
par(las=1)
plot(p, lin2$gen, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("gen'")),
     xlab="", xaxt="n", cex.lab=1.2)
axis(1, at=c(1e-04,1e-02,1,1e02,1e04), labels=c(expression(10^-4),expression(10^-2),1,expression(10^2),expression(10^4)))


#lines(p, lin2, col="red")
for (i in 1:(nrow(toplo)-1)){
  segments(x0=toplo$p[i], y0=lin2$gen[i], x1=toplo$p[i+1], y1=lin2$gen[i+1], col=mycol[i], lwd=5)
}
abline(v=1, lty=2)
abline(h=0.5, lty=2)


#Wc
par(las=1)
plot(p, lin3, log="x", type="l",
     ylab="", yaxp=c(0,1,2), main=expression(italic("Wc'")),
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
text(x=0.5, y=0.5, labels="Degree of generalization", cex=1.3, srt=90)

dev.off()