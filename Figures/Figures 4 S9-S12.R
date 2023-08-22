################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure 4, and Figures S9 to S12


######################### 1. SETTINGS ##########################################

## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.

if(!require(purrr)){
  install.packages("purrr")
}
library(purrr)

if(!require(scales)){
  install.packages("scales")
}
library(scales)


if(!require(grid)){
  install.packages("grid")
}
library(grid)


## Source the functions.
source("Code/alpha_PDI.R")


## Load the empirical flea-mammal networks (matrices) 
## and its resource abundance distributions (vectors)

consumers<-readRDS("Data/Fleas.RDS")
resources<-readRDS("Data/resource_abundances.RDS")
Rres<-lapply(resources, length)
Rres<-map2(Rres, consumers, ~rep(.x, nrow(.y)))
Rres<-do.call(c, Rres)
######################### 2. CALCULATIONS ######################################


## 2.1. Calculating αPDI

results<-map2(consumers,resources, ~alpha_PDI(.x,.y,corrected=T)$corrected_aPDI)
length(results)


class.gen<-sapply(results, function (x){ sum(x>0.5)/length(x)})# Generalists
class.spe<-sapply(results, function (x){ sum(x<0.5)/length(x)})# Specialists
class.mid<-sapply(results, function (x){ sum(x==0.5)/length(x)})# Midpoint

summ_class<-data.frame(prop=c(class.gen, class.spe), class=c(rep(1, length(class.gen)),rep(2, length(class.spe)) ))



n.inter<-lapply(consumers, rowSums)
length(n.inter)

summ_samp<-data.frame(aPDI=do.call(c, results), NI=do.call(c, n.inter))
nam<-strsplit(rownames(summ_samp), "\\.")
head(nam)
summ_tab<-cbind(summ_samp, do.call(rbind, nam))
names(summ_tab)[3:4]<-c("site","species")
names(summ_tab)

## 2.2. Correlations
c1<-cor.test(summ_tab$aPDI, summ_tab$NI, method="spearman", exact=F)

summ_samp$R<-Rres

c2<-cor.test(summ_tab$aPDI, summ_samp$R, method="spearman", exact=F)



######################### 3. PLOTTING ######################################

## 3.1. Color settings
color_group<-function(gen_fleas){
  cor.cla<-rep(NA, length(gen_fleas))
  for (i in 1:length(gen_fleas)){
    ifelse(gen_fleas[i]<=0.1,cor.cla[i]<-1,
           ifelse (gen_fleas[i]>0.1 & gen_fleas[i]<=0.2, cor.cla[i]<-2,
                   ifelse(gen_fleas[i]>0.2 & gen_fleas[i]<=0.3, cor.cla[i]<-3, 
                          ifelse(gen_fleas[i]>0.3 & gen_fleas[i]<=0.4, cor.cla[i]<-4, 
                                 ifelse(gen_fleas[i]>0.4 & gen_fleas[i]<0.5, cor.cla[i]<-5,
                                        ifelse(gen_fleas[i]==0.5, cor.cla[i]<-6, 
                                               ifelse(gen_fleas[i]>0.5 & gen_fleas[i]<=0.6, cor.cla[i]<-7,
                                                      ifelse(gen_fleas[i]>0.6 & gen_fleas[i]<=0.7,cor.cla[i]<-8,
                                                             ifelse(gen_fleas[i]>0.7 & gen_fleas[i]<=0.8, cor.cla[i]<-9, 
                                                                    ifelse(gen_fleas[i]>0.8 & gen_fleas[i]<=0.9, cor.cla[i]<-10, 
                                                                           cor.cla[i]<-11))))))))))
    
    
  }
  return(cor.cla)
}


color_group(summ_samp[,1])

colfunc <- colorRampPalette(c("#00918d","#f1f1f1", "#df4f00"))
spar_col<-as.raster(matrix(colfunc(11), nrow=1))
#spar_col<-rep(rep(spar_col, each=1), 3)

spar_col[c(color_group(summ_samp[,1]))]

## 3.2. Plotting Figure 4
png(filename="Figures/Exported/Figure4.png", width=5000, height=1500, res=600)
layout(matrix(c(1,2,3), ncol=3))


par(mar=c(4,5,2,2),las=1)
plot(summ_samp$NI, summ_samp$aPDI,log="x",  pch=21, main="A",
     ylab=expression("Degree of generalization ("*alpha*italic("PDI')")), yaxt="n",
     xlab=expression("Sampling intensity (log("*italic(n)*"))"),
     bg=alpha(spar_col[c(color_group(summ_samp[,1]))],0.4), 
     col=alpha(spar_col[c(color_group(summ_samp[,1]))],0.4), cex=1.2 )
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c("0.0","",0.5,"","1.0"))
text(x=1e+04, y=1,labels=bquote(italic("rho")*" = "*.(round(c1$estimate,2))*", P < 0.001"))

plot(summ_samp$R, summ_samp$aPDI,  pch=21, main="B",
     ylab=expression("Degree of generalization ("*alpha*italic("PDI')")), yaxt="n",
     xlab=expression("Number of potential resources ("*italic(R)*")"),
     bg=alpha(spar_col[c(color_group(summ_samp[,1]))],0.4), 
     col=alpha(spar_col[c(color_group(summ_samp[,1]))],0.4), cex=1.2 )
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c("0.0","",0.5,"","1.0"))
text(x=22.5, y=1,labels=bquote(italic("rho")*" = "*.(round(c2$estimate,2))*", P = "*.(round(c2$p.value,3))))

x1<-jitter(rep(1, length(class.gen)), factor=6)
x2<-jitter(rep(2, length(class.spe)), factor=6)
plot( x=NULL, y=NULL, xlim=c(0.5,2.5), ylim=c(0,1),
      main="C", ylab="Proportion of species in each network", xlab="", xaxt="n", yaxt="n")
axis(1, at=c(1,2), labels=c("Generalists", "Specialists"))
for (i in 1:length(class.gen)){
  lines(x=c(x1[i],x2[i]), c(class.gen[[i]],class.spe[[i]]),col=alpha("gray",0.7), lty=2  )
}
points(x=x1, class.gen, bg=alpha(spar_col[11],0.4), col=alpha(spar_col[11],0.4), pch=24, cex=1.2)
points(x=x2, class.spe,bg=alpha(spar_col[1],0.4), col=alpha(spar_col[1],0.4), pch=24, cex=1.2)
axis(2, at=c(0,0.25,0.5,0.75,1), labels=c("0.0","",0.5,"","1.0"))

dev.off()



## 3.3. Plotting Figure S9

png(filename="Figures/Exported/FigureS9.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))

for (i in 1:20){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
 
  
}
par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels=expression("Sampling intensity ("*italic("n")*")"), srt=90, cex=1.5)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI)*"’"), cex=1.5)

dev.off()

## 3.4. Plotting Figure S10

png(filename="Figures/Exported/FigureS10.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))



for (i in 21:40){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  
  
}
par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels=expression("Sampling intensity ("*italic("n")*")"), srt=90, cex=1.5)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI)*"’"), cex=1.5)
dev.off()

## 3.5. Plotting Figure S11

png(filename="Figures/Exported/FigureS11.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))



for (i in 41:60){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  
  
}
par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels=expression("Sampling intensity ("*italic("n")*")"), srt=90, cex=1.5)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI)*"’"), cex=1.5)
dev.off()


## 3.6. Plotting Figure S12

png(filename="Figures/Exported/FigureS12.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(17,1,2,3,4,17,5,6,7,8,17,9,10,11,12,17,13,14,15,16,17,18,18,18,18,17,18,18,18,18),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))



for (i in 61:74){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  
  
}

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")

par(mar=c(0,0,0,0))
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.65, labels=expression("Sampling intensity ("*italic("n")*")"), srt=90, cex=1.5)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.95, labels=expression(alpha*italic(PDI)*"’"), cex=1.5)
dev.off()

