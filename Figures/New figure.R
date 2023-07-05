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

if(!require(emdbook)){
  install.packages("emdbook")
  }
library(emdbook)

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

if(!require(grid)){
  install.packages("grid")
}
library(grid)

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/wcfun.R")
source("Code/QNM.R")

## Generate a list of matrices and vectors generated using the quantitative 
## niche model of Fründ et al. (2016)
#if (T) {
#  Nbee <- 5
#  nsim<-1
#  MaUn<-NULL
#  MaEv<-NULL
  #spen<-c(seq(0.1, 60, by=0.5))
  #spen<-lseq(0.01, 50, length =200)#Specialization parameter
#  spen<-seq(0.01, 50, length=200)
  #spen<-c(seq(0.1, 60, length=20)) #Specialization parameter
#  length(spen)
#  lisUn<-NULL
#  lisEv<-NULL
#  spelisUn<-NULL
#  spelisEv<-NULL
#  counter<-1
#  lisnam<-NULL
#  for (Nplant in c(5, 10, 50)){ # Number of potential resources
    
#    for (spe in spen){
      
#      for (i in 1:nsim){
#        MaUn[[i]]<-gen_uneven(Nbee,Nplant, spe, samp=F, make="random" )
        
#      }
#      lisnam[[counter]]<-Nplant
#      lisUn[[counter]]<-MaUn
#      counter<-counter+1
#    }
    
#  }
#  mat1<-lisUn
#}



## For each consumer in mat1 there are four vectors/matrices: 
##   (1) the consumer abundance distribution ($con_abun)
##   (2) the resource abundance distribution ($res_abun)
##   (3) the true preferences ($preference) 
##   (4) the current resource use pattern ($current).


######################### 2. CALCULATIONS ######################################
mat1<-readRDS("Data/sim_data.rds")

if (T) {
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
  
}


######################### 3. PLOTTING ##########################################


#spar<-rep(rep(spen, each=Nbee), 3)

msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}

#colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
#expec_col<-colfunc(10)

#colorin<-function(expected){

#  expec_vec_col<- ifelse(expected>0.9, expec_col[1],
 #                            ifelse(expected>0.8, expec_col[2],
  #                                  ifelse(expected>0.7, expec_col[3],
   #                                        ifelse(expected>0.6,  expec_col[4],
    #                                              ifelse(expected>0.5,  expec_col[5],
     #                                                    ifelse(expected>0.4,  expec_col[6],
      #                                                          ifelse(expected>0.3,  expec_col[7],
       #                                                                ifelse(expected>0.2,  expec_col[8],
        #                                                                      ifelse(expected>0.1,  expec_col[9], expec_col[10]
         #                                                                            )))))))))
  #return(expec_vec_col)
#}

#colorin(lisObv$aPDI)

png(filename="Figures/Exported/Figure2.2.png", width=5000, height=3600, res=600)
if (TRUE){
  spen<-length(mat1)/3
  Nbee<-ncol(mat1[[1]][[1]]$preference)
  par(mar= c(2,3,1,0),las=1, xpd=F)
  layout(matrix(c(14,14,14,14,12,
                  11,1,2,3,12,
                  11,4,5,6,12,
                  11,7,8,9,12,
                  13,10,10,10,12
  ), ncol=5, byrow=T),
  widths=c(5,30,30,30,15),
  heights=c(5,30,30,30,5))
  #layout.show(14)
  
  l5<-length(lisObv$aPDI)/3
  l10<-length(lisObv$aPDI)/3*2
  l50<-length(lisObv$aPDI)
  
  
  x1_xit<-jitter(rep(1, l5), factor=6)
  x2_xit<-jitter(rep(2, l5), factor=6)
  x3_xit<-jitter(rep(3, l5), factor=6)
  
  x_xit<-c(x1_xit,x2_xit,x3_xit)
  colfunc <- colorRampPalette(c("#df4f00","#f1f1f1","#00918d"))
  #colfunc<- colour_ramp(c("#df4f00","#f1f1f1","#00918d"))
  spar_col<-as.raster(matrix(colfunc(spen), nrow=1))
  spar_col<-rep(rep(spar_col, each=Nbee), 3)
  #legend_image <- as.raster(matrix(colfunc(20), nrow=1))
  
  ###aPDI
  
  error1<-data.frame(R5=lisObv$aPDI[1:l5]-lisExv$aPDI[1:l5], 
                     R10=lisObv$aPDI[(l5+1):l10]-lisExv$aPDI[(l5+1):l10], 
                     R50=lisObv$aPDI[(l10+1):l50]-lisExv$aPDI[(l10+1):l50])
  
  boxplot(error1, ylim=c(-1,1), pch=8, xaxt="n", main=expression(alpha*italic(PDI)), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$aPDI)
  points(x=x_xit, y=c(error1$R5,error1$R10,error1$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  msqe<-round(apply(as.matrix(error1), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  par(new=T)
  boxplot(error1, ylim=c(-1,1), pch=8, xaxt="n", main=expression(alpha*italic(PDI)), border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  ###Bs
  
  errorBS<-data.frame("R5"=lisObv$Bs[1:l5]-lisExv$Bs[1:l5],
                      "R10"=lisObv$Bs[(l5+1):l10]-lisExv$Bs[(l5+1):l10],
                      "R50"=lisObv$Bs[(l10+1):l50]-lisExv$Bs[(l10+1):l50])
  
  boxplot(errorBS, ylim=c(-1,1), pch=8, xaxt="n", yaxt="n", main=expression(italic(B[s])), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$Bs)
  points(x=x_xit, y=c(errorBS$R5,errorBS$R10,errorBS$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errorBS), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  par(new=T)
  boxplot(errorBS, ylim=c(-1,1), pch=8, xaxt="n", yaxt="n", main=expression(italic(B[s])) , border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  
  ###B
  
  errorBprime<-data.frame(R5=lisObv$`B'`[1:l5]-lisExv$`B'`[1:l5], 
                          R10=lisObv$`B'`[(l5+1):l10]-lisExv$`B'`[(l5+1):l10],
                          R50=lisObv$`B'`[(l10+1):l50]-lisExv$`B'`[(l10+1):l50])
  
  boxplot(errorBprime, ylim=c(-1,1), pch=8, xaxt="n", yaxt="n", main=expression(italic("B'")), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$`B'`)
  points(x=x_xit, y=c(errorBprime$R5,errorBprime$R10,errorBprime$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  
  msqe<-round(apply(as.matrix(errorBprime), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  par(new=T)
  boxplot(errorBprime, ylim=c(-1,1), pch=8, xaxt="n", yaxt="n", main=expression(italic("B'")) , border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  ###W
  
  errorW<-data.frame(R5=lisObv$W[1:l5]-lisExv$W[1:l5], 
                     R10=lisObv$W[(l5+1):l10]-lisExv$W[(l5+1):l10], 
                     R50=lisObv$W[(l10+1):l50]-lisExv$W[(l10+1):l50])
  
  boxplot(errorW, ylim=c(-1,1), pch=8, xaxt="n", main=expression(italic(W)), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$W)
  points(x=x_xit, y=c(errorW$R5,errorW$R10,errorW$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errorW), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  par(new=T)
  boxplot(errorW, ylim=c(-1,1), pch=8, xaxt="n", main=expression(italic(W)), border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  ###PS
  
  errorPS<-data.frame(R5=lisObv$PS[1:l5]-lisExv$PS[1:l5], 
                      R10=lisObv$PS[(l5+1):l10]-lisExv$PS[(l5+1):l10], 
                      R50=lisObv$PS[(l10+1):l50]-lisExv$PS[(l10+1):l50])
  
  boxplot(errorPS, ylim=c(-1,1),  pch=8, xaxt="n", yaxt="n", main=expression(italic(PS)), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$PS)
  points(x=x_xit, y=c(errorPS$R5,errorPS$R10,errorPS$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errorPS), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  par(new=T)
  boxplot(errorPS, ylim=c(-1,1),  pch=8, xaxt="n", yaxt="n", main=expression(italic(PS)), border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  ###FT
  
  errorFT<-data.frame(R5=lisObv$FT[1:l5]-lisExv$FT[1:l5], 
                      R10=lisObv$FT[(l5+1):l10]-lisExv$FT[(l5+1):l10], 
                      R50=lisObv$FT[(l10+1):l50]-lisExv$FT[(l10+1):l50])
  
  boxplot(errorFT, ylim=c(-1,1),pch=8, xaxt="n", yaxt="n", main=expression(italic(FT)), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$FT)
  points(x=x_xit, y=c(errorFT$R5,errorFT$R10,errorFT$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errorFT), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=F)
  axis(2, at=seq(-1,1, length=5), labels=F)
  par(new=T)
  boxplot(errorFT, ylim=c(-1,1),pch=8, xaxt="n", yaxt="n", main=expression(italic(FT)),border=alpha("black",0.4), col=alpha("gray",0))
  par(new=F)
  
  ###1-d
  
  errord<-data.frame(R5=lisObv$`1-d'`[1:l5]-lisExv$`1-d'`[1:l5], 
                     R10=lisObv$`1-d'`[(l5+1):l10]-lisExv$`1-d'`[(l5+1):l10], 
                     R50=lisObv$`1-d'`[(l10+1):l50]-lisExv$`1-d'`[(l10+1):l50])
  
  boxplot(errord, ylim=c(-1,1), pch=8, xaxt="n", main=expression(1 - italic("d'")), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$`1-d`)
  points(x=x_xit, y=c(errord$R5,errord$R10,errord$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errord), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(1, at=c(1,2,3), labels=c(5,10,50))
  par(new=T)
  boxplot(errord, ylim=c(-1,1), pch=8, xaxt="n", main=expression(1 - italic("d'")), border=alpha("black",0.4),  col=alpha("gray",0))
  par(new=F)
  
  ###gen
  
  errorgen<-data.frame(R5=lisObv$gen[1:l5]-lisExv$gen[1:l5], 
                       R10=lisObv$gen[(l5+1):l10]-lisExv$gen[(l5+1):l10], 
                       R50=lisObv$gen[(l10+1):l50]-lisExv$gen[(l10+1):l50])
  
  boxplot(errorgen, ylim=c(-1,1), pch=8, yaxt="n", xaxt="n",  main=expression(italic(gen)), col=alpha("gray",0))
  #spar_col<-colorin(lisExv$gen)
  points(x=x_xit, y=c(errorgen$R5,errorgen$R10,errorgen$R50),
         bg= alpha(spar_col,0.2), pch=21, col=alpha(spar_col,0.2), cex=0.8)
  
  msqe<-round(apply(as.matrix(errorgen), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(2, at=seq(-1,1, length=5), labels=F)
  axis(1, at=c(1,2,3), labels=c(5,10,50))
  par(new=T)
  boxplot(errorgen, ylim=c(-1,1), pch=8, yaxt="n", xaxt="n",  main=expression(italic(gen)),border=alpha("black",0.4),  col=alpha("gray",0))
  par(new=F)
  
  ###Wc
  
  errorWc<-data.frame(R5=lisObvwc[1:l5]-lisExvwc[1:l5], 
                      R10=lisObvwc[(l5+1):l10]-lisExvwc[(l5+1):l10], 
                      R50=lisObvwc[(l10+1):l50]-lisExvwc[(l10+1):l50])
  
  
  boxplot(errorWc, ylim=c(-1,1), pch=8, yaxt="n", xaxt="n", main=expression(italic(Wc)), col=alpha("gray",0))
  #spar_col<-colorin(lisExvwc)
  points(x=x_xit, y=c(errorWc$R5,errorWc$R10,errorWc$R50),
         bg= alpha(spar_col,0.3), pch=21, col=alpha(spar_col,0.3), cex=0.8)
  
  errorWc<-na.omit(errorWc) ## Wc may me unable to estimate the value for some consumers, and it will produce NAs
  msqe<-round(apply(as.matrix(errorWc), 2, msqerrorfun),3)
  msqe<-sprintf(msqe, fmt='%#.3f')
  text(x=1, y=-1, labels=bquote("MSE = "*.(msqe[1])), cex=0.7)
  text(x=2, y=-1, labels=bquote("MSE = "*.(msqe[2])), cex=0.7)
  text(x=3, y=-1, labels=bquote("MSE = "*.(msqe[3])), cex=0.7)
  abline(h=0, lty=2)
  axis(2, at=seq(-1,1, length=5), labels=F)
  axis(1, at=c(1,2,3), labels=c(5,10,50))
  par(new=T)
  boxplot(errorWc, ylim=c(-1,1), pch=8, yaxt="n", xaxt="n", main=expression(italic(Wc)),border=alpha("black",0.4),  col=alpha("gray",0))
  par(new=F)
  
  ### X Label
  
  par(mar= c(0,0,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.535, y=0.7, labels="Number of potential resources", cex=1)
  
  ### Y Label
  
  par(mar= c(0,0.5,0,0),las=1)
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.5, labels="Difference between the observed and expected degree of generalization", cex=1, srt=90)
  
  ###### Legend
  
  par(mar= c(4,2,4,0),las=1)
  plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialization\nparameter",
       cex.main=1, ylim=c(0,1))
  
  #rasterImage(legend_image, xleft=0, ybottom=0, xright=1, ytop=1)
  
  
  legend_image <- as.raster(rev(colfunc(20)))
  grid.raster(legend_image, width=0.04, height = 0.8, x = unit(0.93, "npc"))
  par(xpd=T)
  abline(v=-0.2, lty=3)
  text(y=c(-0.02,1.03), x =0.7 , labels = c(0.01, "50.0"), cex=1)
  par(xpd=F)
  #par(mar= c(0,0,2,0),las=1)
  #plot(y=c(-1,2),x=c(-1,1),type = 'n', axes = F,xlab = '', ylab = '', main= "Specialization\nparameter",
  #     cex.main=1)
  #text(y=c(0.0,1,1.9), x =0.5 , labels = c(0.01,"","50.0"), cex=1)
  #rasterImage(legend_image, xleft=0, ybottom=0, xright=1, ytop=1, angle=90)
}


dev.off()

