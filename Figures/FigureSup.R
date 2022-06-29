#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
#             Boris R. Krasnov, Marco A. R. Mello

#### See README for further info 
#    https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme

#### This script reproduces Figure S1 to S36.


#_______________________________________________________________________________


############### SUMMARY ###############

#           1. SETTING
#           2. CREATING FUNCTIONS
#           4. CALCULATING
#           5. PLOTTING AND EXPORTING

#######################################

# 1. SETTING

source("Code/alpha_PDI.R")
lisEv<-readRDS("Data/vectors432.RDS")

## For each consumer in lisEv there are five vectors: 
##   (1) the resource abundance distribution ($res_abun),
##   (2) the true preferences ($preference), 
##   (3) the current pattern of resource use ($current),
##   (4) the observed pattern of resource use in a case with 10^6 observations ($large), and
##   (5) the observed pattern of resource use in a case with 100 observations ($small)
## Here we will focus on the $large vector

# 2. CREATING FUNCTIONS

# 2.1. Sampling simulation function

samp.sim<-function(prob, n.obs, n.rep){
  
  prob<-as.vector(prob)
  prob<-prob/sum(prob)
  sam.seq<-array(rep(0, length(prob)*(n.obs+1)), dim=c(n.obs+1,length(prob),n.rep))
  for (k in 1:n.rep){
    
    for (i in 1:n.obs){
      new.int<-sample(1:length(prob), 1, replace=F, prob=as.vector(prob))
      new.vec<-sam.seq[i,,k]
      new.vec[new.int]<-new.vec[new.int]+1
      sam.seq[i+1,,k]<-new.vec
    }
    
    
  }
  
  return(sam.seq)
}



# 2.2. Function to find the minimum number of observations needed for an accurate estimation
#      This function is based on Fründ et al. (2016)

findmin<-function(corr.apdi, expected, accuracy=0.05){
  
  difseq<-abs(corr.apdi-expected) 
  
  if (all(difseq < accuracy,na.rm=T)) {
    nsat <- 1
    if(all(is.na(difseq))) nsat <- NA
  } else {
   
    nsat.mx<-apply(difseq, 2, function(x){ifelse(any(which(x >= 0.05))==0,1,max(which(x >= 0.05))+1)})
    
    nsat<-mean(nsat.mx)
  } 
  return(nsat)
}

# 4. CALCULATING

# 4.1. Simulate the sampling process for each consumer
n.rep<-10
n.obs<-100
fga<-lapply(lisEv, function(x){lapply(x, function(x){samp.sim(x$large,n.obs,n.rep)})})

# 4.2. Calculate alpha PDI for each step of the sampling process
#      This calculation may take several hours based on the 
#      number of repetitions and observations used. For the current repetitions and observations
#      it took up to 6 hours in a PC with the following specifications:
#      Windows 10 with AMD A12 2.70GHz, RAM 12 GB 


respp<-lapply(fga, function(x){lapply(x, function(x){apply(x, 3, function(x){alpha_PDI(x[-1,],rep(10, NCOL(x)))$corrected_aPDI})})})



# 4.3. Calculate expected values based on a vector of 10^6 observations ($large)
lisEexv<-lapply(lisEv, function(x){lapply(x, function(x){alpha_PDI(t(x$large), rep(1,NROW(x$large)))$raw_aPDI})})
lisEexv<-unlist(lisEexv)


# 5. PLOTTING AND EXPORTING
colw<-c("#00ceff", "#078ab5","#004c6d")

for (i in 1:length(respp)){
  
  svg(filename=paste0("FigSup", i, ".png"), width=8, height=9)
  layout(matrix(c(13,1,2,3,13,4,5,6,13,7,8,9,13,10,11,12,13,14,14,14), byrow=T, ncol=4), 
         widths=c(10,30,30,30), heights=c(22,22,22,22,12))
  layout.show(14)
  par(mar=c(2,3,1,1), las=1)
 
  for (j in 1:length(respp[[i]])){
    plot(x=seq(1,100), y=respp[[i]][[j]][,1], type="l", ylim=c(0,1), xlab="",ylab=i, col="#0c0c0c")
    
    for (k in 2:n.rep){
      lines(x=seq(1,100), y=respp[[i]][[j]][,k], col="#0c0c0c")
      
    }
    abline(v=findmin(respp[[i]][[j]], lisEexv[(i-1)*length(respp[[i]])+j]), col="#ee9b00", lwd=2)
    
    if (i<13){
      points(x=101, y=lisEexv[(i-1)*length(respp[[i]])+j], col=colw[3], pch=21, cex=2, bg=colw[3])
    } else if (12<i&i<25){
      points(x=101, y=lisEexv[(i-1)*length(respp[[i]])+j], col=colw[2], pch=23, cex=2, bg=colw[2])
    } else if (24<i){
      points(x=101, y=lisEexv[(i-1)*length(respp[[i]])+j], col=colw[1], pch=22, cex=2, bg=colw[1])
    }
    
  }
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.55, labels=expression("Estimated generality   "*alpha*italic(PDI)), srt=90, cex=1.6)
  
  plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
  text(x=0.5, y=0.5, labels="Sampling intensity (number of observations)", cex=1.6)
  dev.off()
}







      