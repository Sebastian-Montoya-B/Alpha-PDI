################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script sets the alpha_PDI function for further calculations.


alpha_PDI<-function(data, abun){
  
  q<-abun
  
  if (class(data)[1]!="matrix"){
    data<-matrix(data, ncol=length(data))
  }
  
  if (length(q)!=NCOL(data)){
    stop("The number of resources and their abundance vector do not match")
  }
  

  q<-q/sum(q)
  

  aPDIfun<- function(p,q){
  #This function calculates aPDI
    p<-(p/sum(p))/q #Manly alpha
    p<-p/sum(p)
    sort.p <- sort(p, decreasing = TRUE)
    P.max <- sort.p[1]
    aPDI <- 1-sum(P.max - sort.p[-1])/(length(q) - 1)
    return(aPDI)
  }
  

  aPDImaxFind <- function(x, q) {
    #For the given number of observations find the maximum possible value of aPDI
    #This function is similar to the one used in the bipartite package
    # (Dormann et al. 2008) to find dmin and calculate d'.
    expec <- floor(q * (sum(x)))
    restuse <- sum(x) - sum(expec)
    x.new <- expec
    i.vec <- 1:length(x)
    
    if (restuse!=0){
      for (j in 1:restuse) {
        
        aPDI.check <- numeric(0)
        xsum <- sum(x.new)
        
        p.1 <- x.new/(xsum + 1)/q
        
        for (i in i.vec) {
          pi.1 <- (x.new[i] + 1)/(xsum + 1)/q[i]
          pen<-c(pi.1, p.1[-i])
          pen<-pen/sum(pen)
          sort.p.1 <- sort(pen, decreasing = TRUE)
          P1.max <- sort.p.1[1]
          aPDI.check[i] <- 1-sum(P1.max - sort.p.1[-1])/(length(q) - 1)
        }
        
        i.best <- which.max(aPDI.check)[1]
        x.new[i.best] <- x.new[i.best] + 1
      }
      
    }
    
    return(aPDIfun(x.new, q))
    
  }
  
  aPDI.raw<-apply(data, 1, aPDIfun, q)
  aPDI.max<-apply(data, 1, aPDImaxFind, q)
  aPDI.max<-ifelse(aPDI.raw > aPDI.max, aPDI.raw, aPDI.max)
  aPDI.unb<-aPDI.raw/aPDI.max
  return(list(corrected_aPDI=aPDI.unb,raw_aPDI=aPDI.raw, max_aPDI=aPDI.max ))
}

