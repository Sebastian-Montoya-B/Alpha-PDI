################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script sets the alpha_PDI function for further calculations.


alpha_PDI<-function(data, abun, corrected=T, m=1){
  
  q<-abun
  
  if (class(data)[1]!="matrix"){
    data<-matrix(data, ncol=length(data))
  }
  
  if (length(q)!=NCOL(data)){
    stop("The number of resources and their abundance vector do not match")
  }
  
  if (corrected==T & any(data%%1!=0)){
    stop("You must provide non-negative integer values to calculate the corrected version of alpha PDI")
  }
  
  
  q<-q/sum(q)
  
  aPDIfun<- function(p,q,m){
    #This function calculates aPDI
    #aPDI<-1-(sum(1-(p/(q*max(p/q)))))/(length(q)-1)
    aPDI<-1-((sum((1-(p/(q*max(p/q))))^m))/(length(q)-1))
    return(aPDI)
  }
  
  if (corrected==T){
    aPDImaxFind <- function(x, q, m) {
      #For the given number of observations find the maximum possible value of aPDI
      #This function is similar to the one used in the bipartite package (Dormann et al. 2008) to find dmin and calculate d'
      #expec <- floor(q * (sum(x)-length(x)))
      #expec<-expec + 1
      
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
            pen<-pen/max(pen)#
            sort.p.1 <- sort(pen, decreasing = TRUE)
            P1.max <- sort.p.1[1]
            aPDI.check[i] <- 1-sum((P1.max - sort.p.1[-1])^m)/(length(q) - 1)
          }
          
          i.best <- which.max(aPDI.check)[1]
          x.new[i.best] <- x.new[i.best] + 1
        }
        
      }
      
      return(aPDIfun(x.new, q, m))
      
    }
    
    aPDImaxFind2 <- function(x, q, m) {
      #For the given number of observations find the maximum possible value of aPDI
      
      
      expec <- floor(q * (sum(x)-length(x)))
      expec<-expec + 1
      
      #expec <- floor(q * (sum(x)))
      restuse <- sum(x) - sum(expec)
      
      if (any(expec<0)){
        restuse<-restuse+sum(expec[which(expec<0)])
        expec[which(expec<0)]<-0
      }
      
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
            pen<-pen/max(pen)#
            sort.p.1 <- sort(pen, decreasing = TRUE)
            P1.max <- sort.p.1[1]
            aPDI.check[i] <- 1-sum((P1.max - sort.p.1[-1])^m)/(length(q) - 1)
          }
          
          i.best <- which.max(aPDI.check)[1]
          x.new[i.best] <- x.new[i.best] + 1
        }
        
      }
      
      return(aPDIfun(x.new, q, m))
      
      
      
      
    }
    
    aPDI.raw<-apply(data, 1, aPDIfun, q, m)
    aPDI.max1<-apply(data, 1, aPDImaxFind, q, m)
    aPDI.max2<-apply(data, 1, aPDImaxFind2, q, m)
    aPDI.max<-ifelse(aPDI.max1 > aPDI.max2, aPDI.max1, aPDI.max2)
    aPDI.max<-ifelse(aPDI.raw > aPDI.max, aPDI.raw, aPDI.max)
    aPDI.unb<-ifelse(aPDI.max==0,0, aPDI.raw/aPDI.max)
    
    return(list(corrected_aPDI=aPDI.unb,raw_aPDI=aPDI.raw, max_aPDI=aPDI.max ))
    
  } else{
    aPDI.raw<-apply(data, 1, aPDIfun, q, m)
    return(raw_aPDI=aPDI.raw)
  }
  
}
