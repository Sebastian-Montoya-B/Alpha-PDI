################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script generates the simulated data for furhter analysis.


######################### 1. SETTINGS ##########################################


## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.

if(!require(emdbook)){
  install.packages("emdbook")
}
library(emdbook)

if(!require(purrr)){
  install.packages("purrr")
}
library(purrr)

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

## Source the functions.
source("Code/QNM.R")

## Generate a list of vectors using the quantitative niche model of 
## Fründ et al. (2016)

if (T) {
  Nbee <- 1
  #nsim<-1
  MaUn<-NULL
  MaEv<-NULL
  #spen<-lseq(0.01, 3000, length =2000)# Specialization parameter
  spen<-lseq(0.1, 60, length =2000)# Specialization parameter
  lisUn<-NULL
  lisEv<-NULL
  spelisUn<-NULL
  spelisEv<-NULL
  counter<-1
  lisnam<-NULL
  for (Nplant in c(5, 15, 55)){ # Number of potential resources
    
    for (spe in spen){
      MaUn<-gen_uneven2(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,100,500),maxsamp=1000000, make="spread" )

      lisnam[[counter]]<-Nplant
      lisUn[[counter]]<-MaUn
      counter<-counter+1
    }
    
  }
  mat1<-lisUn
  saveRDS(mat1, "Data/sim_data.RDS")
  ## For each consumer in mat 1 there are eight vectors: 
  ##   (1) the resource abundance distribution ($res_abun)
  ##   (2) the true preferences ($preference)
  ##   (3) the current pattern of resource use ($current)
  ##   (4) the observed pattern of resource use in a case with 10^6 observations ($large)
  ##   (5) the observed pattern of resource use in a case with 10 observations ($small1)
  ##   (6) the observed pattern of resource use in a case with 50 observations ($small2)
  ##   (7) the observed pattern of resource use in a case with 100 observations ($small3)
  ##   (8) the observed pattern of resource use in a case with 500 observations ($small4)
  
  
  rere<- list()
  for (i in 1: length(mat1)){
    rere[[i]]<-mat1[[i]]$res_abun
    
  }

  ## Generating matrices for Wc calculation
  Nbee <- 20
  MaUn2<-NULL
  spen2<-sample(lseq(0.1,60, length =2000), length(mat1), replace=T)
  length(spen)
  nres<-length(mat1)/3
  resvec<-c(rep(5,nres),rep(15,nres),rep(55,nres))

  for (i in 1:(nres*3)){
    print(i)
    while(T){
      newvec<-try(gen_uneven3(Nbee,Nplant=resvec[i], spen2[i], samp=T,minsamp=c(10,50,100,500),
                              maxsamp=1000000, make="random", res_abun = rere[[i]]), silent=T)
      if(!is(newvec, "try-error")) break
    }
    MaUn2[[i]]<-newvec
  }
  
  
  
  
  mat2<-list()
  for (i in 1: length(mat1)){ 
    mat2[[i]]<-map2(mat1[[i]], MaUn2[[i]], ~ cbind(.x,.y))
    
  }

  saveRDS(mat2, "Data/sim_data2.RDS")
  
  
  
  
}



