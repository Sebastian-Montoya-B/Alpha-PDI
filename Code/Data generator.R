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
if(!require(zCompositions)){
  install.packages("zCompositions")
  library(zCompositions)
}

if(!require(bipartite)){
  install.packages("bipartite")
  library(bipartite)
}


if(!require(emdbook)){
  install.packages("emdbook")
}
library(emdbook)

if(!require(dplyr)){
  install.packages("dplyr")
  library(dplyr)
}

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/wcfun.R")
source("Code/QNM.R")

## Generate a list of vectors using the quantitative niche model of 
## Fründ et al. (2016)

if (T) {
  Nbee <- 1
  #nsim<-1
  MaUn<-NULL
  MaEv<-NULL
  #spen<-c(seq(0.1, 60, by=0.5))
  spen<-lseq(0.01, 3000, length =2000)#Specialization parameter
  #spen<-seq(0.01, 400, length=500)
  #spen<-c(seq(0.1, 60, length=20)) #Specialization parameter
  length(spen)
  lisUn<-NULL
  lisEv<-NULL
  spelisUn<-NULL
  spelisEv<-NULL
  counter<-1
  lisnam<-NULL
  for (Nplant in c(5, 15, 55)){ # Number of potential resources
    
    for (spe in spen){
      
      #for (i in 1:nsim){
        #MaUn[[i]]<-gen_uneven(Nbee,Nplant, spe, samp=F, make="random" )
        MaUn<-gen_uneven2(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),maxsamp=1000, make="spread" )
      #}
      lisnam[[counter]]<-Nplant
      lisUn[[counter]]<-MaUn
      counter<-counter+1
    }
    
  }
  mat1<-lisUn
}

length(mat1)

saveRDS(mat1, "Data/sim_data.RDS")

#####
if (T) {
  Nbee <- 2
  #nsim<-1
  MaUn<-NULL
  MaEv<-NULL
  #spen<-c(seq(0.1, 60, by=0.5))
  spen<-lseq(0.01, 3000, length =2000)#Specialization parameter
  #spen<-seq(0.01, 400, length=500)
  #spen<-c(seq(0.1, 60, length=20)) #Specialization parameter
  length(spen)
  lisUn<-NULL
  lisEv<-NULL
  spelisUn<-NULL
  spelisEv<-NULL
  counter<-1
  lisnam<-NULL
  for (Nplant in c(5, 15, 55)){ # Number of potential resources
    
    for (spe in spen){
      
      #for (i in 1:nsim){
      #MaUn[[i]]<-gen_uneven(Nbee,Nplant, spe, samp=F, make="random" )
      MaUn<-gen_uneven2(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),maxsamp=1000, make="spread" )
      #}
      lisnam[[counter]]<-Nplant
      lisUn[[counter]]<-MaUn
      counter<-counter+1
    }
    
  }
  mat2<-lisUn
}

length(mat2)

rere<- list()
for (i in 1: length(mat2)){
  rere[[i]]<-mat2[[i]]$res_abun
  
}

rere
saveRDS(mat2, "Data/sim_data2.RDS")

#####

if (T) {
  Nbee <- 1
  #nsim<-1
  MaUn<-NULL
  MaEv<-NULL
  #spen<-c(seq(0.1, 60, by=0.5))
  #spen<-lseq(0.01, 400, length =500)#Specialization parameter
  spen2<-seq(0.1, 50, length=2000)
  #spen<-c(seq(0.1, 60, length=20)) #Specialization parameter
  length(spen)
  lisUn<-NULL
  lisEv<-NULL
  spelisUn<-NULL
  spelisEv<-NULL
  counter<-1
  lisnam<-NULL
  nres<-length(mat2)/3
  volt<-1
  ini<-1
  for (Nplant in c(5, 15, 55)){ # Number of potential resources
    
    
    for (spe in spen2){
      
      
        for (i in ini:(nres*volt)){
          #MaUn[[i]]<-gen_uneven(Nbee,Nplant, spe, samp=F, make="random" )
          MaUn[[i]]<-gen_uneven3(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),
                                 maxsamp=1000, make="random", res_abun = rere[[i]] )
        }
    
      
      lisnam[[counter]]<-Nplant
      lisUn[[counter]]<-MaUn
      counter<-counter+1
      ini<-(nres*volt)+1
      volt<-volt+1
    }
    
  }
  mat3<-lisUn
}

head(rere[[i]])
rere[[6000]]
length(rere)
str(rere)


gen_uneven3(Nbee,Nplant, spe, samp=T,minsamp=c(10,50,80,100),
            maxsamp=1000, make="random", res_abun = rere[[6000]] )
is.vector(rere[[1]])
class(rere[[1]])
mat1[[1]]
length(mat1)
web_p <- makeweb(specpar=spe, Nbee, Nplant) 

mat2[[1]]
length(mat2)
length(mat2[[2]][[1]])
