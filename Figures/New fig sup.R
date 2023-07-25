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

if(!require(purrr)){
  install.packages("purrr")
}
library(purrr)

if(!require(grid)){
  install.packages("grid")
}
library(grid)

if(!require(dplyr)){
  install.packages("dplyr")
}
library(dplyr)

## Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")



######################### 2. CALCULATIONS ######################################
msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}

lisUn<-mat1
if (T){
  
  lisUexv2<-lapply(lisUn, function(x){ genfun(t(x$large), x$res_abun)})
  lisUexv2<-dplyr::bind_rows(lisUexv2)
  
  lisUexv2_10<-lapply(lisUn, function(x){ genfun(t(x$small1), x$res_abun)})
  lisUexv2_10<-unlist(lisUexv2_10)
  
  lisUexv2_50<-lapply(lisUn, function(x){ genfun(t(x$small2), x$res_abun)})
  lisUexv2_50<-unlist(lisUexv2_50)
  
  lisUexv2_100<-lapply(lisUn, function(x){ genfun(t(x$small3), x$res_abun)})
  lisUexv2_100<-unlist(lisUexv2_100)
  
  lisUexv2_500<-lapply(lisUn, function(x){ genfun(t(x$small4), x$res_abun)})
  lisUexv2_500<-unlist(lisUexv2_500)
  
}
