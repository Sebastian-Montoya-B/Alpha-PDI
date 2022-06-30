################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces the Spearman correlations calculated between the
### specialization parameter and the generality indices.


######################### 1. SETTINGS ##########################################


## 1.1 Clean the environment.
rm(list= ls())

## 1.2 Source the functions.
source("Code/alpha_PDI.R")
source("Code/genfun.R")
lisEv<-readRDS("Data/vectors1.RDS")


######################### 2. CALCULATIONS ######################################


##  2.1. Generality 

tt<-lapply(lisEv, function(x){lapply(x, function(x){alpha_PDI(t(x$preference), rep(1,NROW(x$preference)), corrected=F)})})
tt<-unlist(tt)

tt2<-lapply(lisEv, function(x){lapply(x, function(x){genfun(t(x$preference*100), rep(1,NROW(x$preference*100)))})})
tt2<-dplyr::bind_rows(tt2)


## 2.1. Spearman correlations between generality and the specialization parameter

cor.test(spen, tt, method="spearman", exact=F) # alpha PDI
cor.test(spen, tt2$Bs, method="spearman", exact=F) # Bs
cor.test(spen, tt2$`B'`, method="spearman", exact=F) # B'
cor.test(spen, tt2$W, method="spearman", exact=F) # W
cor.test(spen, tt2$PS, method="spearman", exact=F) # PS
cor.test(spen, tt2$FT, method="spearman", exact=F) # FT
cor.test(spen, tt2$`1-d'`, method="spearman", exact=F) #1-d'
cor.test(spen, tt2$gen, method="spearman", exact=F) #gen