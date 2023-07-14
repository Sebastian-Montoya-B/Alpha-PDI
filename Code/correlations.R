################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces the Spearman correlations calculated between the
### specialization parameter and the generalization indices.


######################### 1. SETTINGS ##########################################


## Clean the environment
rm(list= ls())

## Check the required packages, install them if necessary, and load them.
if(!require(emdbook)){
  install.packages("emdbook")
}
library(emdbook)

## Source the functions
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/QNM.R")

## Load the data
mat1<-readRDS("Data/sim_data.rds")
mat2<-readRDS("Data/sim_data2.rds")

spen<-lseq(0.1, 60, length =2000)
######################### 2. CALCULATIONS ######################################


## Generalization 
### The calculation of d' requires integers in order to estimate d_min.
### Therefore, the preference vectors have been multiplied by 1000000.
### Since indices work with proportions, this multiplication does not change
### their results, but allows us to use d'.

tt<-sapply(mat1, function(x){alpha_PDI(t(x$preference),rep(1,NROW(x$preference)), corrected=F)})
tt2<-lapply(mat1, function(x){genfun(t(x$preference*100000),rep(1,NROW(x$preference)))})
tt2<-dplyr::bind_rows(tt2)



## Spearman correlations between generalization indices and the
## specialization parameter

# Plot the correlations
par(mfrow = c(2,4), mar = c(5, 4, 1, 1), las=1)
plot(rep(spen,3), tt, xlab = NA, ylab = "αPDI")
plot(rep(spen,3), tt2$Bs, xlab = NA, ylab = "Bs")
plot(rep(spen,3), tt2$`B'`, xlab = NA, ylab = "B'")
plot(rep(spen,3), tt2$W, xlab = NA, ylab = "W")
plot(rep(spen,3), tt2$PS, xlab = NA, ylab = "PS")
plot(rep(spen,3), tt2$FT, xlab = NA, ylab = "FT")
plot(rep(spen,3), tt2$`1-d'`, xlab = NA, ylab = "1-d'")
plot(rep(spen,3), tt2$gen, xlab = NA, ylab = "gen")
mtext("Specialization parameter", side = 1, line = -2, outer = T)
par(mfrow=c(1,1))

# Calculate the correlations
## when the number of resources is small (e.g., 5), the value of specialization 
## parameter needed for extreme specialization (using exclusively a resource)
## is not as high as the maximum used for all simulations (60). This means that
## indices get to 0 (their minimum) before the specialization parameter gets to 
## 60. For a correlation, this means that the vector of index values has a big 
## sequence of 0s in their tail (e.g., tt[1:2000]). 
## That sequence of 0s disrupts the correlation analysis, indicating that 
## variables are less correlated than what they really are. Therefore, for those 
## cases, correlations were re-calculated removing the sequence of 0 at the end 
## of the vectors (e.g., tt[1:1610]).

## αPDI
cor.test(spen, tt[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1610], tt[1:1610], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt[4001:6000], method="spearman", exact=F) # 55 resources

## Bs
cor.test(spen, tt2$Bs[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1614], tt2$Bs[1:1614], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$Bs[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$Bs[4001:6000], method="spearman", exact=F) # 55 resources

## B'
cor.test(spen, tt2$`B'`[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1614], tt2$`B'`[1:1614], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$`B'`[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$`B'`[4001:6000], method="spearman", exact=F) # 55 resources

## W
cor.test(spen, tt2$W[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1627], tt2$W[1:1627], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$W[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$W[4001:6000], method="spearman", exact=F) # 55 resources

# PS
cor.test(spen, tt2$PS[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1613], tt2$PS[1:1613], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$PS[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$PS[4001:6000], method="spearman", exact=F) # 55 resources

# FT
cor.test(spen, tt2$FT[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1723], tt2$FT[1:1723], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$FT[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$FT[4001:6000], method="spearman", exact=F) # 55 resources

#1-d'
cor.test(spen, tt2$`1-d'`[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1627], tt2$`1-d'`[1:1627], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$`1-d'`[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$`1-d'`[4001:6000], method="spearman", exact=F) # 55 resources

#gen
cor.test(spen, tt2$gen[1:2000], method="spearman", exact=F) # 5 resources
cor.test(spen[1:1622], tt2$gen[1:1622], method="spearman", exact=F) # 5 resources removing 0s at the end.
cor.test(spen, tt2$gen[2001:4000], method="spearman", exact=F) # 15 resources
cor.test(spen, tt2$gen[4001:6000], method="spearman", exact=F) # 55 resources
