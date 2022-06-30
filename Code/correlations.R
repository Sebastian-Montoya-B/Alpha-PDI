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
lisEv<-readRDS("Data/vectors1.RDS") #See README for details
spen<-c(0.5,1, seq(2,60,by=2)) #Specialization parameter used to generate lisEv


######################### 2. CALCULATIONS ######################################


##  2.1. Generality 

tt<-lapply(lisEv, function(x){
  lapply(x, function(x){
    alpha_PDI(t(x$preference),
              rep(1,NROW(x$preference)),
              corrected=F)})})
tt<-unlist(tt)

tt2<-lapply(lisEv, function(x){
  lapply(x, function(x){
    genfun(t(x$preference*100),
           rep(1,NROW(x$preference*100)))})})
tt2<-dplyr::bind_rows(tt2)


## 2.1. Spearman correlations between generality indices and the
## specialization parameter

# Plot the correlations
par(mfrow = c(2,4), mar = c(5, 4, 1, 1))
plot(spen, tt, xlab = NA, ylab = "αPDI")
plot(spen, tt2$Bs, xlab = NA, ylab = "Bs")
plot(spen, tt2$`B'`, xlab = NA, ylab = "B'")
plot(spen, tt2$W, xlab = NA, ylab = "W")
plot(spen, tt2$PS, xlab = NA, ylab = "PS")
plot(spen, tt2$FT, xlab = NA, ylab = "FT")
plot(spen, tt2$`1-d'`, xlab = NA, ylab = "1-d'")
plot(spen, tt2$gen, xlab = NA, ylab = "gen")
mtext("Specialization parameter", side = 1, line = -2, outer = T)
par(mfrow=c(1,1))

# Calculate the correlations
cor.test(spen, tt, method="spearman", exact=F) # αPDI
cor.test(spen, tt2$Bs, method="spearman", exact=F) # Bs
cor.test(spen, tt2$`B'`, method="spearman", exact=F) # B'
cor.test(spen, tt2$W, method="spearman", exact=F) # W
cor.test(spen, tt2$PS, method="spearman", exact=F) # PS
cor.test(spen, tt2$FT, method="spearman", exact=F) # FT
cor.test(spen, tt2$`1-d'`, method="spearman", exact=F) #1-d'
cor.test(spen, tt2$gen, method="spearman", exact=F) #gen
