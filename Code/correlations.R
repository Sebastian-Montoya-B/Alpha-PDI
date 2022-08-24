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

## Source the functions
source("Code/alpha_PDI.R")
source("Code/genfun.R")
source("Code/QNM.R")

## Generate vectors using the quantitative niche model of 
## Fründ et al. (2016)

Nbee <- 1
nsim<-1
MaUn<-NULL
MaEv<-NULL
spen<-c(seq(0.1, 60, by=0.1)) # Specialization parameter
length(spen)
lisUn<-NULL
lisEv<-NULL
spelisUn<-NULL
spelisEv<-NULL
counter<-1
lisnam<-NULL
Nplant<-51 # Fixed number of resources
for (spe in spen){
  
  for (i in 1:nsim){
    
    MaEv[[i]]<-gen_even(Nbee,Nplant, spe,samp=F, make="spread")
    
  }
  lisnam[[counter]]<-Nplant
  
  lisEv[[counter]]<-MaEv
  counter<-counter+1
}




######################### 2. CALCULATIONS ######################################


##  Generalization 

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


## Spearman correlations between generalization indices and the
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
