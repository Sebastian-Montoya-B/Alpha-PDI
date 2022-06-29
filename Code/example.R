################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### An example to help you use the alpha_PDI, genfun, and wcfun functions.


######################### 1. SETTINGS ##########################################


## 1.1. Load the packages (not necessary for alpha_PDI)

if (!require(zCompositions)) install.packages("zCompositions")
if (!require(bipartite)) install.packages("bipartite")

library(bipartite) ## necessary for genfun
library(zCompositions) ## necessary for wcfun

## 1.2. Load the functions
source("Code/alpha_PDI.R") ## alpha PDI
source("Code/genfun.R") ## other common generality indices
source("Code/wcfun.R") ## generality of Pierotti et al. (2017)

## 1.3. Load your data set or create your own

# Your data must correspond to a vector of resource use (counts or proportions)
# or a matrix of with multiple consumers in the rows and resources in the columns.
data <- matrix(round(runif(30, min=1, max = 100)), byrow = T, ncol=10)
row.names(data) <- c("Consumer1","Consumer2","Consumer3") #This example contains 3 consumers
dim(data)
head(data)
tail(data)

# The resource abundance distribution must be imported as a vector.
# It can use raw frequencies or relative frequencies.
abun <- round(runif(10, min=1, max = 100))
abun


######################### 2. CALCULATIONS ######################################


## 2.1. Alpha PDI

# Use alpha_PDI with the argument "corrected = F" to calculate the raw version
# of the index

alpha_PDI(data, abun, corrected=F) ## You will obtain a value for each consumer

## Use alpha_PDI with the argument "corrected = T" to calculate 
## the raw version of the index along with its values corrected
## by the maximum possible value given the number of interactions.
## The corrected version of alpha PDI cannot be calculated if your
## data are not counts.

alpha_PDI(data, abun, corrected=T) 
## For each consumer, you will obtain a list of three vectors:
## (1) $corrected_aPDI: The value of alpha PDI corrected by the maximum possible
##                      value ($corrected_aPDI/$max_aPDI).
## (2) $raw_aPDI: Uncorected value of alpha PDI.
## (3) $max_aPDI: Maximum possible value of alpha PDI given the total 
##                number of interactions made by the consumer.

## If alpha PDI < 0.5 the consumer is a specialist
## If alpha PDI > 0.5 the consumer is a generalist

### 2.2. Other common indices of generality

genfun(data,abun)

### 2.2. Generality of Pierotti et al (2017)

wcfun(data,abun) #Only works for matrices. May not work for some matrices.

