################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### A tutorial with an example to help you use the alpha_PDI, genfun, and wcfun
### functions. If you want to aplly this analysis to your own data, follow the
### same formats of the data frames and vectors used in this example.


######################### 1. SETTINGS ##########################################


## 1.1 Clean the environment.
rm(list= ls())

## 1.2. Check the required packages, install them if necessary, and load them.
if(!require(zCompositions)){
  install.packages("zCompositions")
  library(zCompositions)
  }

if(!require(bipartite)){
  install.packages("bipartite")
  library(bipartite)
}

if(!require(corrgram)){
  install.packages("corrgram")
  library(corrgram)
}

## 1.3. Source the functions.
source("Code/alpha_PDI.R") # αPDI
source("Code/genfun.R") # other common generality indices
source("Code/wcfun.R") # generality of Pierotti et al. (2017)

## 1.4. Create a random data set or load your own data.
# In this example we create a random data set.
# If you prefer to build the object "data" using your own data, your data
# must correspond to a vector of resource use, which contains either raw
# frequencies (counts) or relative frequencies (proportions). It may also 
# correspond to a matrix with multiple consumers in the rows and multiple
# resources in the columns.
data <- matrix(round(runif(30, min=1, max = 100)), byrow = T, ncol=10)
row.names(data) <- c("Consumer1","Consumer2","Consumer3") #Example with 3 consumers

# Check the format of the random data created in this example.
data
class(data)
str(data)

# Visualize your data as smooth niche lines.
data2 <- as.data.frame(t(data))
plot(data2$Consumer1, type = "p",
     col = "white",
     xlab = "Resource",
     ylab = "Intensity of use",
     ylim = c(min(data2),max(data2)))
lines(spline(data2$Consumer1, n=1000), col = "blue")
lines(spline(data2$Consumer2, n=1000), col = "yellow")
lines(spline(data2$Consumer3, n=1000), col = "green")
legend(nrow(data2)-2, max(data2), c(colnames(data2)),
       lty = 1, col = c("blue", "yellow", "green"))

# The resource abundance distribution must be entered as a vector. It may
# contain raw frequencies (counts) or relative frequencies (proportions).
# If you prefer to build the object "abun" using your own data, your data
# must correspond to a vector of resource abundances.
abun <- round(runif(10, min=1, max = 100))

# Check the format of the random data created in this example.
abun
class(abun)
str(abun)


######################### 2. CALCULATIONS ######################################


## 2.1. Alpha PDI.

# Run the function alpha_PDI with the argument "corrected = F" to calculate the
# raw version of the index.

aPDIraw <- alpha_PDI(data, abun, corrected = F) # You get a value for each consumer
aPDIraw

# Use the function alpha_PDI with the argument "corrected = T" to calculate the
# raw version of the index together with its values corrected by the maximum 
# possible value given the number of interactions. The corrected version of
# αPDI cannot be calculated if your data are not counts (integers).

aPDIcor <- alpha_PDI(data, abun, corrected = T) 
aPDIcor

# For each consumer, you get a list with three vectors:
# (1) $corrected_aPDI: The value of αPDI corrected by the maximum possible
#                      value ($raw_aPDI/$max_aPDI).
# (2) $raw_aPDI: Uncorrected value of αPDI.
# (3) $max_aPDI: Maximum possible value of αPDI given the total 
#                number of interactions made by the consumer.

# If αPDI < 0.5 the consumer is a specialist.
# If αPDI > 0.5 the consumer is a generalist.


## 2.2. Other common indices of generality

indices <- genfun(data, abun)
indices

## 2.3. Wc: Generality of Pierotti et al. (2017)

wcfun(data, abun) #Only works for matrices. May not work for some matrices.


######################### 3. CORRELATIONS ######################################


## Check the correlations between all generality indices calculated in
## this example.
correlations <- cbind(indices, aPDIcor$corrected_aPDI)
names <- colnames(correlations)
colnames(correlations) <- c((names[1:(length(names)-1)]), "aPDI")
correlations <- correlations[,-1]
correlations

corrgram(correlations,
         cor.method = "pearson",
         lower.panel = panel.cor,
         upper.panel = panel.shade,
         diag.panel = panel.density)


######################### 4. REFERENCES ########################################


# Pierotti, M. E. R., Martín-Fernández, J. A., & Barceló-Vidal, C. (2017). The
# peril of proportions: robust niche indices for categorical data. Methods in
# Ecology and Evolution, 8(2), 223–231. 
# doi: https://doi.org/10.1111/2041-210X.12656
