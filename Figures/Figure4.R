################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script reproduces Figure 4, and Figures S10 to S12 of the Appendix S1.


######################### 1. SETTINGS ##########################################

## Clean the environment.
rm(list= ls())

## Check the required packages, install them if necessary, and load them.

if(!require(purrr)){
  install.packages("purrr")
  library(purrr)
}

if(!require(scales)){
  install.packages("scales")
  library(scales)
}

if(!require(igraph)){
  install.packages("igraph")
  library(igraph)
}

if(!require(Ternary)){
  install.packages("Ternary")
  library(Ternary)
}

## Source the functions.
source("Code/alpha_PDI.R")


## Load the empirical flea-mammal networks (matrices) 
## and its resource abundance distributions (vectors)

consumers<-readRDS("Data/Fleas.RDS")
resources<-readRDS("Data/resource_abundances.RDS")

######################### 2. CALCULATIONS ######################################


## 2.1. Calculating αPDI

results<-map2(consumers,resources, ~alpha_PDI(.x,.y,corrected=T)$corrected_aPDI)
length(results)

n.inter<-lapply(consumers, rowSums)
length(n.inter)

## 2.2. Proportion of generalists, specialists, and consumers 
##      with few observations in each network

suff.inter<-lapply(n.inter, function (x){which(x>=80)})

suff.res<-NULL
for (i in 1:length(results)){
  suff.res[[i]]<-results[[i]][suff.inter[[i]]]
}

por.gen<-NULL # generalists
por.spe<-NULL # specialists
por.left<-NULL # few observations
for (i in 1:length(suff.res)){
  # % of generalists
  por.gen[[i]]<-length(which(suff.res[[i]]>0.5))/length(results[[i]])
  # % of specialists
  por.spe[[i]]<-length(which(suff.res[[i]]<0.5))/length(results[[i]])
  # % of consumers without enough observations
  por.left[[i]]<-1-(por.gen[[i]]+por.spe[[i]])
}

obs.per<-cbind(Generalists=unlist(por.gen),Specialists=unlist(por.spe), 
               few_obs=unlist(por.left))
row.names(obs.per)<-names(results)
obs.per


## Calcuations for the chosen network

gen_fleas<-alpha_PDI(consumers$Turkmenistan, 
                     resources$Turkmenistan)$corrected_aPDI


 #write.csv(cbind(aPDI=gen_fleas,obs=rowSums(consumers$Turkmenistan)),
  #         "Turkmenistan_aPDI.csv")
######################### 3. PLOTTING ##########################################


## 3.1. Preparing the figures
color_group<-function(gen_fleas){
  cor.cla<-rep(NA, length(gen_fleas))
  for (i in 1:length(gen_fleas)){
    ifelse(gen_fleas[i]<=0.1,cor.cla[i]<-1,
           ifelse (gen_fleas[i]>0.1 & gen_fleas[i]<=0.2, cor.cla[i]<-2,
                   ifelse(gen_fleas[i]>0.2 & gen_fleas[i]<=0.3, cor.cla[i]<-3, 
                          ifelse(gen_fleas[i]>0.3 & gen_fleas[i]<=0.4, cor.cla[i]<-4, 
                                 ifelse(gen_fleas[i]>0.4 & gen_fleas[i]<0.5, cor.cla[i]<-5,
                                        ifelse(gen_fleas[i]==0.5, cor.cla[i]<-6, 
                                               ifelse(gen_fleas[i]>0.5 & gen_fleas[i]<=0.6, cor.cla[i]<-7,
                                                      ifelse(gen_fleas[i]>0.6 & gen_fleas[i]<=0.7,cor.cla[i]<-8,
                                                             ifelse(gen_fleas[i]>0.7 & gen_fleas[i]<=0.8, cor.cla[i]<-9, 
                                                                    ifelse(gen_fleas[i]>0.8 & gen_fleas[i]<=0.9, cor.cla[i]<-10, 
                                                                           cor.cla[i]<-11))))))))))
    
    
  }
  return(cor.cla)
}


color_group(gen_fleas)
frame_col<-data.frame(consumer=rownames(consumers$Turkmenistan), 
                      aPDI=gen_fleas, color.group=color_group(gen_fleas))


imat<- graph_from_incidence_matrix(consumers$Turkmenistan, 
                                   directed = F,weighted = TRUE)
imat

algoritmo<-layout.lgl(imat)

## The myrhombus function was created by other author 
## see https://stackoverflow.com/questions/53886154/rhombus-shape-igraph-node

myrhombus <- function(coords, v = NULL, params) {
  vertex.color <- params("vertex", "color")
  if (length(vertex.color) != 1 && !is.null(v)) {
    vertex.color <- vertex.color[v]
  }
  vertex.size <- 1/200 * params("vertex", "size")
  if (length(vertex.size) != 1 && !is.null(v)) {
    vertex.size <- vertex.size[v]
  }
  
  symbols(x = coords[, 1], y = coords[, 2], bg = vertex.color,
          stars = cbind(1.2*vertex.size, vertex.size, 1.2*vertex.size, vertex.size),
          add = TRUE, inches = FALSE)
}

add_shape("rhombus", clip = shapes("circle")$clip,
          plot = myrhombus)




colores<-c("#00918d","#51a4a0","#7cb7b4","#a4cac8","#cadddc","#f1f1f1",
           "#f7d1c1","#f7b294","#f39267", "#ea723b","#df4f00")


V(imat)$color = V(imat)$type
V(imat)$color[which(V(imat)$color==FALSE)]<-colores[frame_col$color.group]
V(imat)$color = gsub(TRUE,"gray",V(imat)$color)
V(imat)$color[which(rowSums(consumers$Turkmenistan)<80)]<-alpha("#658c1f",0.7)## Consumers with less than 100 interactions
V(imat)$size<-V(imat)$type
V(imat)$size[which(V(imat)$size==FALSE)]<-rep(9,nrow(frame_col))
V(imat)$size[which(V(imat)$size==TRUE)]<-log(resources$Turkmenistan)
V(imat)$shape <- c("circle", "rhombus")[V(imat)$type+1]

E(imat)$width <- log(E(imat)$weight/5000)

colfunc <- colorRampPalette(c("#00918d", "#ffffff","#df4f00"))
legend_image <- as.raster(matrix(colfunc(20), nrow=1))

## 3.2. Exporting Figure 4

png(filename="Figures/Exported/Figure4.png", width=6000, height=4000, res=600)

layout(matrix(c(1,1,2,3,5,4),ncol=2, byrow=T), width = c(1,1),height = c(1,5,1))

### 3.2.1. First plot
par(mar=c(0,0,0,0))
plot(y=c(0,1),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '')
text(x=c(0.23,0.78), y=0.5, cex=2,
     labels=c("(A)\nPercentage of consumers in each network",
              "(B)\nTurkmenistan network"))

### 3.2.2. Second plot
par(mar=c(3,2,1,2))
TernaryPlot(alab="Consumers\nwith few observations \u2192", blab="Generalists \u2192", 
            clab="\u2190 Specialists", lab.col=c("#658c1f","#df4f00","#00918d"),
            grid.lines = 4, lab.cex=1.8, axis.cex=1.2,
            grid.minor.lines = 1, grid.minor.lty = "dotted")
TernaryPoints(obs.per[,c(3,1,2)], pch=16, cex=1.6,
              col="black", bg=alpha("black",0.7))
TernaryText(c(0.355,1-(0.355+0.295), 0.295) , labels=expression(""%<-%"(B)"), cex=1.3, col="red")


### 3.2.3. Third plot
par(las=1,mar=c(0,0,0,0))
plot(imat,vertex.label.cex=1,vertex.color = V(imat)$color, 
     vertex.size = V(imat)$size, edge.curved=.3,layout = algoritmo, rescale=T, 
     vertex.label.color="black",vertex.label=NA)


### 3.2.4. Fourth plot
par(mar=c(2,2,2,2))

plot(y=c(0,2),x=c(0,1),type = 'n', axes = F,xlab = '', ylab = '', 
     cex.main=2.2 ,main = expression(alpha*italic(PDI[Corrected])))
text(x=seq(0,1,l=3), y =1.5 , labels = seq(0,1,l=3), cex=1.5)
rasterImage(legend_image, 0, 0, 1,1)

dev.off()


## 3.3. Plotting Figure S10

png(filename="Figures/Exported/FigureS10.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))
layout.show(22)
for (i in 1:20){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  abline(h=80, lty=2)
  
}

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels="Number of observations", srt=90, cex=1.6)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI[Corrected])), cex=1.6)

dev.off()

## 3.4. Plotting Figure S11

png(filename="Figures/Exported/FigureS11.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))
layout.show(22)


for (i in 21:40){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  abline(h=80, lty=2)
  
}

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels="Number of observations", srt=90, cex=1.6)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI[Corrected])), cex=1.6)
dev.off()

## 3.5. Plotting Figure S12

png(filename="Figures/Exported/FigureS12.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(21,1,2,3,4,21,5,6,7,8,21,9,10,11,12,21,13,14,15,16,21,17,18,19,20,21,22,22,22,22),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))
layout.show(22)


for (i in 41:60){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  abline(h=80, lty=2)
  
}

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.55, labels="Number of observations", srt=90, cex=1.6)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.6, labels=expression(alpha*italic(PDI[Corrected])), cex=1.6)
dev.off()


## 3.6. Plotting Figure S13

png(filename="Figures/Exported/FigureS13.png", width=4000, height=4600, res=600)
par(las=1, mar=c(2,2,2,2))
layout(matrix(c(17,1,2,3,4,17,5,6,7,8,17,9,10,11,12,17,13,14,15,16,17,18,18,18,18,17,18,18,18,18),
              ncol=5, byrow = T), widths=c(10,45/2,45/2,45/2,45/2), heights=c(18,18,18,18,18,10))
layout.show(18)


for (i in 61:74){
  plot(results[[i]], n.inter[[i]], xlim=c(0,1), xlab="", ylab="", 
       main=names(results)[i], cex.main=0.9, xaxp=c(0,1,2))
  abline(v=0.5, lty=2)
  abline(h=80, lty=2)
  
}

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")


plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.65, labels="Number of observations", srt=90, cex=1.6)

plot(x=NULL, y=NULL, ann=F,xaxt="n", yaxt="n", xlim=c(0,1), ylim=c(0,1), type="n", bty="n")
text(x=0.5, y=0.95, labels=expression(alpha*italic(PDI[Corrected])), cex=1.6)
dev.off()

