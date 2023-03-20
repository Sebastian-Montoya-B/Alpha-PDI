################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### This script sets the genfun function for further calculations.


genfun<-function(data, abun){
  
  if (class(data)[1]!="matrix"){
    data<-matrix(data, ncol=length(data))
  }
  
  require(bipartite)
  
  ######---Schoener (1974) Bs---##### 
  
  Bs<-function(data, abun){
    
    p<-t(apply(data, 1, function(x) {x/sum(x)}))
    q<-abun/sum(abun)
    p2<-p*p
    q2<-q*q
    bs<-1/rowSums(t(apply(p2, 1, function (x){x/q2})))
    
    return(bs)
  }
  
  
  
  ######---HULBERT (1978) B'---##### 
  
  Bprime<-function(data, abun){
    
    p<-t(apply(data, 1, function(x) {x/sum(x)}))
    q<-abun/sum(abun)
    p2<-p*p
    bprime<-1/rowSums(t(apply(p2, 1, function (x){x/q})))
    
    return(bprime)
  }
  
  
  ######---Petraitis (1979) W---##### 
  
  WW<-function(data, abun){
    
    p<-t(apply(data, 1, function(x) {x/sum(x)}))
    q<-abun/sum(abun)
    
    lnw<-(-1)*rowSums(t(apply(p, 1, function(x){ifelse(x!=0,x*log(x/q),0)})))
    
    return(exp(lnw))
  }
  
  
  #####--- Feinsinger et al. (1981) PS ----
  
  PS<-function(data, abun){
    
    p<-t(apply(data, 1, function(x) {x/sum(x)}))
    q<-abun/sum(abun)
    ps<-1-0.5*t(apply(p, 1, function(x) {sum(abs(x-q))}))
    
    return(as.vector(ps))
  }
  
  #####---Smith (1982) FT---####
  
  FT<-function(data, abun){
    
    abun<-abun/sum(abun)
    data_prop<-t(apply(data, 1, function(x) {x/sum(x)}))
    ft<-t(apply(data_prop, 1, function(x) {sum(sqrt(x*abun))}))
    
    return(as.vector(ft))
  }
  
  #####----Fort et al. (2016) gen----
  
  generalization<-function(data, abun){
    gener<-1-bipartite::dfun(data, abuns=abun)$d/log(sum(data))
    
    return(gener)
  }
  


  #####----Compilation----
  
  nbs<-Bs(data, abun)
  nbp<-Bprime(data, abun)
  nww<-WW(data, abun)
  nps<-PS(data, abun)
  nft<-FT(data, abun)
  nd<-1-dfun(data, abuns = abun)$dprime
  ng<-generalization(data, abun) 

    
  minq<-min(abun/sum(abun))
  mind<-1-log(1/minq)/log(sum(data))
  
  #Normalization
  besest<-(nbs-minq^2)/(sum((abun/sum(abun))^2)-minq^2)
  beprimast<-(nbp-minq)/(1-minq)
  wpetst<-(nww-minq)/(1-minq)
  pesest<-(nps-minq)/(1-minq)
  feitst<-(nft-minq^(1/2))/(1-minq^(1/2))
  ngst<-(ng-mind)/(1-mind)
  eq<-seq(1,nrow(data))
  

    summ<-as.data.frame(cbind(Consumer=eq,Bs=besest, 
                              "B'"=beprimast, W=wpetst, PS=pesest, FT=feitst,
                              "1-d'"=nd, gen=ngst))
  
  return(summ)
  
}  


###################### REFERENCES ##############################################


# Feinsinger, P., Spears, E., & Poole, R. (1981). A Simple Measure of Niche. 
# Ecology, 62(1), 27–32. https://www.jstor.org/stable/1936664

# Fort, H., Vázquez, D. P., & Lan, B. L. (2016). Abundance and generalisation
# in mutualistic networks: Solving the chicken-and-egg dilemma. Ecology Letters,
# 19(1), 4–11. http://dx.doi.org/10.1111/ele.12535

# Hurlbert, S. (1978). The Measurement of Niche Overlap and Some Relatives.
# Ecology, 59(1), 67–77. Retrieved from https://www.jstor.org/stable/1936632

# Petraitis, P. S. (1979). Likelihood Measures of Niche Breadth and Overlap.
# Ecology, 60(4), 703–710. https://doi.org/10.2307/1936607

# Schoener, T. W. (1974). Some Methods for Calculating Competition Coefficients
# from Resource-Utilization Spectra. The American Naturalist, 108(961), 332–340.
# http://dx.doi.org/10.1086/282911

# Smith, E. P. (1982). Niche breadth, resource availability, and inference. 
# Ecology, 63(6), 1675–1681. http://dx.doi.org/10.2307/1940109