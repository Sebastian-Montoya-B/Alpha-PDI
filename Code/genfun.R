#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, Boris R. Krasnov, Marco A. R. Mello

#### See README for further info

#### This script sets the genfun function for further calculations.

#### genfun function:

genfun<-function(data, abun){
  
  
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
  
  #####----Fort et al. (2016) g----
  
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
  

    summ<-as.data.frame(cbind(Consumer=eq,Bs=besest, "B'"=beprimast, W=wpetst, PS=pesest, FT=feitst,  "1-d'"=nd, gen=ngst))
  
  return(summ)
  
}  

