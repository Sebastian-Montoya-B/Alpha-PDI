################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


### Most of this code was written by Pierotti et al. (2017).
### We added some new parts and synthesized the code to a single function.
### This script sets the wcfun function for further calculations.


wcfun<-function(data, abun){
  require(zCompositions)
  
  abun<-abun/sum(abun)

  # Pierotti et al. (2017) Wc
  
  fndAit<-function(x,y){
    # x,y: two vector of counts
    clrx<-fnclr(x)
    clry<-fnclr(y)
    d<-sqrt(t(clrx-clry)%*%(clrx-clry))
    return(d)
    
    ## fnclr: clr-transformation
    fnclr<-function(x){
      lx<-log(x)
      clrx<-lx-mean(lx)
      return(clrx)
    }
    
  }
  
  # # fnclr: clr-transformation
  fnclr<-function(x){
    # x: vector of counts
    lx<-log(x)
    clrx<-lx-mean(lx)
    return(clrx)
  }
  
  # fnWc: Niche Width
  fnWc<-function(x,r=rep(1/length(x),length(x))){
    # x: vector of counts (D: number of components)
    # r: vector of resources proportion (Default=(1/D,...,1/D))
    if ((is.matrix(x))) stop("x must be a vector")
    if ((is.matrix(r))) stop("r must be a vector")
    if (length(x)!=length(r)) stop("The number of components in x and r do not agree")
    #if (any(x==0)) stop("Zeros in x are not allowed in function fnWc. Use package zCompositions to replace the zeros and rerun function fnWc")
    if (any(r==0)) stop("Zeros in r are not allowed in function fnWc. Use package zCompositions to replace the zeros and rerun function fnWc")
    
    
    w=exp(-fndAit(x,r)^2)
    return(w)
    
    
    # fndAit: Aitchison distance function
    fndAit<-function(x,y){
      # x,y: two vector of counts
      clrx<-fnclr(x)
      clry<-fnclr(y)
      d<-sqrt(t(clrx-clry)%*%(clrx-clry))
      return(d)
    }
    ## fnclr: clr-transformation
    fnclr<-function(x){
      lx<-log(x)
      clrx<-lx-mean(lx)
      return(clrx)
    }
    
  }
  
  if (any(data==0)){
    zdata<-cmultRepl(data, method="GBM")
  } else{
    zdata<-data
  }
  
  nwc<-apply(zdata, 1, function(x) {fnWc(x, abun)})
  
  return(nwc)
}


###################### REFERENCES ##############################################


# Pierotti, M. E. R., Martín-Fernández, J. A., & Barceló-Vidal, C. (2017). The
# peril of proportions: robust niche indices for categorical data. Methods in
# Ecology and Evolution, 8(2), 223–231. 
# doi: https://doi.org/10.1111/2041-210X.12656