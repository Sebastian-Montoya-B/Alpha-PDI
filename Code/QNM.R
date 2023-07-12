################################################################################
#### Ecological Synthesis Lab (SintECO): https://marcomellolab.wordpress.com

#### Authors: Sebastian Montoya-Bustamante, Carsten F. Dormann, 
####          Boris R. Krasnov, Marco A. R. Mello

#### See README for further info:
#### https://github.com/Sebastian-Montoya-B/Alpha-PDI#readme
################################################################################


## This script creates the functions used to generate the theoretical matrices
## used in our analyses.

### Most of this code was written by Fründ et al. (2016), who used it to create
### their quantitative niche model.
### We added some new parts and synthesized the code.


######################## 1. makeweb function ###################################


## A function to generate true preference vectors for a given specialization
## parameter (specpar), number of consumers (Nbee) and resources (Nplant), and
## a selected shape for the true generalization curve (nicheshape).
## Within this function, trait matching between consumers and resources is 
## random.
makeweb <- function(specpar = 1, Nbee=4, Nplant=4, nicheshape="normal",
                    planttrait_n="multi"){

  fun_pref <- function(traitdif){
    if (nicheshape == "normal") {
      prefs <- dnorm(traitdif,mean=0,sd=1/specpar)
    }
    if (nicheshape == "parabolic"){
      prefs <- 1-specpar*(traitdif)^2
      prefs[prefs < 0] <- 0
    }
    if (nicheshape == "power") {
      prefs <- (1-traitdif)^specpar
    }
    if (nicheshape == "logistic") {
      prefs <- plogis(traitdif-specpar,0,0.05*(1-specpar))
    }
    prefs
  }
  
  if (planttrait_n == "single"){
    planttrait <- runif(Nplant)
    web <- replicate(Nbee, fun_pref(abs(runif(1) - planttrait)))
  }
  if (planttrait_n == "multi"){
    web <- replicate(Nbee, fun_pref(abs(runif(1) - runif(Nplant))))
  }
  web <- web / matrix(colSums(web),nrow=Nplant,ncol=Nbee,byrow=TRUE) # standardize link weights to probability;
  web
}


######################## 2. get_skewabuns function #############################


## A function to generate resource abundance distributions based on a log-normal
## distribution. It excludes 0 and Inf values.
get_skewabuns <- function(myN, abun_mean=10, abun_sdlog=1.5){
  abun <- qlnorm(seq(0, 1, length.out=myN+2), log(abun_mean), abun_sdlog)[-c(1,myN+2)] 
  #abun <- sort(abun, decr=TRUE)
  abun<- sample(abun) # We added this line so the most abundant resource is not always the first
  abun * abun_mean / mean(abun)
}


######################## 3. make_currentweb function ###########################


## A function to generate the vectors (or matrices) of current patterns of 
## resource use. It requieres the true preference vector (web_p), the resource 
## abundance distribution (plantabun), and the consumer abundance distribution 
## (beeabun). 
make_currentweb <- function(web_p, plantabun, beeabun){
  web_relabun <- (plantabun %*% t(beeabun)) / mean(plantabun)
  web_p * web_relabun
}


######################## 4. sampleweb function #################################


## A function to simulate the process of sampling from a vector (or matrix) of 
## the current pattern of resource use. It requires the vector of current 
## pattern of resource use (web), and the number of observations (obsperbee).
sampleweb <- function(web,obsperbee=NULL,method='perbee'){
  
  if(method=='perbee') {
    if (length(obsperbee)==1){obsperbee <- rep(obsperbee,ncol(web))}
    if (length(obsperbee)!=ncol(web)){stop("length of obsperbee neither 1 nor matching number of species")}
    sampledweb <- sapply(1:ncol(web), function(j) {
      web.j <- web[,j]
      table(sample(factor(1:length(web.j)),obsperbee[j],prob=web.j,replace=TRUE))
      
    })
  }
  
  if(method=='perweb'){
    if (length(obsperbee)!=1){stop("obsperbee must have length 1 with method 'perweb'")}
    Nobs <- round(obsperbee * ncol(web))
    sampledweb.vect <- sample(as.factor(1:length(web)), size=Nobs, prob=as.numeric(web), replace=TRUE)
    sampledweb <- matrix(as.numeric(table(sampledweb.vect)), nrow=nrow(web))
  }
  
  if(method=='rarefy'){
    if(is.null(obsperbee)) obsperbee <- min(colSums(web))
    if(obsperbee<1){stop("cannot rarefy to less than 1 observation")}
    sampledweb <- sapply(1:ncol(web), function(j) {
      web.j <- web[,j]
      table(sample(rep(factor(1:length(web.j)),web.j),obsperbee,replace=FALSE))
    })
  }
  dimnames(sampledweb) <- dimnames(web)
  sampledweb
}


######################## 5. makeweb2 function ##################################


## A function to generate true preference vectors for a given specialization
## parameter (specpar), number of consumers (con) and resources (res), for a
## normal-shaped true generalization curve. Unlike the makeweb function that 
## uses random trait matching, makeweb2 generates vectors of trait matching 
## between consumers and resources with an an equidistant sequence from -1 to 1 
## of length equal to the number of potential resources. When using an odd 
## number of resources, there will always be a resource with highest trait
## matching (0).

makeweb2 <- function(specpar = 1, con=2, res=4){
  traitdif<-replicate(con, seq(-1,1,length=res))
  web <- dnorm(traitdif,mean=0,sd=1/specpar)
  web <- web / matrix(colSums(web),nrow=res,ncol=con,byrow=TRUE)
  web
}


######################## 6. gen_even function ##################################


## It compiles all the previous generated functions to generate vectors 
## (or matrices) of resource use assuming an even resource abundance 
## distribution. The output is a list with the resource abundance distribution
## used ($res_abnun), the true preference vector ($preference), and the vector 
## of the current pattern of resource use $current. If samp=T, the output list
## will also contain observed patterns of resource use vectors of two sizes, 
## $small and $large, based on minsamp and maxsamp number of observations, 
## respectively.

gen_even<-function(Nbee, Nplant,spe,samp=F,minsamp=5,maxsamp=1000, make=c("random", "spread")){
  
  # for even webs:
  beeabun <- rep(10, Nbee) 
  plantabun <- rep(10, Nplant) 
  
  
  if (make=="random"){
    web_p <- makeweb(specpar=spe, Nbee, Nplant) 
  } else if (make=="spread"){
    web_p <- makeweb2(specpar=spe, Nbee, Nplant)
  }
  
  web_current <- make_currentweb(web_p, plantabun=plantabun, beeabun=beeabun)
  
  if (samp==T){
    web_smallsamp <- sampleweb(web_current, obsperbee=minsamp, method='perweb')
    web_largesamp <- sampleweb(web_current, obsperbee=maxsamp, method='perweb')
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current, small=web_smallsamp, large=web_largesamp)
  } else{
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current)
  }
  
  return(sim_spec)
}

######################## 7. gen_uneven function ################################


## It compiles all the previous generated functions to generate vectors 
## (or matrices) of resource use assuming an uneven resource abundance 
## distribution. The output is a list with the resource abundance distribution
## used ($res_abnun), the true preference vector ($preference), and the vector 
## of the current pattern of resource use $current. If samp=T, the output list
## will also contain sampled vectors of two sizes, $small and $large, based on
## minsamp and maxsamp number of observations, respectively.
gen_uneven<-function(Nbee, Nplant, spe,samp=F, minsamp=5,maxsamp=1000, make=c("random", "spread")){
  
  beeabun <- rep(10, Nbee)
  plantabun <- get_skewabuns(Nplant)
  if (make=="random"){
    web_p <- makeweb(specpar=spe, Nbee, Nplant) 
  } else if (make=="spread"){
    web_p <- makeweb2(specpar=spe, Nbee, Nplant)
  }
  
  web_current <- make_currentweb(web_p, plantabun=plantabun, beeabun=beeabun) 
  
  if (samp==T){
    web_smallsamp <- sampleweb(web_current, obsperbee=minsamp, method='perweb')
    web_largesamp <- sampleweb(web_current, obsperbee=maxsamp, method='perweb')
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current, small=web_smallsamp, large=web_largesamp)
  } else{
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current)
  }
}


######################## 8. gen_even2 function ##################################


## It compiles all the previous generated functions to generate vectors 
## (or matrices) of resource use assuming an even resource abundance 
## distribution. The output is a list with the resource abundance distribution
## used ($res_abnun), the true preference vector ($preference), and the vector 
## of the current pattern of resource use $current. If samp=T, the output list
## will also contain observed patterns of resource use vectors of two sizes, 
## $small and $large, based on minsamp and maxsamp number of observations, 
## respectively.

gen_even2<-function(Nbee, Nplant,spe,samp=F,minsamp=5,maxsamp=1000, make=c("random", "spread")){
  
  # for even webs:
  beeabun <- rep(10, Nbee) 
  plantabun <- rep(10, Nplant) 
  
  
  if (make=="random"){
    web_p <- makeweb(specpar=spe, Nbee, Nplant) 
  } else if (make=="spread"){
    web_p <- makeweb2(specpar=spe, Nbee, Nplant)
  }
  
  web_current <- make_currentweb(web_p, plantabun=plantabun, beeabun=beeabun)
  
  if (samp==T){
    web_smallsamp1 <- sampleweb(web_current, obsperbee=minsamp[1], method='perweb')
    web_smallsamp2 <- sampleweb(web_current, obsperbee=minsamp[2], method='perweb')
    web_smallsamp3 <- sampleweb(web_current, obsperbee=minsamp[3], method='perweb')
    web_smallsamp4 <- sampleweb(web_current, obsperbee=minsamp[4], method='perweb')
    web_largesamp <- sampleweb(web_current, obsperbee=maxsamp, method='perweb')
    sim_spec<-list(res_abun=plantabun, preference=web_p, 
                   current=web_current, small10=web_smallsamp1, 
                   small50=web_smallsamp2, small80=web_smallsamp3, 
                   small100=web_smallsamp4, 
                   large=web_largesamp)
  } else{
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current)
  }
  
  return(sim_spec)
}



######################## 9. gen_uneven2 function ################################


## It compiles all the previous generated functions to generate vectors 
## (or matrices) of resource use assuming an uneven resource abundance 
## distribution. The output is a list with the resource abundance distribution
## used ($res_abnun), the true preference vector ($preference), and the vector 
## of the current pattern of resource use $current. If samp=T, the output list
## will also contain sampled vectors of two sizes, $small and $large, based on
## minsamp and maxsamp number of observations, respectively.
gen_uneven2<-function(Nbee, Nplant, spe,samp=F, minsamp=c(10,50,80,100),maxsamp=1000, make=c("random", "spread")){
  
  beeabun <- rep(10, Nbee)
  plantabun <- get_skewabuns(Nplant)
  if (make=="random"){
    web_p <- makeweb(specpar=spe, Nbee, Nplant) 
  } else if (make=="spread"){
    web_p <- makeweb2(specpar=spe, Nbee, Nplant)
  }
  
  web_current <- make_currentweb(web_p, plantabun=plantabun, beeabun=beeabun) 
  
  if (samp==T){
    web_smallsamp1 <- sampleweb(web_current, obsperbee=minsamp[1], method='perweb')
    web_smallsamp2 <- sampleweb(web_current, obsperbee=minsamp[2], method='perweb')
    web_smallsamp3 <- sampleweb(web_current, obsperbee=minsamp[3], method='perweb')
    web_smallsamp4 <- sampleweb(web_current, obsperbee=minsamp[4], method='perweb')
    web_largesamp <- sampleweb(web_current, obsperbee=maxsamp, method='perweb')
    sim_spec<-list(res_abun=plantabun, preference=web_p, 
                   current=web_current, small1=web_smallsamp1, 
                   small2=web_smallsamp2, small3=web_smallsamp3, 
                   small4=web_smallsamp4, 
                   large=web_largesamp)
  } else{
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current)
  }
  
  return(sim_spec)
}


######################## 10. gen_uneven3 function ################################


## Almost the same as gen_uneven2. It was exclusively created to make networks so
## we can calculate the Wc index. Instead of generating new resource abundance
## distributions, it uses existing vectors of abundance (res_abun argument).

gen_uneven3<-function(Nbee,Nplant,  spe,samp=F, minsamp=c(10,50,80,100),maxsamp=1000, make=c("random", "spread"),res_abun){
  
  beeabun <- rep(10, Nbee)
  plantabun <- res_abun
  #plantabun <- get_skewabuns(Nplant)
  if (make=="random"){
    web_p <- makeweb(specpar=spe, Nbee, Nplant) 
  } else if (make=="spread"){
    web_p <- makeweb2(specpar=spe, Nbee, Nplant)
  }
  
  web_current <- make_currentweb(web_p, plantabun=plantabun, beeabun=beeabun) 
  
  if (samp==T){
    web_smallsamp1 <- sampleweb(web_current, obsperbee=minsamp[1], method='perweb')
    web_smallsamp2 <- sampleweb(web_current, obsperbee=minsamp[2], method='perweb')
    web_smallsamp3 <- sampleweb(web_current, obsperbee=minsamp[3], method='perweb')
    web_smallsamp4 <- sampleweb(web_current, obsperbee=minsamp[4], method='perweb')
    web_largesamp <- sampleweb(web_current, obsperbee=maxsamp, method='perweb')
    sim_spec<-list(res_abun=plantabun, preference=web_p, 
                   current=web_current, small10=web_smallsamp1, 
                   small50=web_smallsamp2, small80=web_smallsamp3, 
                   small100=web_smallsamp4, 
                   large=web_largesamp)
  } else{
    sim_spec<-list(res_abun=plantabun, preference=web_p, current=web_current)
  }
  
  return(sim_spec)
}

######################## REFERENCES ############################################


## Fründ, J., Mccann, K. S., & Williams, N. M. (2016). Sampling bias is a 
## challenge for quantifying specialization and network structure: lessons 
## from a quantitative niche model. Oikos, 502–513. 
## doi: https://doi.org/10.1111/oik.02256