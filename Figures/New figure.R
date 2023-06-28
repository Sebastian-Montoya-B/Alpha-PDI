msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}


par(mar= c(4,2,2,1),las=1)
layout(matrix(seq(1,9), ncol=3, byrow=T))
layout.show(9)
###aPDI

error1<-data.frame(R5=lisObv$aPDI[1:l5]-lisExv$aPDI[1:l5], 
           R10=lisObv$aPDI[(l5+1):l10]-lisExv$aPDI[(l5+1):l10], 
           R50=lisObv$aPDI[(l10+1):l50]-lisExv$aPDI[(l10+1):l50])

boxplot(error1, ylim=c(-1,1))

msqe<-round(apply(as.matrix(error1), 2, msqerrorfun),2)

###Bs

errorBS<-data.frame("R5"=lisObv$Bs[1:l5]-lisExv$Bs[1:l5],
                  "R10"=lisObv$Bs[(l5+1):l10]-lisExv$Bs[(l5+1):l10],
                  "R50"=lisObv$Bs[(l10+1):l50]-lisExv$Bs[(l10+1):l50])

boxplot(errorBS, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorBS), 2, msqerrorfun),2)

###B

errorBprime<-data.frame(R5=lisObv$`B'`[1:l5]-lisExv$`B'`[1:l5], 
                        R10=lisObv$`B'`[(l5+1):l10]-lisExv$`B'`[(l5+1):l10],
                        R50=lisObv$`B'`[(l10+1):l50]-lisExv$`B'`[(l10+1):l50])

boxplot(errorBprime, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorBprime), 2, msqerrorfun),2)

###W

errorW<-data.frame(R5=lisObv$W[1:l5]-lisExv$W[1:l5], 
                   R10=lisObv$W[(l5+1):l10]-lisExv$W[(l5+1):l10], 
                   R50=lisObv$W[(l10+1):l50]-lisExv$W[(l10+1):l50])

boxplot(errorW, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorW), 2, msqerrorfun),2)

###PS

errorPS<-data.frame(R5=lisObv$PS[1:l5]-lisExv$PS[1:l5], 
                   R10=lisObv$PS[(l5+1):l10]-lisExv$PS[(l5+1):l10], 
                   R50=lisObv$PS[(l10+1):l50]-lisExv$PS[(l10+1):l50])

boxplot(errorPS, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorPS), 2, msqerrorfun),2)

###FT

errorFT<-data.frame(R5=lisObv$FT[1:l5]-lisExv$FT[1:l5], 
                    R10=lisObv$FT[(l5+1):l10]-lisExv$FT[(l5+1):l10], 
                    R50=lisObv$FT[(l10+1):l50]-lisExv$FT[(l10+1):l50])

boxplot(errorFT, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorFT), 2, msqerrorfun),2)

###1-d

errord<-data.frame(R5=lisObv$`1-d'`[1:l5]-lisExv$`1-d'`[1:l5], 
                    R10=lisObv$`1-d'`[(l5+1):l10]-lisExv$`1-d'`[(l5+1):l10], 
                    R50=lisObv$`1-d'`[(l10+1):l50]-lisExv$`1-d'`[(l10+1):l50])

boxplot(errord, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errord), 2, msqerrorfun),2)

###gen

errorgen<-data.frame(R5=lisObv$gen[1:l5]-lisExv$gen[1:l5], 
                   R10=lisObv$gen[(l5+1):l10]-lisExv$gen[(l5+1):l10], 
                   R50=lisObv$gen[(l10+1):l50]-lisExv$gen[(l10+1):l50])

boxplot(errorgen, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorgen), 2, msqerrorfun),2)

###Wc

errorWc<-data.frame(R5=lisObvwc[1:l5]-lisExvwc[1:l5], 
                     R10=lisObvwc[(l5+1):l10]-lisExvwc[(l5+1):l10], 
                     R50=lisObvwc[(l10+1):l50]-lisExvwc[(l10+1):l50])


msqe<-round(apply(as.matrix(errorWc), 2, msqerrorfun),2)

errorWc<-na.omit(errorWc) ## Wc may me unable to estimate the value for some consumers, and it will produce NAs

boxplot(errorWc, ylim=c(-1,1))

msqe<-round(apply(as.matrix(errorWc), 2, msqerrorfun),2)
