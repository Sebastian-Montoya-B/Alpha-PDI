msqerrorfun<-function(error){
  return(sum((error)^2)/length(error))
}

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
