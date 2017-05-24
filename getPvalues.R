######################################

##Copyright (c) 2008,2009,2010,2011 Hyunjin Shin, Tao Liu <taoliu@jimmy.harvard.edu>

#This code is free software; you can redistribute it and/or modify it
#under the terms of the BSD License (see the file COPYING included with
#the distribution).

#This function is originally a Python function converted to R.

LSTEP = 200
EXPTHRES = exp(LSTEP)
EXPSTEP  = exp(-LSTEP)

poissonSmall<-function(k,a){

    if(k < 0){
        return(1) 
    }                       # special cases
    seur = exp( -a )

    if(floor(k)>1){
       for(i in 1:floor(k)){
           last = seur
           seur = last * a / i
       }
    }

    cdf = 0
    i = k+1

    while(seur>0){
        last = seur
        seur = last * a / i
        cdf =cdf+ seur
        i=i+1
    }
    return(cdf)
}

poissonLarge<-function(k,a){
    
    if(k < 0){
        return(1)       
    }                 # special cases
    num_parts = floor(a/LSTEP)
    last_part = a-floor(a/LSTEP)*LSTEP
    lastexp = exp(-last_part)
    seur = EXPSTEP
    num_parts =num_parts- 1

    for(i in 1:floor(k)){
        last = seur
        seur = last * a / i
        if(seur > EXPTHRES){
           if(num_parts>=1){
               seur = seur*EXPSTEP
               num_parts = num_parts-1
           }
           else{
               cdf = cdf*lastexp
               lastexp = 1
           }
       }
    }
    cdf = 0
    i = k+1
    while(seur >0){
        last = seur
        seur = last * a / i
        cdf =cdf+seur
        i=i+1
        if(seur > EXPTHRES || cdf > EXPTHRES){
           if(num_parts>=1){
               cdf =cdf* EXPSTEP
               seur= seur*EXPSTEP
               num_parts =num_parts- 1
           }
           else{
               cdf =cdf* lastexp
               lastexp = 1
           }
        }
    }
    if(floor(num_parts)>1){
       for(i in 0:(floor(num_parts)-1)){
           cdf =cdf* EXPSTEP
       }
    }
    cdf =cdf* lastexp
    return(cdf) 
}

###################################################
poisson<-function(k,a){

   if(k<0.1){
      return(1)
    }
   if(a>700){
      poissonLarge(k,a)
   }
   else{
      poissonSmall(k,a)
   }
}

#################################################
test<-function(signal, input, scaleChip, scaleInput){

if(scaleChip>=scaleInput){
  scale<-scaleChip/scaleInput
  for(m in 1:3){
    input[[m]]<-scale*input[[m]]
   }
}
if(scaleInput>scaleChip){
  scale<-scaleInput/scaleChip
  signal<-scale*signal
}

pvalues<-matrix(ncol=ncol(signal), nrow=nrow(signal))

for(i in 1:nrow(signal)){
   a<-which(peaks[,4]==strsplit(rownames(signal)[i], split="\t")[[1]][1])
   length<-peaks[a,3]-peaks[a,2]
   for(j in 1:ncol(signal)){
      lambda<-max(scaleInput/2.7e9, input[[1]][i,j]/1000, input[[2]][i,j]/5000, input[[3]][i,j]/10000) 
      pvalues[i,j]<-poisson(signal[i,j],lambda*length)  

 
   }
}

return(pvalues)
}

#####################################################
args<-commandArgs(TRUE)

signal<-read.table(args[1])

input<-c()
for(m in 1:3){
  input[[m]]<-read.table(args[m+1])
}

peaks<-read.table(args[5], header=T)
data<-read.table(args[6])
include<-c()
for(i in 1:nrow(data)){
   include[i]<-which(peaks[,4]==as.character(data[i,1]))
}
peaks<-peaks[include,]


cells<-2
signal<-matrix(as.matrix(signal), ncol=cells)
rownames(signal)<-data[,1]
for(m in 1:3){
   input[[m]]<-matrix(as.matrix(input[[m]]), ncol=2)
   rownames(input[[m]])<-data[,1]
}



#scale of Signal
scale<-read.table(args[7]);scale<-c(scale[[1]])
sizeChip<-max(scale)

#scale of Input
scale<-read.table(args[8]);scale<-c(scale[[1]])
sizeInput<-max(scale)


pvalues<-test(signal, input, sizeChip, sizeInput)
pvalues<-cbind(peaks,pvalues)
colnames(pvalues)<-c("chromosome", "start", "end", "name", paste("pvalue_cell", seq(1:cells), sep="")) 
write.table(pvalues, args[9], quote=F, row.names=F)
