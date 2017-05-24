##############################################################
#############LIKELIHOOD#######################################
likelihood<-function(x){
   prob<-p
   prob<-as.numeric(prob)
   max_x<-matrix(x, ncol=2)
   likelihood<-0
      for(j in 1:ncol(data)){
         likelihood<-likelihood+data[,j]*log(scale[j]*(as.numeric(prob[j])*as.numeric(max_x[,1])+as.numeric(1-prob[j])*as.numeric(max_x[,2])))-scale[j]*as.numeric(prob[j])*as.numeric(max_x[,1])-scale[j]*as.numeric(1-prob[j])*as.numeric(max_x[,2])-lgamma(data[,j]+1)
     }

  return(-sum(likelihood))
}

##############################################################
#############POSTERIOR########################################

posterior<-function(x){
   post<-0
   post<-likelihood(x)
   return(post)
}

##############################################################
#############GRADIENT OF POSTERIOR############################
gradient_poisson<-function(x,i,t){
   prob<-p[,t]
   prob<-as.numeric(prob)
   max_x<-matrix(x, ncol=celltype_n)
   pois<-0
   for(j in 1:sample_n){
      pois=pois+data[i,j]/(scale[j]*prob[j]*max_x[i,t]+scale[j]*(1-prob[j])*max_x[i,-t])*scale[j]*prob[j]-scale[j]*prob[j]
   }
   return(-pois)
}



gradient<-function(x){
   as.numeric(c(gradient_poisson(x,,1), gradient_poisson(x,,2)))
}
##############################################
# install.packages("gtools",repos='http://cran.us.r-project.org')
# install.packages("nnls",repos='http://cran.us.r-project.org')
library(gtools)
library(nnls)
args<-commandArgs(TRUE)


wt=as.numeric(args[2])
scale<-read.table(args[4]);scale<-scale[[1]]/max(scale[[1]])
data<-read.table(args[1])
rownames(data)<-make.unique(as.character(data[,1]));data<-data[,-1]
sample_n<-ncol(data)


pEstimated<-as.matrix(read.table(args[3]))
celltype_n<-ncol(pEstimated)+1
p<-matrix(ncol=celltype_n, nrow=sample_n)
for(i in 1:nrow(pEstimated)){
      p[i,1]<-as.numeric(pEstimated[i])
      p[i,2]<-1-as.numeric(pEstimated[i])
}

temp<-c();mean<-matrix(ncol=2, nrow=nrow(data))
data_temp<-as.matrix(t(t(data)/scale))
for (i in 1:nrow(data)){
 temp<- nnls(as.matrix(p),as.numeric(data_temp[i,]))
  for (j in 1:celltype_n){
     mean[i,j]<-temp$x[j]
 }
}



tmpfX=Inf;tmpX<-c();ai<-c();bi<-c();ci<-c()
lower<-as.numeric(args[5]);upper<-as.numeric(args[6])
maxIter<-as.numeric(args[7])
for(m in 1:maxIter){
   for (i in 1:celltype_n){
     ai<-c(ai,mean[,i]+runif(nrow(data), -30,50))
    }
   ai[which(ai<lower)]<-lower;ai[which(ai>upper)]<-upper
   start<-ai
   lowerALL=rep(lower, length(start));upperALL=rep(upper, length(start))
   mini<-optim(start,posterior, gr=gradient,upper=upperALL, lower=lowerALL,method="L-BFGS-B", control=list(maxit=1000000))
   variables<-mini$par;fX<-mini$value
   if(as.numeric(tmpfX) > as.numeric(fX)){
        tmpX = variables
        tmpfX = fX
   }

}

write.table(tmpX, args[8])





