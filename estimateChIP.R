##############################################################
#############DIRICHLET PRIOR##################################

dirichlet<-function(x, alfa){
   return(ddirichlet(c(x,1-x),c(wt*alfa,wt*(1-alfa))))
}


prior_dir<-function(x){
   prob<-x[1:sample_n]
   max_x<-matrix(x[(sample_n+1):length(x)], ncol=2)
   dir<-0
   for(i in 1:sample_n){
      dir=dir+log(dirichlet(prob[i], p[i,1]))
  }
   return(-dir)
}

##############################################################
#############LIKELIHOOD#######################################

likelihood<-function(x){
   prob<-x[1:sample_n]
   prob<-as.numeric(prob)
   max_x<-matrix(x[(sample_n+1):length(x)], ncol=celltype_n)
   likelihood<-0
      for(j in 1:ncol(data)){
         likelihood<-likelihood+data[,j]*log(scale[j]*(as.numeric(prob[j])*as.numeric(max_x[,1])+as.numeric(1-prob[j])*as.numeric(max_x[,2])))-scale[j]*as.numeric(prob[j])*as.numeric   (max_x[,1])-scale[j]*as.numeric(1-prob[j])*as.numeric(max_x[,2])-lgamma(data[,j]+1)
     }

  return(-sum(likelihood))
}
##############################################################
#############POSTERIOR########################################

posterior<-function(x){
   post<-0
   post<-prior_dir(x)+likelihood(x)
   return(post)
}
##############################################################
#############GRADIENT OF POSTERIOR############################
gradient_dir<-function(x,s){
   alfa<-p[,1]
   prob<-x[1:sample_n]
   prob<-as.numeric(prob)
   max_x<-matrix(x[(sample_n+1):length(x)], ncol=2)
   dir<-(wt*alfa-1)/prob-(wt*(1-alfa)-1)/(1-prob)
   for(j in 1:ncol(data)){
      dir[j]=dir[j]+sum((data[,j]/(scale[j]*prob[j]*max_x[,1]+scale[j]*(1-prob[j])*max_x[,2])*scale[j]*(max_x[,1]-max_x[,2])-scale[j]*max_x[,1]+scale[j]*max_x[,2]))
   }
   return(-dir)
}

gradient_poisson<-function(x,i,t){
  if(celltype_n==t){
     prob<-1-x[1:sample_n]
   }else{
     prob<-x[1:sample_n]
   }
   prob<-as.numeric(prob)
   max_x<-matrix(x[(sample_n+1):length(x)], ncol=2)
   pois<-0
   for(j in 1:sample_n){
      pois=pois+data[i,j]/(scale[j]*prob[j]*max_x[i,t]+scale[j]*(1-prob[j])*max_x[i,-t])*scale[j]*prob[j]-scale[j]*prob[j]
   }
   return(-pois)
}

gradient<-function(x){
   as.numeric(c(gradient_dir(x),gradient_poisson(x,,1), gradient_poisson(x,,2)))
}


##############################################################
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

p<-read.table(args[3])
celltype_n<-ncol(p)
p0<-matrix(ncol=celltype_n, nrow=sample_n)
for(i in 1:nrow(p)){
   for(j in 1:ncol(p)){
      p0[i,j]<-as.numeric(p[i,j])
   }
}

apu<-c();mean<-matrix(ncol=2, nrow=nrow(data))
data_temp<-as.matrix(t(t(data)/scale))
for (i in 1:nrow(data)){
 apu<- nnls(as.matrix(p),as.numeric(data_temp[i,]))
  for (j in 1:celltype_n){
     mean[i,j]<-apu$x[j]
 }
}


tmpfX=Inf;tmpX<-c();ai<-c();bi<-c();ci<-c()
lower<-as.numeric(args[5]);upper<-as.numeric(args[6])
maxIter<-as.numeric(args[7])
for(m in 1:maxIter){
  ai<-c()
  for (i in 1:celltype_n){
   ai<-c(ai,mean[,i]+runif(nrow(data), -30,50))
  }
  ai[which(ai<lower)]<-lower;ai[which(ai>upper)]<-upper
  ci<-p[,1]+runif(sample_n, -0.2, 0.2)
  ci[which(ci<0.05)]<-0.05; ci[which(ci>0.95)]<-0.95
  start<-c(ci,ai)
  lowerALL=lowerALL=c(rep(0.01, sample_n), rep(lower, nrow(data)*2));upperALL=c(rep(0.99, sample_n), rep(upper, nrow(data)*2))
  mini<-optim(start,posterior, gr=gradient,upper=upperALL, lower=lowerALL,method="L-BFGS-B", control=list(maxit=1000000))
  variables<-mini$par;fX<-mini$value
  if(as.numeric(tmpfX) > as.numeric(fX)){
     tmpX = variables
     tmpfX = fX
  }

}

prob<-tmpX[1:(sample_n*(celltype_n-1))]; max_x<-tmpX[-(1:(sample_n*(celltype_n-1)))]
write.table(prob, args[8]); write.table(max_x, args[9])

