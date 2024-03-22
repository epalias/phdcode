#create a d*d Walsh-Hadamard matrix
WHmat=function(d){
  A=matrix(NA,d,d)
  for(i in 1:d){
    for(j in 1:d){
      A[i,j]=(-1)^((as.numeric(intToBits(i-1))%*%as.numeric(intToBits(j-1)))%%2)/sqrt(d)
    }
  }
  A
}

d=100 #ambient dimension
n=2000 #number of instances
A=diag(sqrt(c(Surrogate::RandVec(s=100,n=100)[[1]]))) #generate eigenvalues with fixed stable rank
A=diag((d:1)/d) #matrix that determines the ellipsoid
set=uniformly::runif_on_sphere(n,d,1)%*%A #ellipsoid; set a seed for this and the below
w=rnorm(d) #vector that determines linear separability; set a seed for this and the above
label=as.factor(sign(set%*%w)) #determine the labels by linear boundary
set=set[order(label),] #order set by labels
label=label[order(label)] #order labels
s=sort(sample(1:nrow(set),.8*nrow(set))) #training set
t=10 #number of random compressions per k
e=0 #standard deviation of noise

#embed set to d dimensions
set=set%*%t(pracma::randortho(d)[,1:ncol(set)])
#add Gaussian noise to set
set=set+mvtnorm::rmvnorm(nrow(set),sigma=diag(e,ncol(set)))

#iterations
for(k in (1:10)*5){ #projection dimension
  a1=rep(NA,t) #test errors of R1 with metric
  a2=rep(NA,t) #test errors of R2 with metric
  a3=rep(NA,t) #test errors of R3 with metric
  a4=rep(NA,t) #test errors of R4 with metric
  a5=rep(NA,t) #test errors of R5 with metric
  b1=rep(NA,t) #test errors of R1 without metric
  b2=rep(NA,t) #test errors of R2 without metric
  b3=rep(NA,t) #test errors of R3 without metric
  b4=rep(NA,t) #test errors of R4 without metric
  b5=rep(NA,t) #test errors of R5 without metric
  for(i in 1:t){ #number random compressions
    R1=matrix(rnorm(ncol(set)*k),k)/sqrt(k) #Gaussian random matrix (Dasgupta and Gupta, 2003)
    R2=pracma::randortho(ncol(set))[1:k,] #random orthogonal
    R3=matrix(sample(-1:1,ncol(set)*k,T,c(1/6,2/3,1/6)),k) #JL with binary coins II (Achlioptas, 2003)
    R4=min(log(nrow(set))^2/ncol(set),1); R4=(matrix(sample(0:1,ncol(set)*k,T,c(R4,1-R4)),k)*matrix(rnorm(ncol(set)*k),k))%*%WHmat(ncol(set))%*%diag(sample(c(-1,1),ncol(set),T)) #fast JL transform (Ailon and Chazelle, 2009)
    R5=rnorm(ncol(set)); R5=pracma::Toeplitz(R5[c(1,ncol(set):(ncol(set)-k+2))],R5)%*%diag(sample(c(-1,1),ncol(set),T))/sqrt(k) #circulant transform (Hinrichs and Vybíral, 2011) and (Freksen and Larsen, 2020)
    x1=set%*%t(R1) #compress original data set with R1
    xtrain1=x1[s,] #training set of R1
    xtest1=x1[-s,] #test set of R1
    x2=set%*%t(R2) #compress original data set with R2
    xtrain2=x2[s,] #training set of R2
    xtest2=x2[-s,] #test set of R2
    x3=set%*%t(R3) #compress original data set with R3
    xtrain3=x3[s,] #training set of R3
    xtest3=x3[-s,] #test set of R3
    x4=set%*%t(R4) #compress original data set with R4
    xtrain4=x4[s,] #training set of R4
    xtest4=x4[-s,] #test set of R4
    x5=set%*%t(R5) #compress original data set with R5
    xtrain5=x5[s,] #training set of R5
    xtest5=x5[-s,] #test set of R5
    ytrain=label[s] #training labels
    ytest=label[-s] #test labels
    M1=mlpack::lmnn(data.frame(xtrain1,ytrain))$output #train metric on the R1 sample
    a1[i]=mean(class::knn(xtrain1%*%t(M1),xtest1%*%t(M1),ytrain)!=ytest) #run 1-NN on the R1 set with metric
    b1[i]=mean(class::knn(xtrain1,xtest1,ytrain)!=ytest) #run 1-NN on the R1 set without metric
    M2=mlpack::lmnn(data.frame(xtrain2,ytrain))$output #train metric on the R2 sample
    a2[i]=mean(class::knn(xtrain2%*%t(M2),xtest2%*%t(M2),ytrain)!=ytest) #run 1-NN on the R2 set with metric
    b2[i]=mean(class::knn(xtrain2,xtest2,ytrain)!=ytest) #run 1-NN on the R2 set without metric
    M3=mlpack::lmnn(data.frame(xtrain3,ytrain))$output #train metric on the R3 sample
    a3[i]=mean(class::knn(xtrain3%*%t(M3),xtest3%*%t(M3),ytrain)!=ytest) #run 1-NN on the R3 set with metric
    b3[i]=mean(class::knn(xtrain3,xtest3,ytrain)!=ytest) #run 1-NN on the R3 set without metric
    M4=mlpack::lmnn(data.frame(xtrain4,ytrain))$output #train metric on the R4 sample
    a4[i]=mean(class::knn(xtrain4%*%t(M4),xtest4%*%t(M4),ytrain)!=ytest) #run 1-NN on the R4 set with metric
    b4[i]=mean(class::knn(xtrain4,xtest4,ytrain)!=ytest) #run 1-NN on the R4 set without metric
    M5=mlpack::lmnn(data.frame(xtrain5,ytrain))$output #train metric on the R5 sample
    a5[i]=mean(class::knn(xtrain5%*%t(M5),xtest5%*%t(M5),ytrain)!=ytest) #run 1-NN on the R5 set with metric
    b5[i]=mean(class::knn(xtrain5,xtest5,ytrain)!=ytest) #run 1-NN on the R5 set without metric
  } #end of i loop
  print(k) #print which iteration we are
  write.table(cbind(d=ncol(set),n=n,set=1,k=k,sigma=e,
                    am1=mean(a1),asd1=sd(a1),bm1=mean(b1),bsd1=sd(b1),
                    am2=mean(a2),asd2=sd(a2),bm2=mean(b2),bsd2=sd(b2),
                    am3=mean(a3),asd3=sd(a3),bm3=mean(b3),bsd3=sd(b3),
                    am4=mean(a4),asd4=sd(a4),bm4=mean(b4),bsd4=sd(b4),
                    am5=mean(a5),asd5=sd(a5),bm5=mean(b5),bsd5=sd(b5)
  ),
  row.names=F,col.names=F,quote=F,append=T,file="/rds/homes/e/exp093/Workfiles/results.txt") #write the results
} #end of k loop

###find the k nearest neighbours of same labels (Weinberger and Saul (2009), Section 3.1)
#x is the set whose rows are the instances sorted by labels; for this use x=x[order(y),];y=y[order(y)]
#y is the set of labels in class 'factor'
#k is how many nearest neighbours to consider

KNN=function(x,y,k){
  a=NULL
  b=c(0,as.vector(table(y)))
  for(i in 1:(length(b)-1)){
    a=rbind(a,dbscan::kNN(x[y==levels(y)[i],],k)$id+b[i])
  }
  a
}

###loss function of metric (Weinberger and Saul (2009), Section 3.2)

#x is the set whose rows are the instances
#y is the set of labels in class 'factor'
#k is how many nearest neighbours to consider
#M is the Mahalanobis metric; not necessarily positive-definite
#r is the regularisation parameter

metloss=function(x,y,k,M,r=.5){
  x=x%*%t(M)
  KN2=KNN(x,y,k)
  l=0
  for(i in 1:nrow(x)){
    for(j in KN2[i,]){
      sim=(x[i,]-x[j,])%*%(x[i,]-x[j,])
      l=l+(1-r)*sim
      dis=c(apply(sweep(x[y!=y[i],],2,x[i,])^2,1,sum)) #dissimilar distances
      for(m in 1:length(dis)){
        l=l+r*max(1+sim-dis[m],0)
      }
    }
  }
  l
}

#run metric learning after fixing the set

testerror=NULL
R=matrix(rnorm(20*100),ncol=20)/20
x=set%*%R #compress instances
for(t in (1:8)*200){
  M=mlpack::lmnn(data.frame(x[sample(s,t),],label[sample(s,t)]))$output
  M=M/svd(M)$d[1]
  testerror=c(testerror,metloss(x[-s,],label[-s],1,M))
}
testerror
