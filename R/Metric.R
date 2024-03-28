#' Out-of-sample error rate for metric under different compression methods.
#' @param x matrix whose rows are the instances
#' @param y vector of labels
#' @param k projection dimension
#' @param s vector of indices to be used as training set
#' @param t number of random comrpessions to average
#' @param W which type of compression to use
#' 
#' @export
Metric=function(x,y,k=ncol(x),s=sort(sample(1:nrow(x),.8*nrow(x))),t=10,W=c('G','O','B','F','C')){
  ytrain=y[s] #training labels
  ytest=y[-s] #test labels
  if('G' %in% W){
    a1=rep(NA,t) #test errors of R1 with metric
    b1=rep(NA,t) #test errors of R1 with metric
    a2=rep(NA,t) #test errors of R2 with metric
    b2=rep(NA,t) #test errors of R2 with metric
    a3=rep(NA,t) #test errors of R3 with metric
    b3=rep(NA,t) #test errors of R3 with metric
    a4=rep(NA,t) #test errors of R4 with metric
    b4=rep(NA,t) #test errors of R4 with metric
    a5=rep(NA,t) #test errors of R5 with metric
    b5=rep(NA,t) #test errors of R5 with metric
    for(i in 1:t){
      R1=matrix(stats::rnorm(ncol(x)*k),k)/sqrt(k) #Gaussian random matrix (Dasgupta and Gupta, 2003)
      x1=x%*%t(R1) #compress original data set with R1
      xtrain1=x1[s,] #training set of R1
      xtest1=x1[-s,] #test set of R1
      M1=mlpack::lmnn(data.frame(xtrain1,ytrain))$output #train metric on the R1 sample
      a1[i]=mean(class::knn(xtrain1%*%t(M1),xtest1%*%t(M1),ytrain)!=ytest) #run 1-NN on the R1 set with metric
      b1[i]=mean(class::knn(xtrain1,xtest1,ytrain)!=ytest) #run 1-NN on the R1 set without metric
    }
  }
  if('O' %in% W){
    for(i in 1:t){
      R2=pracma::randortho(ncol(x))[1:k,] #random orthogonal
      x2=x%*%t(R2) #compress original data set with R2
      xtrain2=x2[s,] #training set of R2
      xtest2=x2[-s,] #test set of R2
      M2=mlpack::lmnn(data.frame(xtrain2,ytrain))$output #train metric on the R2 sample
      a2[i]=mean(class::knn(xtrain2%*%t(M2),xtest2%*%t(M2),ytrain)!=ytest) #run 1-NN on the R2 set with metric
      b2[i]=mean(class::knn(xtrain2,xtest2,ytrain)!=ytest) #run 1-NN on the R2 set without metric
    }
  }
  if('B' %in% W){
    for(i in 1:t){
      R3=matrix(sample(-1:1,ncol(x)*k,T,c(1/6,2/3,1/6)),k) #JL with binary coins II (Achlioptas, 2003)
      x3=x%*%t(R3) #compress original data set with R3
      xtrain3=x3[s,] #training set of R3
      xtest3=x3[-s,] #test set of R3
      M3=mlpack::lmnn(data.frame(xtrain3,ytrain))$output #train metric on the R3 sample
      a3[i]=mean(class::knn(xtrain3%*%t(M3),xtest3%*%t(M3),ytrain)!=ytest) #run 1-NN on the R3 set with metric
      b3[i]=mean(class::knn(xtrain3,xtest3,ytrain)!=ytest) #run 1-NN on the R3 set without metric
    }
  }
  if('F' %in% W){
    for(i in 1:t){
      R4=min(log(nrow(x))^2/ncol(x),1); R4=(matrix(sample(0:1,ncol(x)*k,T,c(R4,1-R4)),k)*matrix(stats::rnorm(ncol(x)*k),k))%*%WHmat(ncol(x))%*%diag(sample(c(-1,1),ncol(x),T)) #fast JL transform (Ailon and Chazelle, 2009)
      x4=x%*%t(R4) #compress original data set with R4
      xtrain4=x4[s,] #training set of R4
      xtest4=x4[-s,] #test set of R4
      M4=mlpack::lmnn(data.frame(xtrain4,ytrain))$output #train metric on the R4 sample
      a4[i]=mean(class::knn(xtrain4%*%t(M4),xtest4%*%t(M4),ytrain)!=ytest) #run 1-NN on the R4 set with metric
      b4[i]=mean(class::knn(xtrain4,xtest4,ytrain)!=ytest) #run 1-NN on the R4 set without metric
    }
  }
  if('C' %in% W){
    for(i in 1:t){
      R5=stats::rnorm(ncol(x)); R5=pracma::Toeplitz(R5[c(1,ncol(x):(ncol(x)-k+2))],R5)%*%diag(sample(c(-1,1),ncol(x),T))/sqrt(k) #circulant transform (Hinrichs and Vybíral, 2011) and (Freksen and Larsen, 2020)
      x5=x%*%t(R5) #compress original data set with R5
      xtrain5=x5[s,] #training set of R5
      xtest5=x5[-s,] #test set of R5
      M5=mlpack::lmnn(data.frame(xtrain5,ytrain))$output #train metric on the R5 sample
      a5[i]=mean(class::knn(xtrain5%*%t(M5),xtest5%*%t(M5),ytrain)!=ytest) #run 1-NN on the R5 set with metric
      b5[i]=mean(class::knn(xtrain5,xtest5,ytrain)!=ytest) #run 1-NN on the R5 set without metric
    }
  }
  cbind(d=ncol(x),n=nrow(x),k=k,
        am1=mean(a1),asd1=stats::sd(a1),bm1=mean(b1),bsd1=stats::sd(b1),
        am2=mean(a2),asd2=stats::sd(a2),bm2=mean(b2),bsd2=stats::sd(b2),
        am3=mean(a3),asd3=stats::sd(a3),bm3=mean(b3),bsd3=stats::sd(b3),
        am4=mean(a4),asd4=stats::sd(a4),bm4=mean(b4),bsd4=stats::sd(b4),
        am5=mean(a5),asd5=stats::sd(a5),bm5=mean(b5),bsd5=stats::sd(b5))
}
#'
#' Finds the k nearest neighbours of same labels (Weinberger and Saul (2009), Section 3.1).
#' @param x is the set whose rows are the instances sorted by labels
#' @param y is the set of labels in class 'factor'
#' @param k is how many nearest neighbours to consider
#' 
#' @export
KNN=function(x,y,k){
  a=NULL
  b=c(0,as.vector(table(y)))
  for(i in 1:(length(b)-1)){
    a=rbind(a,dbscan::kNN(x[y==levels(y)[i],],k)$id+b[i])
  }
  a
}
