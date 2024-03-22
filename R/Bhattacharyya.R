###functions

#sample uniformly in d dimensions given largest eigenvalue and effective rank
#d is the number of dimensions
#l is the largest eigenvalue
#r is the effective rank

reff=function(d,l,r){
  if(r>=d){
    rep(l,d)
  }else if(r<=1){
    c(l,rep(0,d-1))
  }else{
    c(l,rev(sort(Surrogate::RandVec(b=l,s=(r-1)*l,n=d-1)[[1]])))
  }
}

###CODE FOR THE QDA-FINITE SAMPLE PAPER STARTS HERE###

#m0 is the true mean of Gaussian 0
#S0 is the true covariance of Gaussian 0
#m1 is the true mean of Gaussian 1
#S1 is the true covariance of Gaussian 1
#m0hat is the estimation of m0
#S0hat is the estimation of S0
#m1hat is the estimation of m1
#S1hat is the estimation of S1
#p0 is the prior of Gaussian 0
#p1 is the prior of Gaussian 1
#l is the parameter of the Chernoff bound
#n is the training sample size to estimate means and covariances
#N is the test sample size from Gaussian 0
#t is the number of iterations


###Bhattacharyya bound for finite-sample size

bhatbound=function(m0,S0,m0hat,S0hat,m1hat,S1hat,p0=.5,p1=.5,l=.5){
  S=solve(solve(S0)-l*solve(S0hat)+l*solve(S1hat))
  m=S%*%(-solve(S0)%*%m0+l*solve(S0hat)%*%m0hat-l*solve(S1hat)%*%m1hat)
  c((p1/p0)^l*exp(-.5*log((det(S0)*det(S1hat)^l)/(det(S)*det(S0hat)^l))-.5*(t(m0)%*%solve(S0)%*%m0-l*t(m0hat)%*%solve(S0hat)%*%m0hat+l*t(m1hat)%*%solve(S1hat)%*%m1hat-t(m)%*%solve(S)%*%m)))
}

###QDA error rate when sampling from only one Gaussian

qdaresult=function(m0,S0,m1,S1,p0=.5,p1=.5,l=.5,n=10,N=100,t=1000){
  a=rep(0,t)
  b=rep(0,t)
  c=rep(0,t)
  rep(bhatbound(m0,S0,m0,S0,m1,S1,p0,p1,l),t)
  for(i in 1:t){
    x0=mvtnorm::rmvnorm(n,m0,S0) #training sample
    x1=mvtnorm::rmvnorm(n,m1,S1) #training sample
    m0hat=apply(x0,2,mean)
    m1hat=apply(x1,2,mean)
    S0hat=cov(x0)
    S1hat=cov(x1)
    a[i]=bhatbound(m0,S0,m0hat,S0hat,m1hat,S1hat,p0,p1,l)
    b[i]=bhatbound(m0hat,S0hat,m0hat,S0hat,m1hat,S1hat,p0,p1,l)
    x0test=mvtnorm::rmvnorm(N,m0,S0) #test sample
    for(j in 1:N){
      c[i]=c[i]+(mvtnorm::dmvnorm(x0test[j,],m0hat,S0hat)<mvtnorm::dmvnorm(x0test[j,],m1hat,S1hat))/N
    }
  }
  cbind(B=rep(bhatbound(m0,S0,m0,S0,m1,S1,l),t),Bt=a,Be=b,P=c)
}

#use the following to run everything after setting m0,S0,m1,S1,n

apply(qdaresult(m0,S0,m1,S1,n),2,mean,na.rm=T)

###CODE FOR THE QDA-FINITE SAMPLE PAPER ENDS HERE###

###no projection

noprojexp=function(d,m,l0,l1,s0,s1,n){
  m0=rep(0,d)
  m1=c(uniformly::runif_on_sphere(1,d,m)) #mean distance of m
  A=pracma::randortho(d) #random orthogonal matrix of dimension d
  S0=diag(reff(d,l0,s0)) #largest eigenvalue is l0
  S1=A%*%diag(reff(d,l1,s1))%*%t(A) #largest eigenvalue is l1
  qdaerror(m0,S0,m1,S1,n) #empirical error
}


###random projection

#high-probability Bhattacharyya bound under random projection
#probability is at least 1-e

rpbb=function(m0,S0,m1,S1,k,e){
  S=(S0+S1)/2 #mean covariance matrix
  l0=eigen(S0)$values[1] #largest eigenvalue of S0
  l1=eigen(S1)$values[1] #largest eigenvalue of S1
  l=eigen(S)$values[1] #largest eigenvalue of S
  t=sum(diag(S)) #trace of S
  t0=sum(diag(S0)) #trace of S0
  t1=sum(diag(S1)) #trace of S1
  exp(-(1-e)*k*sum((m0-m1)^2)/(8*(sqrt(t)+sqrt(k*l)+e)^2))*
    (sqrt((sqrt(t0)+sqrt(k*l0)+e)*(sqrt(t1)+sqrt(k*l1)+e))/
       max(sqrt(t)-sqrt(k*l)-e,0))^k #bound
}

#experiments with random projection

rpexp=function(d,m,l0,l1,K,S,n,r){
  alls=NULL
  allk=NULL
  memperr=NULL #mean of empirical error under all random matrices
  semperr=NULL #standard deviation of empirical error under all random matrices
  bound=NULL
  m0=rep(0,d)
  m1=c(uniformly::runif_on_sphere(1,d,m)) #mean distance of m
  for(s in S){ #rank or effective rank of S0 or S1
    A=pracma::randortho(d) #random orthogonal matrix of dimension d
    S0=diag(rep(c(l0,0),c(s,d-s))) #largest eigenvalue is l0
    S1=A%*%diag(rep(c(l1,0),c(s,d-s)))%*%t(A) #largest eigenvalue is l1
    test=mvtnorm::rmvnorm(n,m0,S0) #sample a test set in the ambient space
    for(k in K){ #projection dimension; no higher than rank; for effective rank it can equal s
      emperr=NULL #empirical error for each Gaussian random matrix
      for(j in 1:r){ #for each Gaussian random matrix
        R=matrix(rnorm(k*d),k) #sample a Gaussian random matrix
        x=test%*%t(R) #project the test set
        m0k=R%*%m0 #project m0
        S0k=R%*%S0%*%t(R) #project S0
        m1k=R%*%m1 #project m1
        S1k=R%*%S1%*%t(R) #project S1
        a=0 #empirical error for current Gaussian random matrix
        for(i in 1:n){
          a=a+(mvtnorm::dmvnorm(x[i,],m0k,S0k)<mvtnorm::dmvnorm(x[i,],m1k,S1k)) #QDA decision rule under rp
        }
        emperr=c(emperr,a/n) #empirical error on the projected test set
      }
      memperr=c(memperr,mean(emperr)) #average empirical error
      semperr=c(semperr,sd(emperr)) #5% quantile of the empirical error
      bound=c(bound,rpbb(m0,S0,m1,S1,k,2.862)) #e so that probability is at least 0.95
      alls=c(alls,s) #append rank of S0 or S1
      allk=c(allk,k) #append projection dimension
      print(c(k,s)) #at which iteration we are (optional)
    }
  }
  cbind(s=alls,k=allk,memperr=memperr,semperr=semperr,bound=bound)
}

###random orthogonal projection

#Bhattacharyya bound under random orthogonal projection
#probability is at least 1-e

ropbb=function(m0,S0,m1,S1,k,e){
  S=(S0+S1)/2 #mean covariance matrix
  l0=eigen(S0)$values #largest eigenvalue of S0
  l1=eigen(S1)$values #largest eigenvalue of S1
  l=eigen(S)$values #largest eigenvalue of S
  d=ncol(S) #ambient dimension
  exp(-(1-e)*k*sum((m0-m1)^2)/(8*d*l[1]))*
    prod(sqrt(sqrt(l0[1:k]*l1[1:k])/l[d-k+(1:k)])) #bound
}

#experiments with random orthogonal projection

ropexp=function(d,m,l0,l1,K,S,n,r){
  alls=NULL
  allk=NULL
  memperr=NULL #mean of empirical error under all random matrices
  semperr=NULL #standard deviation of empirical error under all random matrices
  bound=NULL
  m0=rep(0,d)
  m1=c(uniformly::runif_on_sphere(1,d,m)) #mean distance of m
  for(s in S){ #rank or effective rank of S0 or S1
    A=pracma::randortho(d) #random orthogonal matrix of dimension d
    S0=diag(reff(d,l0,r)) #largest eigenvalue is l0
    S1=A%*%diag(reff(d,l1,r))%*%t(A) #largest eigenvalue is l1
    test=mvtnorm::rmvnorm(n,m0,S0) #sample a test set in the ambient space
    for(k in K){ #projection dimension; no higher than rank; for effective rank it can equal s
      emperr=NULL #empirical error for each Gaussian random matrix
      for(j in 1:r){ #for each Gaussian random matrix
        R=matrix(rnorm(k*d),k) #sample a Gaussian random matrix
        R=solve(expm::sqrtm(R%*%t(R)))%*%R #orthonormalize R
        x=test%*%t(R) #project the test set
        m0k=R%*%m0 #project m0
        S0k=R%*%S0%*%t(R) #project S0
        m1k=R%*%m1 #project m1
        S1k=R%*%S1%*%t(R) #project S1
        a=0 #empirical error for current Gaussian random matrix
        for(i in 1:n){
          a=a+(mvtnorm::dmvnorm(x[i,],m0k,S0k)<mvtnorm::dmvnorm(x[i,],m1k,S1k)) #QDA decision rule under rp
        }
        emperr=c(emperr,a/n) #empirical error on the projected test set
      }
      memperr=c(memperr,mean(emperr)) #average empirical error
      semperr=c(semperr,sd(emperr)) #5% quantile of the empirical error
      bound=c(bound,ropbb(m0,S0,m1,S1,k,2.862)) #e so that probability is at least 0.95
      alls=c(alls,s) #append rank of S0 or S1
      allk=c(allk,k) #append projection dimension
      print(c(k,s)) #at which iteration we are (optional)
    }
  }
  cbind(s=alls,k=allk,memperr=memperr,semperr=semperr,bound=bound)
}

###PCA

#Bhattacharyya bound under PCA

pcabb=function(m0,S0,m1,S1,k){
  S=(S0+S1)/2 #mean covariance matrix
  E=eigen(S) #eigenvalue decomposition of S
  U=E$vectors[,1:k] #k principal eigenvectors of S
  l=E$values #eigenvalues of S
  L=diag(1/sqrt(l[1:k])) #diagonal matrix with the inverse square roots of the k largest eigenvalues of S
  l0=eigen(S0)$values #eigenvalues of S0
  l1=eigen(S1)$values #eigenvalues of S1
  exp(-.125*sum((L%*%t(U)%*%(m0-m1))^2))*
    prod(sqrt(sqrt(l0[1:k]*l1[1:k])/l[1:k])) #bound
}

#experiments with PCA

pcaexp=function(d,m,l0,l1,K,S,n){
  allk=NULL
  alls=NULL
  emperr=NULL
  bound=NULL
  m0=rep(0,d)
  m1=c(uniformly::runif_on_sphere(1,d,m))
  for(s in S){ #rank or effective rank of S0 or S1
    A=pracma::randortho(d) #random orthogonal matrix of dimension d
    S0=diag(reff(d,l0,s)) #largest eigenvalue is l0
    S1=A%*%diag(reff(d,l1,s))%*%t(A) #largest eigenvalue is l1
    test=mvtnorm::rmvnorm(n,m0,S0) #sample a test set in the ambient space
    E=eigen((S0+S1)/2)$vectors
    for(k in K){
      Uk=E[,1:k] #matrix for PCA projection
      x=test%*%Uk #project the test set
      m0k=t(Uk)%*%m0 #PCA-projected m0
      S0k=t(Uk)%*%S0%*%Uk #PCA-projected S0
      m1k=t(Uk)%*%m1 #PCA-projected m1
      S1k=t(Uk)%*%S1%*%Uk #PCA-projected S1
      a=0 #misclassification of true Gaussian
      for(i in 1:n){
        a=a+(mvtnorm::dmvnorm(x[i,],m0k,S0k)<mvtnorm::dmvnorm(x[i,],m1k,S1k)) #QDA decision rule under pca
      }
      emperr=c(emperr,a/n) #empirical error of true Gaussian sample
      bound=c(bound,pcabb(m0,S0,m1,S1,k))
      allk=c(allk,k)
      alls=c(alls,s)
      print(c(s,k)) #at which iteration we are (optional)
    }
  }
  cbind(s=alls,k=allk,emperr=emperr,bound=bound) #can use confusion matrix for emperr
}
