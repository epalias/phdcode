#' Finite-sample Bhattacharyya bound between two Gaussians \code{0} and \code{1}.
#'
#' @param m0 true mean vector of Gaussian \code{0}
#' @param S0 true covariance matrix of Gaussian \code{0}
#' @param m0hat estimation of \code{m0}
#' @param S0hat estimation of \code{S0}
#' @param m1hat estimation of \code{m1}
#' @param S1hat estimation of \code{S1}
#' @param p0 prior of Gaussian \code{0}
#' @param p1 prior of Gaussian \code{1}
#' @param l parameter of the Chernoff bound
#' 
#' @return Chernoff distance between the two Gaussians
#' @export
bhatbound=function(m0,S0,m0hat,S0hat,m1hat,S1hat,p0=.5,p1=.5,l=.5){
  S=solve(solve(S0)-l*solve(S0hat)+l*solve(S1hat))
  m=S%*%(-solve(S0)%*%m0+l*solve(S0hat)%*%m0hat-l*solve(S1hat)%*%m1hat)
  c((p1/p0)^l*exp(-.5*log((det(S0)*det(S1hat)^l)/(det(S)*det(S0hat)^l))-.5*(t(m0)%*%solve(S0)%*%m0-l*t(m0hat)%*%solve(S0hat)%*%m0hat+l*t(m1hat)%*%solve(S1hat)%*%m1hat-t(m)%*%solve(S)%*%m)))
}

#' QDA error rate when sampling only from one Gaussian.
#'
#' @param n training sample size to estimate mean vectors and covariance matrices
#' @param N test sample size from Gaussian \code{0}
#' @param t number of iterations
#' @param m1 true mean vector of Gaussian \code{1}
#' @param S1 true covariance matrix of Gaussian \code{1}
#' @inheritParams bhatbound
#' 
#' @return Average classification error
#' @export
qdaresult=function(m0,S0,m1,S1,p0=.5,p1=.5,l=.5,n=10,N=100,t=1000){
  a=rep(0,t)
  b=rep(0,t)
  c=rep(0,t)
  rep(bhatbound(m0,S0,m0,S0,m1,S1,p0,p1,l),t)
  for(i in 1:t){
    x0=rmvnorm(n,m0,S0) #training sample
    x1=rmvnorm(n,m1,S1) #training sample
    m0hat=apply(x0,2,mean)
    m1hat=apply(x1,2,mean)
    S0hat=stats::cov(x0)
    S1hat=stats::cov(x1)
    a[i]=bhatbound(m0,S0,m0hat,S0hat,m1hat,S1hat,p0,p1,l)
    b[i]=bhatbound(m0hat,S0hat,m0hat,S0hat,m1hat,S1hat,p0,p1,l)
    x0test=rmvnorm(N,m0,S0) #test sample
    for(j in 1:N){
      c[i]=c[i]+(dmvnorm(x0test[j,],m0hat,S0hat)<dmvnorm(x0test[j,],m1hat,S1hat))/N
    }
  }
  cbind(B=rep(bhatbound(m0,S0,m0,S0,m1,S1,l),t),Bt=a,Be=b,P=c)
}
