#' Convert a table to a form suitable for \code{matrix plots} in \code{pgfplots}.
#' 
#' @param x matrix whose rows are the points; its first two columns must be categorical and contain all possible pairs of values
#' 
#' @import utils
#' @return a table for use in \code{matrix plot} in \code{pgfplots}
#'
matplot=function(x){
  cx=length(unique(x[,1]))
  cy=length(unique(x[,2]))
  a=x[order(x[,1]),]
  cbind(rep(1:cx,cy),rep(1:cy,each=cx),a[order(a[,2]),3])
}
#'
#' Complex covariance matrix.
#' 
#' @param x matrix whose rows are the points; real or complex
#'
#' @return the sample covariance matrix of \code{x}
#'
ccov=function(x){
  t(sweep(x,2,apply(x,2,mean)))%*%Conj(sweep(x,2,apply(x,2,mean)))/(nrow(x)-1)
}
#'
#' Pseudo-covariance matrix.
#'
#' @param x matrix whose rows are the points; real or complex
#' 
#' @return the sample pseudo-covariance matrix of \code{x}
#' 
pcov=function(x){
  t(sweep(x,2,apply(x,2,mean)))%*%sweep(x,2,apply(x,2,mean))/(nrow(x)-1)
}

#' Sample from multi-variate complex Gaussian distribution.
#' 
#' @param n number of vectors to sample
#' @param m mean vector of the distribution
#' @param C covariance matrix of the distribution
#' @param P pseudo-covariance matrix of the distribution
#' 
#' @import mvtnorm
#' @return matrix whose rows are the sampled vectors
rcmnorm=function(n,m,C,P){
  x=rmvnorm(n,c(Re(m),Im(m)),
            .5*rbind(cbind(Re(C)+Re(P),-Im(C)+Im(P)),cbind(Im(C)+Im(P),Re(C)-Re(P))))
  x[,1:ncol(C)]+1i*x[,(ncol(C)+1):(2*ncol(C))]
}

#' Checks if covariance and pseudo-covariance matrices are compatible.
#'
#' @inheritParams rcmnorm
#' 
#' @return two columns of eigenvalues; must be non-negative to be compatible
#' 
covcheck=function(C,P){
  cbind(covariance=eigen(C)$values,
        expression=eigen(Conj(C)-Conj(P)%*%solve(C)%*%P)$values)
}

#' Linear transformation into identity and diagonal covariances.
#' 
#' @param A Hermitian matrix to be turned into identity
#' @param B Hermitian matrix to be turned into diagonal
#' 
#' @return matrix that must be applied to turn \code{A} into identity and \code{B} into diagonal
#' 
diagtrans=function(A,B){
  U=eigen(A)
  D=diag(1/sqrt(U$values))
  U=U$vectors
  Conj(t(eigen(D%*%Conj(t(U))%*%B%*%U%*%D)$vectors))%*%D%*%Conj(t(U))
}

#' Permutation of vectorization.
#' 
#' @param m number of rows
#' @param n number of columns
#' 
#' @return permutation matrix \code{P} such that \code{P%*%c(A)=c(t(A))}
#' 
vecperm=function(m,n){
  a=c(matrix(1:(m*n),n,byrow=T))
  P=diag(0,m*n)
  for(i in 1:(m*n)){
    P[i,a[i]]=1
  }
  P
}

#' Sample uniformly \code{d} eigenvalues given largest eigenvalue and effective rank.
#' 
#' @param d number of dimensions
#' @param l largest eigenvalue
#' @param r desired effective rank
#' 
#' @import Surrogate
#' @return vector of eigenvalues with desired effective rank
#' @export
#' @examples
#' #generates all ones
#' 
#' reff(5,1,5) 
#' 
#' #generates randomly
#' 
#' reff(5,1,3)
#' 
#' #generates all-but-one zeros
#' 
#' reff(5,1,1)
reff=function(d,l=1,r=d){
  if(r>=d){
    rep(l,d)
  }else if(r<=1){
    c(l,rep(0,d-1))
  }else{
    c(l,rev(sort(RandVec(b=l,s=(r-1)*l,n=d-1)[[1]])))
  }
}

#' Create a \code{d*d} Walsh-Hadamard matrix.
#' 
#' @param d number of dimensions
#'
#' @return a \code{d*d} matrix with the Walsh-Hadamard property
#' @export
#' 
#' @examples
#' # 4*4 Walsh-Hadamard matrix
#' 
#' WHmat(4)
#' 
WHmat=function(d){
  A=matrix(NA,d,d)
  for(i in 1:d){
    for(j in 1:d){
      A[i,j]=(-1)^((as.numeric(intToBits(i-1))%*%as.numeric(intToBits(j-1)))%%2)/sqrt(d)
    }
  }
  A
}
