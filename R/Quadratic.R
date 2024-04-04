###functions

#' Value of quadratic classifier for all instances.
#' @param x matrix whose rows are the instances
#' @param A matrix in the quadratic form
#' @param b vector in the quadratic form
#' @param c constant in the quadratic form
#' 
#' @return vector of values
#' @export
#' 
#' @importFrom Rdpack reprompt
#' 
#' @examples
#' # direct label prediction
#' 
#' sign(qucl(matrix(rnorm(500),5)))
qucl=function(x,A=diag(ncol(x)),b=rep(0,ncol(x)),c=0){
  c(apply(x*x%*%A,1,sum)+x%*%b+c)
}

#' CVXR for projection of matrix onto convex set.
#' @param X matrix to be projected onto convex set
#' @param a norm-constraint of the convex set
#' @param nor which matrix-norm to constrain; defaults to nuclear
#' @param C inverse of the matrix that skews the convex set; defaults to no skew
#' 
#' @return the matrix in the convex set closest in norm \code{nor} to \code{X}
#' @import CVXR
#' @export
qproj=function(X,a=1,nor=c('nuc',1,2,'inf','fro'),C=diag(ncol(X))){
  A=Variable(nrow(X),ncol(X),symmetric=T)
  obj=cvxr_norm(X-A,nor[1])
  constr=list(cvxr_norm(C%*%A%*%C,nor[1])<=a) #any norm can be used
  prob=Problem(Minimize(obj),c(constr))
  result=solve(prob)[[1]] #matrix inside the convex set
  result
}

#' CVXR for projection of vector onto convex set.
#' @param X vector to be projected onto convex set
#' @param a norm-constraint of the convex set
#' @param nor which p-norm to constrain; defaults to infinity
#' @param C inverse of the vector that skews the convex set; defaults to no skew
#' 
#' @return the vector in the convex set closest in norm \code{nor} to \code{X}
#' @export
vproj=function(X,a=1,nor=Inf,C=rep(1,length(X))){
  A=Variable(length(X))
  obj=p_norm(X-A,nor)
  constr=list(p_norm(C*A,nor)<=a) #any norm can be used
  prob=Problem(Minimize(obj),c(constr))
  result=solve(prob)[[1]] #matrix inside the convex set
  result
}
#' Smoothed hinge loss
#'
#' @param x input vector
#' @param g parameter
#' 
#' @return vector of values
#' @export
smhinge=function(x,g=.5){
  (((1-g)<x)&(x<1))*(1-x)^2/(2*g)+(x<=(1-g))*(1-x-g/2)
}

#' Derivative of smoothed hinge loss
#' 
#' @param x vector
#' @param g parameter
#' 
#' @return vector of values
#' @export
smhingeder=function(x,g=.5){
  (((1-g)<x)&(x<1))*(x-1)/g+(x<=(1-g))*(-1)
}
#' Train quadratic classifier under projection with SGD
#' @param x matrix whose rows are the instances
#' @param y vector of labels
#' @param B projection matrix; defaults to no projection
#' @param g parameter of the smoothed hinge loss
#' @param epoch number of epochs in the SGD
#' @param alpha learning rate of the SGD; can be modified to vary for each iteration
#' @param A matrix with which to initialise SGD
#' @param b vector with which to initialise SGD
#' @param c constant with which to initialise SGD
#' @param skew whether to skew the compressed class or not
#' @inheritParams qproj
#' 
#' @return list containing the following elements
#' 
#' \strong{`A`} matrix in the trained quadratic classifier
#' 
#' \strong{`b`} vector in the trained quadratic classifier
#' 
#' \strong{`c`} constant in the trained quadratic classifier
#' 
#' @export
quadsgd=function(x,y,B,g=.5,epoch=2,alpha=.1,A=diag(nrow(B)),b=rep(0,nrow(B)),c=0,a=1,nor=c('nuc',1,2,'inf','fro'),skew=T){
  if(skew==T){
    C=svd(B)$u%*%diag(svd(B)$d)%*%t(svd(B)$u) #matrix that skews quadratic class
  } else {
    C=diag(nrow(B)) #no skew in the quadratic class
  }
  x=x%*%t(B) #compress set
  for(i in 1:epoch){
    for(j in 1:nrow(x)){
      L=alpha*smhingeder(y[j]*qucl(x[j,],A,b,c),g)
      A=qproj(A-L*y[j]*x[j,]%*%t(x[j,]),a,nor,C) #update matrix and project onto quadratic class
      b=b-L*y[j]*x[j,] #update vector; remove line to train a homogeneous classifier
      c=c-L*y[j] #update constant; remove line to train a homogeneous classifier
    }
  }
  list(A=A,b=b,c=c) #output trained quadratic classifier
}
