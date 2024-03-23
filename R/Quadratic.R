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
#' @examples 
#' # direct label prediction
#' 
#' sign(qucl(matrix(rnorm(500),5)))
qucl=function(x,A=diag(ncol(x)),b=rep(0,ncol(x)),c=0){
  c(apply(x*x%*%A,1,sum)+x%*%b+c)
}

#' Convex optimization for projection onto convex set using CVXR.
#' @param X matrix to be projected onto convex set
#' @param C inverse matrix that skews the convex set
#' @param a norm-constraint of the convex set
#' @param nor which norm to constrain
#' 
#' @return the matrix in the convex set closest in norm \code{nor} to \code{X}
#' @export

qproj=function(X,a,C=diag(nrow=nrow(X),ncol=ncol(X)),nor=c('nuc',1,2,'inf','fro')){
  A=CVXR::Variable(nrow(X),ncol(X),symmetric=T)
  obj=CVXR::cvxr_norm(X-A,'nuc') #any norm can be used
  constr=list(CVXR::cvxr_norm(C%*%A%*%C,nor[1])<=a) #any norm can be used
  prob=CVXR::Problem(CVXR::Minimize(obj),c(constr))
  result=CVXR::solve(prob)[[1]] #quadratic classifier matrix
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
  ((1-g<x)&(x<1))*(1-x)^2/(2*g)+(x<1-g)*(1-x-g/2)
}

#' Derivative of smoothed hinge loss
#' 
#' @param x vector
#' @param g parameter
#' 
#' @return vector of values
#' @export
smhingeder=function(x,g=.5){
  ((1-g<x)&(x<1))*(x-1)/g-(x<1-g)
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
#' @param a norm-constraint for the matrix in the quadratic boundary
#' @param nor which norm to constrain

quadsgd=function(x,y,B=diag(ncol(x)),g=.5,epoch=2,alpha=1,A=diag(0,nrow(B)),b=rep(0,nrow(B)),c=0,a=3,nor=c('nuc',1,2,'inf','fro')){
  C=svd(B)$u%*%diag(1/svd(B)$d)%*%t(svd(B)$u) #matrix that skews quadratic class; or C=diag(nrow(B))
  x=x%*%t(B) #compress set
  error=NULL
  for(i in 1:epoch){
    for(j in 1:nrow(x)){
      L=alpha*smhingeder(y[j]*qucl(x[j,],A,b,c),g)
      A=A-L*y[j]*x[j,]%*%t(x[j,]) #update matrix
      b=b-L*y[j]*x[j,] #update vector; remove line to train a homogeneous classifier
      c=c-L*y[j] #update constant; remove line to train a homogeneous classifier
      Norm=npmr::nuclear(C%*%A%*%C) #base::norm(C%*%A%*%C,'') for other norms; remove C for no skew
      if(Norm>a){ #if A is out of feasible set
        A=qproj(A,a,C,nor) #project A onto the feasible class and print the time; remove C for no skew
        cat('Epoch',i,'index',j,'had norm',Norm,'now has',npmr::nuclear(C%*%A%*%C),'\n')
      }
      error=c(error,sum(smhinge(y*qucl(x,A,b,c),g))/nrow(x)) #record error
    }
  }
  list(A=A,b=b,c=c,error=error) #quadratic classifier and sequence of errors
}
