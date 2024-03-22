###functions

#quadratic classifier
#x is the matrix whose rows are the instances
#A is the matrix in the decision boundary
#b is the vector in the decision boundary
#c is the constant in the decision boundary
#returns the vector of labels for each instance either -1 or 1

qucl=function(x,A,b,c){
  c(sign(apply(x*x%*%A,1,sum)+x%*%b+c))
}

#convex optimization for projection onto convex set using CVXR
#X is the k*k matrix to be projected onto Q_B
#C is the inverse matrix from \Cref{th:compressed quadratic class equivalence}
#a is the norm-constraint

qproj=function(X,a,C=diag(ncol(X))){
  A=CVXR::Variable(k,k,symmetric=T)
  obj=CVXR::cvxr_norm(X-A,'nuc') #any norm can be used
  constr=list(CVXR::cvxr_norm(C%*%A%*%C,'nuc')<=a) #any norm can be used
  prob=CVXR::Problem(CVXR::Minimize(obj),c(constr))
  result=CVXR::solve(prob)[[1]] #quadratic classifier matrix
  result
}

#smoothed hinge loss
#g is the parameter in (0,1]

smhinge=function(x,g=.5){
  ((1-g<x)&(x<1))*(1-x)^2/(2*g)+(x<1-g)*(1-x-g/2)
}

#derivative of smoothed hinge loss
#g is the parameter in (0,1]

smhingeder=function(x,g=.5){
  ((1-g<x)&(x<1))*(x-1)/g-(x<1-g)
}

###train quadratic classifier under projection with SGD
#x is the matrix whose rows are the instances
#y is the vector of labels
#B is the projection matrix; use diag(d) for no compression
#g is the parameter in the smoothed hinge loss
#epoch is the number of epochs
#alpha is the learning rate; can be made to vary for each iteration
#A is the matrix with which to initialise SGD
#b is the vector with which to initialise SGD
#c is the constant with which to initialise SGD
#a is the norm-constraint for the matrix in the quadratic boundary

quadsgd=function(x,y,B,g=.5,epoch=2,alpha=1,A=diag(0,nrow(B)),b=rep(0,nrow(B)),c=0,a=3){
  C=svd(B)$u%*%diag(1/svd(B)$d)%*%t(svd(B)$u) #matrix that skews quadratic class; or C=diag(nrow(B))
  x=x%*%t(B) #compress set
  for(i in 1:epoch){
    for(j in 1:nrow(x)){
      L=alpha*smhingeder(y[j]*c(x[j,,drop=F]%*%A%*%x[j,]+b%*%x[j,]+c),g)
      A=A-L*y[j]*x[j,]%*%t(x[j,]) #update matrix
      b=b-L*y[j]*x[j,] #update vector; remove line to train a homogeneous classifier
      c=c-L*y[j] #update constant; remove line to train a homogeneous classifier
      nor=npmr::nuclear(C%*%A%*%C) #base::norm(C%*%A%*%C,'') other norms; remove C for no skew
      if(nor>a){ #if A is out of feasible set
        TIME=Sys.time() #note the time
        A=qproj(A,a,C) #project A onto the feasible class and print the time; remove C for no skew
        cat('Epoch',i,'index',j,'had norm',nor,'and cost',Sys.time()-TIME,'sec\n')
      }
    }
  }
  list(A=A,b=b,c=c) #output trained quadratic classifier
}
