###using functions from file: Quadratic classifier.R
###parameters
d=1e+2 #ambient dimension
t=1e+1 #number of random matrices for each projection dimension
n=1e+3 #number of points
x=uniformly::runif_on_sphere(n,d,1) #instances in ambient space
A0=matrix(rnorm(d^2),d) #decision boundary matrix
b0=rep(0,d) #decision boundary vector
c0=0 #decision boundary constant
y=qucl(x,A0,b0,c0) #original labels according to quadratic decision rule
table(y) #table of labels
s=sample(1:nrow(x),.8*nrow(x))
xtrain=x[s,] #training set
xtest=x[-s,] #test set
ytrain=y[s] #training labels
ytest=y[-s] #test labels

#iterations
emperr=NULL
for(k in (1:10)*5){
  err=0
  for(i in 1:t){
    B=matrix(rnorm(d*k),k)/sqrt(k) #projection matrix
    Q=quadsgd(xtrain,ytrain,B) #train the quadratic classifier
    err=err+mean(qucl(xtest%*%t(B),Q[[1]],Q[[2]],Q[[3]])==ytest)/t #mean 0-1 empirical error
  }
  cat('Only',t-i,'iterations remain.\n')
  emperr=c(emperr,err) #mean empirical error
}
