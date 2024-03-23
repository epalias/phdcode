#' sample uniformly \code{d} eigenvalues given largest eigenvalue and effective rank
#' 
#' @param d The number of dimensions.
#' @param l The largest eigenvalue.
#' @param r The desired effective rank.
#' 
#' @return Vector of eigenvalues with desired effective rank.
#' @export

reff=function(d,l,r){
  if(r>=d){
    rep(l,d)
  }else if(r<=1){
    c(l,rep(0,d-1))
  }else{
    c(l,rev(sort(Surrogate::RandVec(b=l,s=(r-1)*l,n=d-1)[[1]])))
  }
}
