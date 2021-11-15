require(robustbase)

## regression will cellwise weights
#' Title
#'
#' @param x
#' @param group
#'
#' @return
#' @export
#'
#' @examples
covCW = function(x){
  n = dim(x)[1]
  p = dim(x)[2]

  wts = matrix(1, nrow = n, ncol = p)
  wts[!is.na(match(as.vector(x),boxplot(as.vector(x),plot=FALSE)$out))] = 0
  wts = wts%*%diag(sqrt(n/colSums(wts)))
  muhat = colSums(x*wts)/colSums(wts)
  xcenter = x - rep(1,n)%*%t(muhat)
  covmatrix = t(xcenter*wts)%*%(xcenter*wts)/n


  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))
  return(list(covmatrix = covmatrix, cormatrix = cormatrix))

}


## regression will cellwise weights
#' Title
#'
#' @param x
#' @param group
#'
#' @return
#' @export
#'
#' @examples
covBC = function(x, group){
  n = dim(x)[1]
  p = dim(x)[2]

  covmatrix = robcovsel::covf(x, cor.method = "pair",
                              scale.method = "qn", pda.method = F)$covmatrix
  covmatrix[group!=3, group!=3] = robustbase::covMcd(x[,group!=3])$cov
  covmatrix[group==1, group==1] = cov(x[,group==1])
  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))

  return(list(covmatrix = covmatrix, cormatrix = cormatrix))

}








