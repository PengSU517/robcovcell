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
covCW = function(x, group = NULL){
  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(group)){
    grpup = 1:p
  }

  wts = matrix(0, ncol = p, nrow = n)
  for (i in unique(group)) {
    labels = (group==i)
    wts[,labels] = robustbase::covMcd(x[,labels])$mcd.wt
  }

  ww = sqrt(n/colSums(wts))
  wtsa = wts%*%diag(ww)
  xtilde = x*wtsa

  covmatrix = t(xtilde)%*%xtilde/n
  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))
  return(list(covmatrix = covmatrix, cormatrix = cormatrix))

}




