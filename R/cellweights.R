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
covCW = function(x, wts = NULL){
  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(wts)){
    wts= apply(x,2, function(xvec){as.numeric( robustbase::covMcd(xvec)$mcd.wt )})
  }
  wts = wts%*%diag(sqrt(n/colSums(wts)))
  muhat = colSums(x*wts)/colSums(wts)
  xcenter = x - rep(1,n)%*%t(muhat)
  covmatrix = t(xcenter*wts)%*%(xcenter*wts)/n


  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))
  return(list(covmatrix = covmatrix, cormatrix = cormatrix, muhat = muhat))

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

  covmatrix = cellWise::DI(x)$cov

  covmatrix[group==1, group==1] = cov(x[,group==1])
  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))

  return(list(covmatrix = covmatrix, cormatrix = cormatrix))

}


covBC2 = function(x, group){
  n = dim(x)[1]
  p = dim(x)[2]

  covmatrix = cellWise::DI(x)$cov
  covmatrix[group!=3, group!=3] = robustbase::covMcd(x[,group!=3])$cov
  covmatrix[group==1, group==1] = cov(x[,group==1])
  scale = sqrt(diag(covmatrix))
  cormatrix = (diag(1/scale))%*%covmatrix%*%(diag(1/scale))
  return(list(covmatrix = covmatrix, cormatrix = cormatrix))

}


detectf = function(xvec, mu, Sigma, th){
  p = length(xvec)
  add1f = function(alter, basis){
    basisa = c(basis, alter)
    m = mahalanobis(xvec[basisa], mu[basisa], Sigma[basisa, basisa])
    return(m)
  }

  basis = NULL
  repeat{
    alters = setdiff(1:p, basis)
    mvec = sapply(alters, add1f, basis = basis)
    label = alters[which.min(mvec)]
    if(min(mvec) > qchisq(th,df = length(basis)+1)){break}
    #min(mvec) > qchisq(th,df = length(basis)+1)
    basis = sort(c(basis, label))
    if(length(basis)>=p){break}
  }
  wts = rep(0,p)
  wts[basis]=1
  return(wts)

}

DE = function(x, mu = NULL, Sigma = NULL, th = 0.95){
  if(is.null(Sigma)){
    fit = covCW(x)
    mu = fit$muhat
    Sigma = fit$covmatrix
  }

  k = 0
  repeat{
    wts = t(apply(x,1, detectf, th = th, mu = mu, Sigma = Sigma))
    fit = covCW(x, wts = wts)
    mu = fit$muhat
    if((KLdiv(Sigma, fit$covmatrix)<1e-6)|(k>20)){break}
    #KLdiv(Sigma, fit$covmatrix)
    Sigma = fit$covmatrix
    k = k+1
  }
  return(list(mu = mu, Sigma = Sigma, wts = wts))

}












