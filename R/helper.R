#' generating predictors and response randomly using default settings
#'
#' @param n sample size
#' @param p dimension
#' @param e contamination rate
#' @param r correlation coefficent
#' @param group magnitude of outliers
#'
#' @return x, the generated design matrix
#'
#' @export
#'
#' @examples
#'
#' dat = genevar()
#' x = dat$x


genevar = function(n = 100, p = 15, e = 0.2, r = 0.9, gamma = 10, mu, sigma, type){

  if(is.null(mu)){
    mu = rep(10,p)
    sigma = diag(rep(1^2,p))
    for (i in 1:p) {for (j in 1:p) {
      if (i !=j)sigma[i,j] = sqrt(sigma[i,i]*sigma[j,j])*r^abs(i-j)}}
  }

  if(type == "blocked"){

    bi1 = matrix(0, nrow = n, ncol = p-10)
    bi2 = matrix(sample(c(rep(1,e*n),rep(0,((1-e)*n)))), nrow = n, ncol = 5)
    bi3 = apply(matrix(0, nrow = n, ncol = 5), 2,
                function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
    bi = cbind(bi1, bi2, bi3)

    outl = rnorm(n = n*p, mean = gamma, sd = 1)
    rsign = sample(c(-1,1), size = n*p, replace = T)
    outlier = matrix(outl*rsign, nrow = n, ncol=p)
  }

  if(type == "structured"){
    bi = apply(matrix(FALSE, nrow = n, ncol = p), 2,
               function(xvec) {xvec[sample(x = 1:n, size = e*n)] = TRUE; return(xvec)})
    gene = function(labels){
      #labels = bi[1,]
      k = sum(labels)
      outl = rep(0,p)

      if(k!=0){
        muk = mu[labels]
        sigmak = sigma[labels, labels]
        u = eigen(sigmak)$vectors[,k]
        v = gamma*sqrt(k)*t(u)/mahalanobis(u+muk, muk, sigmak)
        outl[labels] = v
      }

      return(outl)
    }

    outlier = t(apply(bi, 1, gene))
    dim(outlier)

  }

  if(type == "cellwise"){
    bi = apply(matrix(0, nrow = n, ncol = p), 2,
                function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
    outl = rnorm(n = n*p, mean = gamma, sd = 1)
    rsign = sample(c(-1,1), size = n*p, replace = T)
    outlier = matrix(outl*rsign, nrow = n, ncol=p)
  }

  {
    xr = mvtnorm::rmvnorm(n = n, mean = mu, sigma = sigma)
    x = xr*(1-bi)+(outlier+mu)*bi
  }


  return(list(x = x, sigma = sigma, bi = bi))

}


KLdiv = function(x0, xhat){
  value = matrixcalc::matrix.trace(xhat%*%solve(x0)) - log(det(xhat%*%solve(x0))) - dim(x0)[1]
  return(value)
}










