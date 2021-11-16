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


genevar = function(n = 100, p = 15, e = 0.05, r = 0.5, gamma = 10, group = NULL){

  if(is.null(group)){
    size = p/3
    group = c(rep(1,size), rep(2,size), rep(3,size))
  }

  {
    mu = rep(10,p)
    sigma = diag(rep(1^2,p))
    for (i in 1:p) {for (j in 1:p) {
      if (i !=j)sigma[i,j] = sqrt(sigma[i,i]*sigma[j,j])*r^abs(i-j)}}
  }

  {
    xr = mvtnorm::rmvnorm(n = n, mean = mu, sigma = sigma)
    bi1 = matrix(0, nrow = n, ncol = size)
    bi2 = matrix(sample(c(rep(1,e*n),rep(0,((1-e)*n)))), nrow = n, ncol = size)
    #bi3 = matrix(sample(c(rep(1,e*n),rep(0,((1-e)*n)))), nrow = n, ncol = size)
    bi3 = apply(matrix(0, nrow = n, ncol = size), 2,
                function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
    bi = cbind(bi1, bi2, bi3)
    outl = rnorm(n = n*p, mean = gamma, sd = 1)
    rsign = sample(c(-1,1), size = n*p, replace = T)
    outlier = matrix(outl*rsign, nrow = n, ncol=p)
    x = xr*(1-bi)+(outlier+10)*bi
  }
  return(list(x = x, group = group, sigma = sigma, bi = bi))

}
