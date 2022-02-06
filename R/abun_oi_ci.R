##' \code{abun_oi_ci} is used to calculate the empirical likelihood ratio confidence interval of abundance
##'
##' @title EL ratio confidence interval of abundance
##'
##' @param object a \code{abun_oi} object.
##' @param level number specifying the nominal level.
##' @param maxN_ci number specifying the largest searching value.
##' @return a vector, the empirical likelihood ratio confidence interval of abundance.
##' @importFrom stats qchisq uniroot
##' @export
##'
abun_oi_ci <- function (object, level = 0.95, maxN_ci = 1e9) {

  y <- object@y
  n <- length(y)
  if ( maxN_ci < n ) stop("value of 'maxN_ci' must be greater than the number of individuals captured")

  like_full <- object@loglikelihood
  hatN <- object@N

  ###########   ELR based CI  ##############
  rn <- function (N) {

    like_null <- switch (object@model,
                         "ztoi" = abun_ztoi ( y = y, K = object@K, x = object@x, z = object@z,
                                              dist = object@dist, maxN = object@maxN,
                                              start_beta = object@start_beta,
                                              start_eta = object@start_eta,
                                              eps = object@epsilon,
                                              maxit = object@maxit, N0 = N ),
                         "oizt" = abun_oizt ( y = y, K = object@K, x = object@x, z = object@z,
                                              dist = object@dist, maxN = object@maxN,
                                              start_beta = object@start_beta,
                                              start_eta = object@start_eta,
                                              eps = object@epsilon,
                                              maxit = object@maxit, N0 = N ),
                         "zt" = abun_zt ( y = y, K = object@K, x = object@x,
                                          dist = object@dist, maxN = object@maxN,
                                          start_beta = object@start_beta,
                                          eps = object@epsilon,
                                          maxit = object@maxit, N0 = N ) )

    2 * (like_full - like_null@loglikelihood) - qchisq(level, 1)

  }

  hatN <- object@N
  ntemp <- c(hatN, 5 * hatN)
  nit <- 2
  ind <- rn(ntemp[nit]) <= 0
  while ( ind ) {
    ntemp <- c(ntemp, ntemp[nit]*5)
    nit <- nit+1
    ind <- rn(ntemp[nit]) <= 0

    if (ntemp[nit] >= maxN_ci) {
      ci_upper <- maxN_ci
      break
    }
  }

  if (!ind & ntemp[nit] < maxN_ci) ci_upper <- uniroot( rn, c(ntemp[nit-1], ntemp[nit]), tol=0.01 )$root


  if ( rn(n) <= 0 ) {
    ci_lower <- n
  } else {
    ci_lower <- uniroot( rn, c(n, hatN), tol=0.01 )$root
  }

  c(ci_lower, ci_upper)

}


