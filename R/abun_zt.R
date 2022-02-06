##' Semiparametric empirical likelihood inference for abundance under zero-truncated models
##'
##' @description \code{abun_zt} is used to implement the maximum empirical likelihood method by fitting zero-truncated count regression models.
##'
##' @param y number of times that individuals were captured
##' @param K number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.
##' @param x vector or matrix containing the individual covariates which have influence on the capture probability.
##' @param dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @param maxN number specifying the largest searching value for \eqn{N} in EM algorithm.
##' @param start_beta vector specifying the starting value of the coefficient \eqn{\beta} in count regression model.
##' @param eps positive convergence tolerance \eqn{\epsilon}. The iterations converge in EM algorithm when the increase of the log-EL is less than \eqn{\epsilon}.
##' @param maxit integer specifying the maximal number of iterations in EM algorithm.
##' @param N0 number specifying the value of abundance. If it is \code{NULL}, the maximum EL estimator of abundance is calculated.
##'
##' @return An \code{abun_oi} object.
##'
##' @references
##'
##' Liu, Y., Li, P., Liu, Y., and Zhang, R. (2021).
##' Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data.
##' \emph{Biometrical Journal}.
##'
##' @importFrom methods new
##' @importFrom stats coef glm optimize nlminb plogis
##'
##' @export
##'
##'
abun_zt <- function ( y, K = NULL, x, dist, maxN = NULL,
                      start_beta = NULL, eps = 1e-5, maxit = 5000, N0 = NULL ) {

  ### preparation
  y <- as.numeric(y)
  n <- length(y)
  if( is.null(maxN) ) maxN <- 100*n
  if ( maxN < n ) stop("value of 'maxN' must be greater than the number of individuals captured")

  x_mat <- as.matrix( x )
  nx <- rbind(x_mat, x_mat)
  ny <- c(y, rep(0,n))
  nk <- rep(K, 2*n)

  ### Step 0
  if (is.null(start_beta)) start_beta <- rep(0, ncol(x_mat))
  beta <- as.matrix( start_beta )
  prob <- rep(1/n, n)

  if (dist == "binomial") {
    g <- plogis( as.numeric(x_mat%*%beta) )
    phi <- (1 - g)^K
  }

  if (dist == "poisson") {
    phi <- as.numeric( exp( - exp( x_mat%*%beta ) ) )
  }

  alpha <- sum(phi * prob)

  N <- ifelse ( is.null(N0), n/( 1 - alpha + 1e-300 ), N0 )

  pars<- c(N, alpha, beta)
  likes <- loglikelihood_zt(N, y, x_mat, beta, alpha, prob, dist, K)


  ### iteration

  err <- 1; nit <- 0

  while (err > eps & nit <= maxit) {

    nit <- nit + 1

    ### calculate wi
    wi <- ( N - n ) * phi * prob/( alpha + 1e-300 )
    nwi <- c(rep(1,n), wi)

    ### update beta
    if (dist == "binomial") {
      out <- glm(cbind(ny, nk-ny) ~ nx - 1, family="binomial", weights=nwi)
      beta <- as.matrix( coef(out) )
      g <- plogis( as.numeric(x_mat%*%beta) )
      phi <- (1 - g)^K
    }

    if (dist == "poisson") {
      out <- glm(ny ~ nx - 1, family="poisson", weights=nwi)
      beta <- as.matrix( coef(out) )
      phi <- as.numeric( exp( - exp( x_mat%*%beta ) ) )
    }

    ### update prob & alpha
    prob <- (wi + 1) / N
    alpha <- sum( phi * prob )


    ### update N
    if ( is.null(N0) ){

      obj_N <- function (nt) sum( log( nt + 1 - c(1:n) ) ) + (nt - n)*log(alpha + 1e-300)

      N <- optimize(obj_N, lower=n, upper=maxN, maximum=TRUE, tol = 0.01)$maximum
    }


    ### calculate the log-likelihood
    pars <- rbind(pars, c(N, alpha, beta))
    likes <- c(likes, loglikelihood_zt(N, y, x_mat, beta, alpha, prob, dist, K))


    ### stopping criterion
    err <- ( likes[nit+1] - likes[nit] )

  }

  AIC <- 2*( - likes[nit+1] + 2 + length(beta))

  rt <- new('abun_oi', model = "zt", dist = dist,
            N = N, beta = as.numeric(beta), alpha = alpha,
            loglikelihood = likes[nit+1], AIC = AIC,
            prob = prob, nit = nit, pars_trace = pars, loglikelihood_trace = likes,
            y = y, x = x_mat, start_beta =  start_beta,
            epsilon = eps, maxN = maxN, maxit = maxit)

  if (dist == "binomial") rt@K <- K
  return(rt)
}


##' EL log-likelihood under zero-truncated models
##'
##' @param N number specifying the value of abundance.
##' @param y number of times that individuals were captured
##' @param x matrix containing the individual covariates which have influence on the capture probability.
##' @param beta vector of regression cofficients in the capture probability model
##' @param alpha number, the probability of never being captured
##' @param prob vector, the probability mass function of covariates
##' @param dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @param K number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.
##'
##' @return The value of EL log-likelihood under zero-truncated models.
##'
##' @references
##'
##' Liu, Y., Li, P., Liu, Y., and Zhang, R. (2021).
##' Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data.
##' \emph{Biometrical Journal}.
##'
##' @importFrom stats plogis
##'
##' @export
##'
##'
loglikelihood_zt <- function (N, y, x, beta, alpha, prob, dist, K) {

  xb <- as.numeric( x%*%beta )
  n <- length(xb)

  if (dist == "binomial") {
    g <- plogis( xb )

    rt <- sum( log( N + 1 - c(1:n) ) ) - sum( log(1:n) ) + (N - n) * log( alpha + 1e-300 ) +
      sum( y * log(g + 1e-300) + (K - y) * log(1 - g + 1e-300) + log(choose(K, y)) ) +
      sum( log( prob + 1e-300 ) )
  }

  if (dist == "poisson") {
    rt <- sum( log( N + 1 - c(1:n) ) ) - sum( log(1:n) ) + (N - n) * log( alpha + 1e-300 ) +
      sum( y * xb - exp(xb) - lfactorial(y) ) + sum( log( prob + 1e-300 ) )
  }
  rt

}

