##' Semiparametric empirical likelihood inference for abundance under zero-truncated one-inflated models
##'
##' @description \code{abun_ztoi} is used to implement the maximum empirical likelihood method by fitting zero-truncated one-inflated count regression models.
##'
##' @param y number of times that individuals were captured
##' @param K number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.
##' @param x vector or matrix containing the individual covariates which have influence on the capture probability.
##' @param z vector or matrix containing the individual covariates which have influence on the probability of one-inflation.
##' @param dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @param maxN number specifying the largest searching value for \eqn{N} in EM algorithm.
##' @param start_beta vector specifying the starting value of the coefficient \eqn{\beta} in count regression model.
##' @param start_eta vector specifying the starting value of the coefficient \eqn{\eta} in one-inflated regression model.
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
abun_ztoi <- function ( y, K = NULL, x, z = rep(1, length(y)),
                        dist, maxN = NULL, start_beta = NULL, start_eta = NULL,
                        eps = 1e-5, maxit = 5000, N0 = NULL ) {

  ### preparation
  f_fun <- function (y, xb, dist) {
    if ( dist == 'binomial' ) {
      prob <- plogis(xb)
      out <- choose(K, y) * prob^y * (1-prob)^{K-y}

    }
    if (dist == 'poisson') {
      lam <- exp(xb)
      out <- lam^{y} * exp(-lam) / factorial(y)
    }
    out
  }


  x_mat <- as.matrix(x)
  z_mat <- as.matrix(z)
  x_mat <- as.matrix( x_mat[order(y),] )
  z_mat <- as.matrix( z_mat[order(y),] )
  y <- sort(y)
  m <- sum(y==1)
  n <- length(y)
  if( is.null(maxN) ) maxN <- 100*n
  if ( maxN < n ) stop("value of 'maxN' must be greater than the number of individuals captured")


  ### Step 0
  if (is.null(start_beta)) start_beta <- rep(0, ncol(x_mat))
  if (is.null(start_eta)) start_eta <- rep(0, ncol(z_mat))

  ### w = logit(eta)
  beta <- as.matrix( start_beta )
  eta <- as.matrix( start_eta )
  prob <- rep(1/n, n)

  xb <- as.numeric( x_mat%*%beta )
  ze <- as.numeric( z_mat%*%eta )
  w <- plogis( ze )
  alpha <- sum( w * f_fun(rep(0, n), xb, dist) * prob )
  N <- ifelse ( is.null(N0), n/( 1 - alpha + 1e-20 ), N0 )

  pars<- c(N, alpha, beta, eta)
  likes <- loglikelihood_ztoi (N, y, x_mat, beta, z_mat, eta, alpha, prob, dist, K)

  ### iteration
  err <- 1; nit <- 0

  while ( err > eps & nit <= maxit ) {

    nit <- nit + 1

    ### calculate omega
    f <- f_fun( rep(1,m), xb[1:m], dist )
    v <- c( w[1:m]*f/(1-w[1:m]+w[1:m]*f), rep(1, n-m) ) # nu
    u <- (N-n) * w * f_fun( rep(0,n), xb, dist ) * prob/( alpha +1e-300 ) # omega

    ### update w
    if ( all(z == 1) ) {
      w <- rep(sum(v + u)/sum(1 + u), n)
      eta <- as.matrix( log(w[1]/(1 - w[1] + 1e-300)) )
    } else {

      nyz <- c(rep(1,n), rep(0,n))
      nz <- rbind(z_mat, z_mat)
      outw <- glm(cbind(nyz, 1-nyz) ~ nz - 1, family = "binomial", weights = c(v+u, 1 - v))
      eta <- as.matrix( coef(outw) )
      w <- as.numeric( plogis(z_mat%*%eta) )
    }

    ### update beta
    ny <- c(y, rep(0,n))
    nx <- rbind(x_mat, x_mat)

    if (dist == 'binomial')
      out <- glm(cbind(ny, K-ny) ~ nx - 1, family=dist, weights=c(v, u))

    if (dist == 'poisson')
      out <- glm(ny ~ nx - 1, family=dist, weights=c(v, u))

    beta <- as.matrix( coef(out) )
    xb <- as.numeric( x_mat%*%beta )


    ### update prob & alpha
    prob <- (u + 1) / N
    alpha <- sum( w * f_fun(rep(0,n), xb, dist) * prob )


    ### update N
    if ( is.null(N0) ) {
      obj <- function (nt) sum(log((nt - n + 1):nt)) + (nt - n)*log(alpha + 1e-20)
      out <- optimize(obj, lower=n, upper=maxN, maximum=TRUE)
      N <- out$maximum
    }

    ### calculate the log-likelihood
    pars <- rbind(pars, c(N, alpha, beta, eta))
    likes <- c(likes, loglikelihood_ztoi (N, y, x_mat, beta, z_mat, eta, alpha, prob, dist, K))

    ### stopping criterion
    err <- ( likes[nit+1] - likes[nit] )

  }

  AIC <- 2*( - likes[nit+1] + 2 + length(beta) + length(eta))

  rt <- new('abun_oi', model = "ztoi", dist = dist,
            N = N, beta = as.numeric(beta),
            eta = as.numeric(eta), alpha = alpha,
            loglikelihood = likes[nit+1], AIC = AIC,
            prob = prob, nit = nit, pars_trace = pars, loglikelihood_trace = likes,
            y = y, x = x_mat, z = z_mat,
            start_beta =  start_beta,
            start_eta = start_eta,
            epsilon = eps, maxN = maxN, maxit = maxit)

  if (dist == "binomial") rt@K <- K
  return(rt)

}


##' EL log-likelihood under zero-truncated one-inflated models
##'
##' @param N number specifying the value of abundance.
##' @param y number of times that individuals were captured
##' @param x matrix containing the individual covariates which have influence on the capture probability.
##' @param beta vector of regression cofficients in the capture probability model
##' @param z matrix containing the individual covariates which have influence on the probability of one-inflation.
##' @param eta vector of regression cofficients in the one-inflated probability model
##' @param alpha number, the probability of never being captured
##' @param prob vector, the probability mass function of covariates
##' @param dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @param K number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.
##'
##' @return The value of EL log-likelihood under zero-truncated one-inflated models.
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
loglikelihood_ztoi <- function (N, y, x, beta, z, eta, alpha, prob, dist, K) {

  xb <- as.numeric( x%*%beta )
  w <- as.numeric( plogis( z%*%eta ) )
  n <- length(xb)

  if (dist == 'binomial') {
    p <- plogis(xb)
    f <-  choose(K, y) * p^y * (1-p)^{K-y}
  }

  if (dist == 'poisson') {
    lam <- exp(xb)
    f <- lam^y * exp(-lam) / factorial(y)
  }

  h <- w * f + (1 - w) * (y==1)

  sum( log( N + 1 - c(1:n) ) ) - sum( log(1:n) ) +
    (N - n)*log( alpha + 1e-300 ) + sum( log( h + 1e-300 ) ) +
    sum( log( prob + 1e-300 ) )

}

