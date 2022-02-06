##' Score test for one-inflation
##'
##' @description Function to calculate the score test statistics for the one-inflation based on zero-truncated one-inflated model and one-inflated zero-truncated model.
##'
##' @param object an abun_oi object produced by \code{\link{abun_oi}} with \code{model} option being \code{"zt"}.
##' @param model character specification of zero-truncated one-inflated model (\code{"ztoi"}, default) or one-inflated zero-truncated model (\code{"oizt"}).
##'
##' @return A number, the score test log-likelihood.
##'
##'
##' @export
##'
scoretest_oi <- function (object, model = c("ztoi", "oizt")) {

  model <- match.arg(model)
  switch(model,
         "ztoi" = scoretest_ztoi(object),
         "oizt" = scoretest_oizt(object))
}


scoretest_ztoi <- function (object) {

  dist = object@dist

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

  ef_fun <- function (xb, dist) {
    if (dist == 'binomial') {
      prob <- plogis(xb)
      out <- K * prob
    }
    if (dist == 'poisson') {
      lam <- exp(xb)
      out <- lam
    }
    out
  }

  varf_fun <- function (xb, dist) {
    if (dist == 'binomial') {
      prob <- plogis(xb)
      out <- K * prob * (1 - prob)
    }
    if (dist == 'poisson') {
      lam <- exp(xb)
      out <- lam
    }
    out
  }

  x_mat = object@x
  y <- object@y
  x_mat <- as.matrix( x_mat[order(y),] )
  y <- sort(y)
  m <- sum(y==1)
  n <- length(y)
  if (dist == "binomial") K <- object@K
  N <- object@N
  beta <- as.matrix(object@beta)
  xb <- as.numeric(x_mat%*%beta)

  alp <- object@alpha
  f0 <- f_fun(rep(0,n), xb, dist)
  f1 <- f_fun(rep(1,n), xb, dist)


  ####  inverse probability weighting
  varphi <- sum((1 - f0)^{-2})/N
  ef <- ef_fun(xb, dist)
  varf <- varf_fun(xb, dist)

  ### Vij's
  V11 <- 1 - 1/alp
  V13 <- 1/alp
  V22 <- t(x_mat) %*% diag(
    (f0 * ef^2/(1 - f0) - varf)/(1 - f0), n, n
  ) %*% x_mat/N

  V23 <- t(x_mat)%*%(f0*ef/(1-f0)^2)/N
  V24 <- (1 - alp)^2 * V23
  V33 <- varphi - 1/alp
  V34 <- (1 - alp)^2 * varphi
  V44 <- (1 - alp)^4 * varphi - (1 - alp)^3

  W <- rbind( cbind(-V11, matrix(0, 1, length(beta)),- V13),
              cbind(matrix(0, length(beta), 1),
                    -V22+V24%*%t(V24)/V44,
                    -V23+V24%*%t(V34)/V44),
              cbind(-V13,
                    -t(V23)+V34%*%t(V24)/V44,
                    -t(V33)+V34%*%t(V34)/V44 ))

  if (all(eigen(W)$values>0)) {
    W_chol <- forwardsolve( t(chol(W)), diag(1, nrow(W)))
    W_inv <- t(W_chol)%*%W_chol
    V_beta <- t(x_mat)%*%( (1 - ef)/ (1-f0))/N

    V_sid <- rbind(2, V_beta, 0)
    V_sid1 <- rbind(1, V_beta, 0)

    V25 <- t(x_mat)%*% (ef/(1 - f0)^2 - 1/(1-f0))/N
    V55 <- sum(1/(f1*(1-f0)))/N - 1
    V_sid2 <- rbind(-1, V25 - V24/(1-alp)^2, 0)
    varr <- t(V_sid1)%*%W_inv%*%V_sid2 + V55

    xbm <- x_mat[1:m,]%*%beta
    f11 <- f_fun(rep(1, m), xbm, dist)
    U <- N - sum(1/f11)

    rt <- list(W_inv = W_inv, V_sid = V_sid, U = U, varr = varr)
  } else{
    xbm <- x_mat[1:m,]%*%beta
    f11 <- f_fun(rep(1, m), xbm, dist)
    U <- N - sum(1/f11)
    rt <- list(W_inv = NA, V_sid = NA, U = NA, varr = -1e9)
  }

  S <- as.numeric( rt$U/sqrt(rt$varr * N) )
  names(S) <- "scoretest_ztoi"
  return(S)

}


scoretest_oizt <- function (object) {

  # object = out
  dist = object@dist

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


  ef_fun <- function (xb, dist) {
    if (dist == 'binomial') {
      prob <- plogis(xb)
      out <- K * prob
    }
    if (dist == 'poisson') {
      lam <- exp(xb)
      out <- lam
    }
    out
  }
  varf_fun <- function (xb, dist) {
    if (dist == 'binomial') {
      prob <- plogis(xb)
      out <- K * prob * (1 - prob)
    }
    if (dist == 'poisson') {
      lam <- exp(xb)
      out <- lam
    }
    out
  }

  x_mat = object@x
  y <- object@y
  x_mat <- as.matrix( x_mat[order(y),] )

  y <- sort(y)
  m <- sum(y==1)
  n <- length(y)
  if (dist == "binomial") K <- object@K
  N <- object@N
  beta <- as.matrix(object@beta)
  xb <- as.numeric(x_mat%*%beta)

  alp <- object@alpha
  f0 <- f_fun(rep(0,n), xb, dist)
  f1 <- f_fun(rep(1,n), xb, dist)
  varphi <- sum((1 - f0)^{-2})/N
  ef <- ef_fun(xb, dist)
  varf <- varf_fun(xb, dist)

  ### Vij's
  V11 <- 1 - 1/alp
  V13 <- 1/alp
  V22 <- t(x_mat) %*% diag(
    (f0 * ef^2/(1 - f0) - varf)/(1 - f0), n, n
  ) %*% x_mat/N

  V23 <- t(x_mat)%*%(f0*ef/(1-f0)^2)/N
  V24 <- (1 - alp)^2 * V23
  V33 <- varphi - 1/alp
  V34 <- (1 - alp)^2 * varphi
  V44 <- (1 - alp)^4 * varphi - (1 - alp)^3

  W <- rbind( cbind(-V11, matrix(0, 1, length(beta)),- V13),
              cbind(matrix(0, length(beta), 1),
                    -V22+V24%*%t(V24)/V44,
                    -V23+V24%*%t(V34)/V44),
              cbind(-V13,
                    -t(V23)+V34%*%t(V24)/V44,
                    -t(V33)+V34%*%t(V34)/V44 ))

  if (all(eigen(W)$values>0)) {

    W_chol <- forwardsolve( t(chol(W)), diag(1, nrow(W)))
    W_inv <- t(W_chol)%*%W_chol

    V_beta <- t(x_mat)%*%( (1 - f0 - ef)/ (1-f0))/N
    V_sid1 <- rbind(0, V_beta, 0)

    V25 <- t(x_mat)%*% ((ef - 1 + f0)/(1-f0))/N
    V55 <- sum((1-f0)/f1)/N - (1 - alp)
    V_sid2 <- rbind(0, V25, 0)
    varr <- t(V_sid1)%*%W_inv%*%V_sid1 + 2*t(V_sid1)%*%W_inv%*%V_sid2 + V55

    xbm <- x_mat[1:m,]%*%beta
    f01 <- f_fun(rep(0, m), xbm, dist)
    f11 <- f_fun(rep(1, m), xbm, dist)
    U <- n - sum((1 - f01)/f11)
    rt <- list(W_inv=W_inv, U=U, varr = varr)
  } else{
    xbm <- x_mat[1:m,]%*%beta
    f11 <- f_fun(rep(1, m), xbm, dist)
    U <- N - sum(1/f11)
    rt <- list(W_inv=NA, V_sid = NA, U=NA, varr = -1e9)
  }

  S <- as.numeric( rt$U/sqrt(rt$varr * N) )
  names(S) <- "scoretest_oizt"
  return(S)
}


