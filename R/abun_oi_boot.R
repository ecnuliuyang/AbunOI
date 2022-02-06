##' Bootstrap-based standard errors and confidence intervals under discrete time CR models
##'
##' @description Function to calculate standard errors and confidence intervals based on a bootstrap procedure under discrete time capture-recapture (CR) models.
##'
##' @param object A \code{abun_oi} object.
##' @param boot number specifying the number of bootstrap samples. Default is 200.
##' @param seed number specifying the random seed. Default is 2021.
##' @return A list with four elements:
##'         \item{se_N}{the standard error of the abundance estimate.}
##'         \item{se_beta}{the standard error of the coefficients in count regression model.}
##'         \item{se_eta}{the standard error of the coefficients in one-inflated regression model.}
##'         \item{se_alpha}{the standard error of the probability of never being captured.}
##'         \item{quant_N}{the quantiles of the sampling distribution of the abundance estimate.}
##'
##' @references
##'
##' Liu, Y., Li, P., and Qin, J. (2017).
##' Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data.
##' \emph{Biometrika} \strong{104}, 527-543.
##'
##'
##' @importFrom stats plogis quantile rbinom rpois sd
##'
##' @export
##'
abun_oi_boot <- function (object, boot = 100, seed = 2021) {

  # object = out; boot = 100; seed = 2021
  rt <- list()
  model <- object@model
  dist <- object@dist

  x_true <- as.matrix( object@x )
  prob_true <- object@prob
  n <- length(prob_true)
  K <- object@K
  N_true <- round(object@N)
  beta_true <- as.matrix( object@beta )

  if ( model %in% c("oizt", "ztoi") ) {
    eta_true <- as.matrix( object@eta )
    z_true <- as.matrix( object@z )
  }

  set.seed(seed)
  N_b <- rep(NA, boot)
  beta_b <- matrix(NA, nrow = length(beta_true), ncol = boot)
  if ( model %in% c("oizt", "ztoi") )  eta_b <- matrix(NA, nrow = length(eta_true), ncol = boot)
  alpha_b <- rep(NA, boot)

  ### generate bootstrap samples
  data_gen <- function () {

    ind <- sample(1:n, N_true, replace = T, prob = prob_true)
    x <- as.matrix( x_true[ind, ] )
    xb <- as.numeric(x%*%beta_true)

    if (model == "zt") {
      y <- switch(dist, "binomial" = rbinom(N_true, K, plogis(xb)),
                  "poisson" = rpois(N_true, exp(xb)))
      rt <- list(y = y[y>0], x = x[y>0,])
    }

    if (model == "ztoi") {

      z <- as.matrix( z_true[ind, ] )
      w <- plogis( as.numeric(z%*%eta_true) )
      v <- rbinom(N_true, 1, w)
      y <- rep(1, N_true)

      for (i in 1:N_true)
        if (v[i] == 1) {
          y[i] <- switch(dist, "binomial" = rbinom(1, K, plogis(xb[i])),
                         "poisson" = rpois(1, exp(xb[i])))
        }

      rt <- list(y = y[y>0], x = x[y>0,], z = z[y>0,])
    }

    if (model == "oizt") {

      z <- as.matrix( z_true[ind, ] )
      w <- plogis( as.numeric(z%*%eta_true) )
      y <- switch(dist, "binomial" = rbinom(N_true, K, plogis(xb)),
                  "poisson" = rpois(N_true, exp(xb)))

      only_one <- rbinom(N_true, 1, 1 - w)
      y[y>0 & only_one==1] <- 1

      rt <- list(y = y[y>0], x = x[y>0,], z = z[y>0,])
    }

    return(rt)
  }

  b <- 1
  while (1) {

    cat("The bootstrap sample", b, "is generated for evaluation of standard error\n")

    dat_b <- data_gen()
    y_b <- dat_b$y
    x_b <- dat_b$x
    z_b <- dat_b$z


    out_b <- switch ( model,
                      "zt" = abun_zt ( y = y_b, K = K, x = x_b, dist = dist, maxN = object@maxN,
                                       start_beta = object@start_beta, eps = object@epsilon,
                                       maxit = object@maxit, N0 = NULL ),
                      "ztoi" = abun_ztoi ( y = y_b, K = K, x = x_b, z = z_b, dist = dist, maxN = object@maxN,
                                           start_beta = object@start_beta, start_eta = object@start_eta,
                                           eps = object@epsilon, maxit = object@maxit, N0 = NULL ),
                      "oizt" = abun_oizt ( y = y_b, K = K, x = x_b, z = z_b, dist = dist, maxN = object@maxN,
                                           start_beta = object@start_beta, start_eta = object@start_eta,
                                           eps = object@epsilon, maxit = object@maxit, N0 = NULL ) )
    N_b[b] <- out_b@N
    beta_b[,b] <- out_b@beta
    if (model == "ztoi" | model == "oizt")  eta_b[,b] <- out_b@eta
    alpha_b[b] <- out_b@alpha

    if (b == boot) break()
    b <- b+1
  }

  rt$se_N <- sd(N_b)
  rt$quant_N <- quantile(N_b, c(0.025, 0.05, 0.1, 0.9, 0.95, 0.975))
  rt$se_alpha <- sd(alpha_b)

  rt$se_beta <- apply(beta_b, 1, sd)
  if (model == "ztoi" | model == "oizt")  rt$se_eta <- apply(eta_b, 1, sd)

  return(rt)

}

