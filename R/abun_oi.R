##' Log-EL under the one-inflated count regression model
##'
##' @description Function to calculate the empirical log-likelihood (log-EL) function under the one-inflated count regression model.
##'
##' @param N number specifying the population size or abundance.
##' @param y vector specifying the number of times that individuals are captured.
##' @param x matrix specifying the individual covariates in count regression model.
##' @param beta vector specifying the coefficients in count regression model.
##' @param z matrix specifying the individual covariates in one-inflated logistic regression model.
##' @param eta vector specifying the regression coefficients in one-inflated logistic regression model.
##' @param alpha number specifying the probability of never being captured.
##' @param prob vector specifying the probability masses of individual covariates.
##' @param dist character specification of the count regression model family, \code{"poisson"} (default) or \code{"binomial"}.
##' @param K number specifying the number of capture occasions.
##' @param model character specification of zero-truncated one-inflated model (\code{"ztoi"}, default) or one-inflated zero-truncated model (\code{"oizt"}).
##'
##' @return A number, the empirical log-likelihood.
##'
##' @export
##'
loglikelihood_oi <- function (N, y, x, beta, z, eta, alpha, prob, dist, K, model) {

  switch (model, "ztoi" = loglikelihood_ztoi(N, y, x, beta, z, eta, alpha, prob, dist, K),
          "oizt" = loglikelihood_oizt(N, y, x, beta, z, eta, alpha, prob, dist, K))

}



##' Control parameters in EM algorithm
##'
##' @description Various parameters that control performing EM algorithm using \code{\link{abun_oi}}.
##'
##' @param epsilon positive convergence tolerance \eqn{\epsilon}. The iterations converge when the increase of the log-EL is less than \eqn{\epsilon}.
##' @param maxit integer specifying the maximal number of iterations.
##' @param maxN number specifying the largest value of \eqn{N} when calculating the estimate of the abundance \eqn{N}. If \code{NULL} the default is 100 times the number of individuals captured.
##' @param maxN_ci number specifying the largest value of \eqn{N} when calculating the confidence interval of \eqn{N}.
##' @param start_beta vector specifying the starting value of the coefficient \eqn{\beta} in count regression model. If \code{NULL} the default is zero.
##' @param start_eta vector specifying the starting value of the coefficient \eqn{\eta} in one-inflated logistic regression model. If \code{NULL} the default is zero.
##'
##' @return A list with the arguments specified.
##'
##' @export
##'
abun_oi.control <- function (epsilon = 1e-5, maxit = 5000, maxN = NULL, maxN_ci = 1e9, start_beta = NULL, start_eta = NULL) {

  if (!is.numeric(epsilon) || epsilon <= 0)
    stop("value of 'epsilon' must be > 0")
  if (!is.numeric(maxit) || maxit <= 0)
    stop("maximum number of iterations must be > 0")
  list(epsilon = epsilon, maxit = maxit, maxN = maxN, maxN_ci = maxN_ci,
       start_beta = start_beta, start_eta = start_eta)
}


##' Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data
##'
##' @description \code{abun_oi} is used to calculate the maximum empirical likelihood estimator and the empirical likelihood ratio confidence interval of abundance by fitting one-inflated count regression model.
##'
##' @param formula symbolic description of the model, see 'Details'.
##' @param data data frame or list containing the variables in the model. If not found in data, the variables are taken from environment(formula).
##' @param model character specification of \code{"zt"} (zero-truncated model without one-inflation), \code{"ztoi"} (zero-truncated one-inflated model), or \code{"oizt"} (one-inflated zero-truncated model).
##' @param dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @param K number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.
##' @param ci logic. If TRUE, the empirical likelihood ratio conficence interval of abundance is calculated.
##' @param level number specifying the nominal level of confidence interval of abundance.
##' @param control list of control arguments in EM algorithm specified via \code{\link{abun_oi.control}}.
##'
##' @details
##'
##' If \code{model = "zt"}, the \code{formula} has the form \code{y ~ x} where \code{y} is the (numeric) vector representing the number of captures and \code{x} is a series of terms which specifies a linear predictor in count regression model.
##'
##' If \code{model = "ztoi"} or \code{model = "oizt"}, the \code{formula} has the form \code{y ~ x|z} where \code{z} is a series of terms which specifies a linear predictor in one-inflated logistic regression model.
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
##' @importFrom stats coef glm optimize plogis model.matrix model.response terms nlminb pnorm printCoefmat
##'
##' @examples
##' ### Estimation results for prinia data
##' ### under zero-truncated binomial regression model without one-inflation
##' (pri_zt <- abun_oi(y ~ x, data = prinia, model = "zt",
##'                    dist = "binomial", K = 17))
##'
##' ### P-values of score tests for one-inflation
##' pri_st_ztoi <- scoretest_oi(pri_zt, model = "ztoi")
##' pri_st_oizt <- scoretest_oi(pri_zt, model = "oizt")
##' pnorm(pri_st_ztoi)
##' pnorm(pri_st_oizt)
##'
##' ### Under zero-truncated one-inflated binomial regression model
##' pri_ztoi <- abun_oi(y ~ x|1, data = prinia, model = "ztoi",
##'                     dist = "binomial", K = 17, ci = TRUE)
##' (pri_ztois <- summary(pri_ztoi, boot = 200))
##' ### Maximum EL estimate of w
##' (w <- plogis(pri_ztois@eta))
##'
##' ### Under one-inflated zero-truncated binomial regression model
##' pri_oizt <- abun_oi(y ~ x|1, data = prinia, model = "oizt",
##'                     dist = "binomial", K = 17, ci = TRUE)
##' round(pri_oizt@N)
##' round(pri_oizt@ci)
##'
##'
##' @export
##'
##'
abun_oi <- function ( formula, data, model = c("zt", "ztoi", "oizt"),
                      dist = c("poisson", "binomial"), K = NULL,
                      ci = FALSE, level = 0.95,
                      control = abun_oi.control() ) {

  model <- match.arg(model)
  dist <- match.arg(dist)
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data"), names(mf), 0L)
  mf <- mf[c(1L, m)]

  # print(model)

  if (model %in% c("ztoi", "oizt")) {

    ### step 1: obtain Y, X, Z

    if ( length(formula[[3]]) == 1) {
      stop("'formula' should be specified as '...|...'")
    } else if (!identical(formula[[3]][[1]], as.name("|"))) {
      stop("'formula' should be specified as '...|...'")
    }

    ### one-inflation depends on covariates
    if ( length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|"))) {
      ff <- formula
      formula[[3]][1] <- call("+")
      mf$formula <- formula
      ffcp <- . ~ .
      ffoi <- ~.
      ffcp[[2]] <- ff[[2]]
      ffcp[[3]] <- ff[[3]][[2]]
      ffoi[[3]] <- ff[[3]][[3]]
      ffoi[[2]] <- NULL
    }

    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())
    mtX <- terms(ffcp, data = data)
    X <- model.matrix(mtX, mf)

    mtZ <- terms(ffoi, data = data)
    Z <- model.matrix(mtZ, mf)

    Y <- as.matrix( model.response(mf, "numeric") )
    if (ncol(Y) > 1) stop("the number of times being captured in 'formula' should be a vector")
    if (dist == "binomial" & is.null(K)) stop("the number of capture occasions 'K' needs to be specified")


    ### step 2: fit data by zero-truncated one-inflated (ztoi) or one-inflated zero-truncated (oizt) model
    fit <- switch(model, "ztoi" = abun_ztoi ( y = Y, K = K, x = X, z = Z, dist = dist, maxN = control$maxN, start_beta = control$start_beta, start_eta = control$start_eta, eps = control$epsilon, maxit = control$maxit, N0 = NULL ),
                  "oizt" = abun_oizt ( y = Y, K = K, x = X, z = Z, dist = dist, maxN = control$maxN, start_beta = control$start_beta, start_eta = control$start_eta, eps = control$epsilon, maxit = control$maxit, N0 = NULL ))

    if ( ci )  fit@ci <- abun_oi_ci(fit, level = level, maxN_ci = control$maxN_ci)

  } else {

    ### step 1: obtain Y, X, Z
    mf$formula <- formula
    ffcp <- . ~ .
    ffcp[[2]] <- formula[[2]]
    ffcp[[3]] <- formula[[3]]
    mf[[1]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mtX <- terms(ffcp, data = data)
    X <- model.matrix(mtX, mf)
    Y <- as.matrix( model.response(mf, "numeric") )
    if (ncol(Y) > 1) stop("the number of times being captured in 'formula' should be a vector")
    if (dist == "binomial" & is.null(K)) stop("the number of capture occasions 'K' needs to be specified")


    ### step 2: fit data by zero-truncated one-inflated (ztoi) or one-inflated zero-truncated (oizt) model
    fit <- abun_zt ( y = Y, K = K, x = X, dist = dist, maxN = control$maxN, start_beta = control$start_beta, eps = control$epsilon, maxit = control$maxit, N0 = NULL )

    if ( ci ) fit@ci <- abun_oi_ci(fit, level = level, maxN_ci = control$maxN_ci)

  }

  fit

}
