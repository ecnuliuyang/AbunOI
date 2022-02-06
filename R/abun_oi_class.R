##' abun_oi class
##'
##' @description A S4 class is used to present the data, model specification, and maximum empirical likelihood (EL) estimation results under one-inflated count regression models.
##'
##' @slot model character specification of \code{"zt"} (zero-truncated model without one-inflation), \code{"ztoi"} (zero-truncated one-inflated model), or \code{"oizt"} (one-inflated zero-truncated model).
##' @slot dist character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.
##' @slot N number representing the maximum EL estimate of abundance \eqn{N}.
##' @slot ci vector representing the EL ratio confidence interval of abundance.
##' @slot beta vector representing the maximum EL estimate of regression coefficients in count regression model.
##' @slot eta vector representing the maximum EL estimate of regression coefficients in logistic regression model for one-inflation.
##' @slot alpha number representing the maximum EL estimate of the probability of never being captured.
##' @slot loglikelihood number representing the log-EL value.
##' @slot AIC number representing the Akaike Information Criterion value.
##' @slot prob vector representing the probability masses of individual covariates.
##' @slot nit number representing the number of iterations of EM algorithm.
##' @slot pars_trace matrix. Row shows the value of parameters in each iteration of EM algorithm.
##' @slot loglikelihood_trace vector. Element represents the value of log-likelihood in each iteration of EM algorithm.
##' @slot y vector specifying the number of times that individuals were captured.
##' @slot K number specifying the number of capture occasions.
##' @slot x matrix specifying the individual covariates in count regression model.
##' @slot z matrix specifying the individual covariates in one-inflated regression model.
##' @slot epsilon positive convergence tolerance \eqn{\epsilon}. The iterations converge in EM algorithm when the increase of the log-EL is less than \eqn{\epsilon}.
##' @slot maxN number specifying the largest searching value for \eqn{N} in EM algorithm.
##' @slot maxN_ci number specifying the largest searching value for calculating the confidence interval of \eqn{N}.
##' @slot level number specifying the nominal level of confidence interval of \eqn{N}.
##' @slot maxit integer specifying the maximal number of iterations in EM algorithm.
##' @slot start_beta vector specifying the starting value of the coefficient \eqn{\beta} in count regression model in EM algorithm.
##' @slot start_eta vector specifying the starting value of the coefficient \eqn{\eta} in one-inflated regression model in EM algorithm.
##' @export
##'
##'
setClass( 'abun_oi', slots =
            list( model = "character",
                  dist = "character",
                  N = "numeric",
                  ci = "numeric",
                  beta = "numeric",
                  eta = "numeric",
                  alpha = "numeric",
                  loglikelihood = "numeric",
                  AIC =  "numeric",
                  prob = "numeric",
                  nit = "numeric",
                  pars_trace = "matrix",
                  loglikelihood_trace = "numeric",
                  y = "numeric",
                  K = "numeric",
                  x = "matrix",
                  z = "matrix",
                  epsilon = "numeric",
                  maxN = "numeric",
                  maxN_ci = "numeric",
                  level = "numeric",
                  maxit = "numeric",
                  start_beta = "numeric",
                  start_eta = "numeric") )



##' Function to show the outcomes of object "abun_oi"
##' @param object an \code{abun_oi} object
setMethod('show', 'abun_oi',
          function(object) {
            digits = max(3L, getOption("digits") - 3L)

            switch(object@model,
                   "ztoi" = cat("\nMaximum EL estimation under zero-truncated one-inflated regression model\n"),
                   "oizt" = cat("\nMaximum EL estimation under one-inflated zero-truncated regression model\n"),
                   "zt" = cat("\nMaximum EL estimation under zero-truncated regression model without one-inflation\n"))

            est_beta <- object@beta
            names(est_beta) <- colnames(object@x)
            if (object@model == "ztoi" | object@model == "oizt")  {
              est_eta <- object@eta
              if (all(object@z == 1)) {
                names(est_eta) <- "(Intercept)"
              } else {
                names(est_eta) <- colnames(object@z)
              }
            }

            est_N <- object@N
            names(est_N) <- "Estimate"

            if (length(object@ci) > 0) {
              est_N <- c( est_N, object@ci )
              names(est_N) <- c( "Estimate", "Lower bound", "Upper bound" )
            }

            cat("\nCoefficients in count regression model:\n")
            print(est_beta)

            if (object@model == "ztoi" | object@model == "oizt") {
              cat("\nCoefficients in one-inflated regression model:\n")
              print(est_eta)
            }

            cat("\nCoefficients in one-inflated regression model:\n")
            print(est_N)
            invisible(object)

          })






##' abun_oi.summary
##' @description a S4 class summarizes the point estimates under discrete time capture-recapture models
##' @slot coefficients a list, which summarizes the estimates of coefficients in count regression model and one-inflated regression model
##' @slot abundance a matrix, which summarizes the estimates of abundance
##' @export
setClass("abun_oi.summary", contains="abun_oi",
         slots=list(coefficients = "list",
                    abundance = "numeric"))


##' Function to summarize the maximum empirical likelihood estimation results
##' @param object an \code{abun_oi} object
##'
##' @param boot integer specifying the number of bootstrap replications when evaluating the variation.
##' @param seed integer specifying the random seed when invoking the bootstrap method \code{\link{abun_oi_boot}}.
##'
##' @importFrom methods new
##' @importFrom stats pnorm
##' @export
setMethod('summary', 'abun_oi',
          function(object, boot = 100, seed = 2021) {
            rt <- new('abun_oi.summary',
                      model = object@model,
                      dist = object@dist,
                      N = object@N,
                      ci = object@ci,
                      beta = object@beta,
                      eta = object@eta,
                      alpha = object@alpha,
                      loglikelihood = object@loglikelihood,
                      AIC =  object@AIC,
                      prob = object@prob,
                      nit = object@nit,
                      pars_trace = object@pars_trace,
                      loglikelihood_trace = object@loglikelihood_trace,
                      y = object@y,
                      K = object@K,
                      x = object@x,
                      z = object@z,
                      epsilon = object@epsilon,
                      maxN = object@maxN,
                      maxN_ci = object@maxN_ci,
                      level = object@level,
                      maxit = object@maxit,
                      start_beta = object@beta,
                      start_eta = object@eta)

            se <- abun_oi_boot(object, boot = boot, seed = seed)
            # print(se)
            rt@coefficients <- list()
            est_beta <- object@beta
            se_beta <- se$se_beta
            zval_beta <- est_beta/abs(se_beta)
            coef_beta <- cbind(est_beta, se_beta,  zval_beta, 2L * pnorm(abs(zval_beta), lower.tail = FALSE))
            coef_beta <- matrix(coef_beta, ncol = 4 )
            colnames(coef_beta) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
            rownames(coef_beta) <- colnames(object@x)  # paste0("beta", 1L:length(est_beta))
            rt@coefficients$beta <- coef_beta

            if (object@model == "ztoi" | object@model == "oizt") {
              est_eta <- object@eta
              se_eta <- se$se_eta
              zval_eta <- est_eta/abs(se_eta)
              coef_eta <- cbind(est_eta, se_eta,  zval_eta, 2L * pnorm(abs(zval_eta), lower.tail = FALSE))
              coef_eta <- matrix(coef_eta, ncol = 4L )
              colnames(coef_eta) <- c('Estimate', 'Std. Error', 'z value', 'Pr(>|z|)')
              if (all(object@z == 1)) {
                rownames(coef_eta) <- "(Intercept)"
              } else {
                rownames(coef_eta) <- colnames(object@z)
              }
              rt@coefficients$eta <- coef_eta
            }

            abundance <- c(round(object@N), round(se$se_N, 2L))
            names(abundance) <- c("Estimate", "Std. Error" )

            if ( length(object@ci) > 0 ) {
              abundance <- c(abundance, object@ci)
              names(abundance) <- c("Estimate", "Std. Error", "ci_lower", "ci_upper")
            }
            rt@abundance <- abundance

            return(rt)
          })


##' Function to show the outcomes of object "abun_oi.summary"
##' @param object an \code{abun_oi} object
##' @importFrom stats printCoefmat
setMethod('show', 'abun_oi.summary',
          function(object) {
            digits = max(3L, getOption("digits") - 3L)
            signif.stars = getOption("show.signif.stars")

            # cat("\n-------------------------------------------------------------------\n")

            switch(object@model,
                   "ztoi" = cat("\nMaximum EL estimation under zero-truncated one-inflated regression model\n"),
                   "oizt" = cat("\nMaximum EL estimation under one-inflated zero-truncated regression model\n"),
                   "zt" = cat("\nMaximum EL estimation under zero-truncated regression model without one-inflation\n"))

            cat("\nCoefficients in count regression model:\n")
            printCoefmat(object@coefficients$beta, digits = digits, signif.stars = signif.stars, na.print = "NA")

            if (object@model == "ztoi" | object@model == "oizt") {
              cat("\nCoefficients in one-inflated regression model:\n")
              printCoefmat(object@coefficients$eta, digits = digits, na.print = "NA")
            }

            cat("\nAbundance:\n")
            print(object@abundance)

            cat ('\nLog-likelihood:', object@loglikelihood,
                 'with AIC of', round(object@AIC, 2L),
                 '\n\n')
            invisible(object)
          })
