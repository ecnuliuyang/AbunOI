% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/abun_oi.R
\name{abun_oi}
\alias{abun_oi}
\title{Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data}
\usage{
abun_oi(
  formula,
  data,
  model = c("zt", "ztoi", "oizt"),
  dist = c("poisson", "binomial"),
  K = NULL,
  ci = FALSE,
  level = 0.95,
  control = abun_oi.control()
)
}
\arguments{
\item{formula}{symbolic description of the model, see 'Details'.}

\item{data}{data frame or list containing the variables in the model. If not found in data, the variables are taken from environment(formula).}

\item{model}{character specification of \code{"zt"} (zero-truncated model without one-inflation), \code{"ztoi"} (zero-truncated one-inflated model), or \code{"oizt"} (one-inflated zero-truncated model).}

\item{dist}{character specification of count regression model family, \code{"poisson"} or \code{"binomial"}.}

\item{K}{number specifying the number of capture occasions when \code{dist} is \code{"binomial"}.}

\item{ci}{logic. If TRUE, the empirical likelihood ratio conficence interval of abundance is calculated.}

\item{level}{number specifying the nominal level of confidence interval of abundance.}

\item{control}{list of control arguments in EM algorithm specified via \code{\link{abun_oi.control}}.}
}
\value{
An \code{abun_oi} object.
}
\description{
\code{abun_oi} is used to calculate the maximum empirical likelihood estimator and the empirical likelihood ratio confidence interval of abundance by fitting one-inflated count regression model.
}
\details{
If \code{model = "zt"}, the \code{formula} has the form \code{y ~ x} where \code{y} is the (numeric) vector representing the number of captures and \code{x} is a series of terms which specifies a linear predictor in count regression model.

If \code{model = "ztoi"} or \code{model = "oizt"}, the \code{formula} has the form \code{y ~ x|z} where \code{z} is a series of terms which specifies a linear predictor in one-inflated logistic regression model.
}
\examples{
### Estimation results for prinia data
### under zero-truncated binomial regression model without one-inflation
(pri_zt <- abun_oi(y ~ x, data = prinia, model = "zt",
                   dist = "binomial", K = 17))

### P-values of score tests for one-inflation
pri_st_ztoi <- scoretest_oi(pri_zt, model = "ztoi")
pri_st_oizt <- scoretest_oi(pri_zt, model = "oizt")
pnorm(pri_st_ztoi)
pnorm(pri_st_oizt)

### Under zero-truncated one-inflated binomial regression model
pri_ztoi <- abun_oi(y ~ x|1, data = prinia, model = "ztoi",
                    dist = "binomial", K = 17, ci = TRUE)
(pri_ztois <- summary(pri_ztoi, boot = 200))
### Maximum EL estimate of w
(w <- plogis(pri_ztois@eta))

### Under one-inflated zero-truncated binomial regression model
pri_oizt <- abun_oi(y ~ x|1, data = prinia, model = "oizt",
                    dist = "binomial", K = 17, ci = TRUE)
round(pri_oizt@N)
round(pri_oizt@ci)


}
\references{
Liu, Y., Li, P., Liu, Y., and Zhang, R. (2021).
Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data.
\emph{Biometrical Journal}.
}
