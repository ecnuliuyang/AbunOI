\name{prinia}
\alias{prinia}
\docType{data}
\title{Yellow-bellied Prinia}
\description{
This data were collected, over 17 weeks from January to April 1993, to estimate the prinia abundance at the Mai Po Bird Sanctuary in Hong Kong, China. It consists of the number of captures, the time of first being captured and wing lengths of the bird species \emph{Prinia flaviventris}. For this dataset, 164 birds were captured at least once.

}
\usage{
data(prinia)
}

\format{
  A data frame consisting 163 observations and the following 4 variables.

\describe{
  \item{i}{a numeric vector, bird identifier.}
  \item{y}{a numeric vector, the number of times the bird was captured.}
  \item{t1}{a numeric vector, the first time of bird being captured.}
  \item{x}{a numeric vector, wing length of captured birds.}
}

}

\details{
This data is used to demonstrate the performance of the empirical likelihood estimation approach associated with the EM algorithm. By fitting a one-inflated binomial regression model to this data, we present the maximum EL estimator and the EL ratio confidence interval of the abundance parameter. The estimation results can be obtained by invoking the function \code{\link{abun_oi}}.

}

\source{
Stoklosa, J., Hwang, W. H., Wu, S. H. and Huggins, R. M. (2011). Heterogeneous capture-recapture models with covariates: A partial likelihood approach for closed populations. \emph{Biometrics} \bold{67}, 1659--1665.

Liu, Y., Li, P. and Qin, J. (2017). Maximum empirical likelihood estimation for abundance in a closed population from capture-recapture data. \emph{Biometrika} \bold{104}, 527--543.
}


\keyword{datasets}
