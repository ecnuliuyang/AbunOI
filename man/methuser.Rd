\name{methuser}
\alias{methuser}
\docType{data}
\title{Female Metamphetamine User}
\description{
This data were collected in 2001, to estimate the number of female metamphetamine users in Bangkok, Thailand. It consists of the number of contacts a drug user has with the treatment institutions and the age information. For this dataset, 274 drug users have contacted the institutions at least once.

}
\usage{
data(methuser)
}

\format{
  A data frame consisting 274 observations and the following 2 variables.

\describe{
  \item{y}{A numeric vector, whose elements represent the number of contacts drug users have with the treatment institutions.}
  \item{age}{A numeric vector, whose elements represent the age of drug users.}
}

}

\details{
This data is used to demonstrate the performance of the empirical likelihood (EL) estimation approach associated with EM algorithms. By fitting a one-inflated poisson regression model to this data, we present the maximum EL estimator and the EL ratio confidence interval of the abundance parameter. The estimation results can be obtained by using the function \code{\link{abun_oi}}.

}

\source{
Bohning D. and van der Heijden, P. G. M. (2009). The covariate adjustment for zero-truncated approaches to estimate the size of hidden and elusive populations. \emph{The Annals of Applied Statistics} \bold{3}, 595--610.

}

\examples{
data(methuser)
head(methuser)
table(methuser$y)
summary(methuser$age)
}

\keyword{datasets}
