# AbunOI

The R package **AbunOI** implements the semiparametric empirical likelihood (EL) inference for abundance via EM algorithms. 

The zero-truncated, zero-truncated one-inflated, and one-inflated zero-truncated binomial and Poisson regression models could be used to fit the one-inflated capture-recapture data. Under these models, this package provides four main functions to make statistical inferences.

+ abun_oi(): used to calculate the maximum EL estimators of abundance and regression coefficients, as well as the EL ratio confidence interval of abundance.

- abun_oi_test(): used to calculate the score test statistics for whether the one-inflation exists.

* abun_oi_boot(): used to calculate the bootstrap-based standard errors for the maximum EL estimators.

+ summary(): used to summarize the maximum EL estimation results.


# Usage

In the **R** software, the following codes are used to install the package:

install.packages("devtools")

library(devtools)

install_github("ecnuliuyang/AbunOI")



As an example, the following codes are provided to show the usage of the main functions in this packag: 

library(AbunOI)

example(abun_oi)


# Reference
Liu, Y., Li, P., Liu, Y., and Zhang, R. (2021). Semiparametric empirical likelihood inference for abundance from one-inflated capture-recapture data. *Biometrical Journal*.


#

For questions, comments or remarks about the code please contact Y. Liu at this email address <liuyangecnu@163.com>.
