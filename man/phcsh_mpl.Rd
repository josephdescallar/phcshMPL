% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/phcsh_mpl.R
\name{phcsh_mpl}
\alias{phcsh_mpl}
\title{Fit Cause Specific Proportional Hazards Regression Model Via MPL}
\usage{
phcsh_mpl(formula, risk, data, control, ...)
}
\arguments{
\item{formula}{a formula object, with the response on the left of a ~
operator, and the terms on the right. The response must be a survival object
as returned by the Surv function (in the Surv, use type = "interval" with
interval censored data).}

\item{risk}{test}

\item{data}{test}

\item{control}{test}

\item{...}{Other arguments. In coxph_mpl, these elements, will be passed to
coxph_mpl_control. In print.coxph_mpl, these elements will be passed to the
print function.}
}
\value{
A character vector.
}
\description{
Simultaneously estimate the regression coefficients and the baseline hazard
function of proportional hazard Cause specific hazards models for competing
risks using maximum penalised likelihood (MPL).
}
