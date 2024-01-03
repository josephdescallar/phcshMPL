#' Fit Cause Specific Proportional Hazards Regression Model Via MPL
#'
#' @description
#' Simultaneously estimate the regression coefficients and the baseline hazard
#' function of proportional hazard Cause specific hazards models for competing
#' risks with or without a cure fraction using maximum penalized likelihood
#' (MPL).
#'
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function (in the Surv, use type = "interval" with
#' interval censored data).
#'
#' @param risk a numeric vector which contains the risk/failure for the
#' i-th observation in the data frame.
#'
#' @param z an optional data frame which contains the covariate to be used when
#' a cure fraction is specified. If left NULL, then no cure fraction is assumed,
#' otherwise parameter estimation will occur for a cure fraction based on
#' the specified variables.
#'
#' @param data a data frame that includes the variables named in the formula
#' argument
#'
#' @param control test
#' @param ... Other arguments. In coxph_mpl, these elements, will be passed to
#' coxph_mpl_control. In print.coxph_mpl, these elements will be passed to the
#' print function.
#' @return A character vector.
#' @importFrom survival Surv
phcsh_mpl <- function(formula, data, z=NULL, risk, control, ...){
  if(is.null(z)){
    phcsh_mpl3(formula, risk, data, control, ...)
  }
  else phcshcf_mpl2(formula, risk, z, data, control, ...)
}
