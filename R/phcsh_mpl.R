#' Fit Cause Specific Proportional Hazards Regression Model
#'  Via MPL
#'
#' Simultaneously estimate the regression coefficients and the baseline hazard
#' function of proportional hazard Cause specific hazards models for competing
#' risks with a cure fraction for at least one risk using maximum penalised
#' likelihood (MPL).
#'
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function (in the Surv, use type = "interval" with
#' interval censored data).
#'
#' @param risk test
#' @param z test
#' @param data test
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
