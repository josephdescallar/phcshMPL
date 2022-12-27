#' Fit Cause Specific Proportional Hazards Regression Model with a cure fraction
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
#' @export
#'
phcshcf_mpl2 <- function(formula, risk, z, data, control, ...){
  mc = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mc), 0)
  mc = mc[c(1, m)]
  if (m[1] == 0){
    stop("A formula argument is required")
  }
  else {
    "-"
  }
  mc[[1]] = as.name("model.frame")
  mc$formula = if (missing(data))
    stats::terms(formula)
  else stats::terms(formula, data = data)
  mf = eval(mc, parent.frame())
  mf_indx <- as.numeric(row.names(mf))
  y = stats::model.extract(mf, "response")
  if (nrow(mf) == 0)
    stop("No (non-missing) observations")
  mt = attr(mf, "terms")
  type = attr(y, "type")
  if (!inherits(y, "Surv")) {
    stop("Response must be a survival object")
  }
}
