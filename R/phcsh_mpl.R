phcsh_mpl <- function(formula, risk, data, control, ...){
  mc = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mc), 0)
  mc = mc[c(1, m)]
  if (m[1] == 0) {
    stop("A formula argument is required")
  }
  data.name = if (m[2] != 0) {
    deparse(match.call()[[3]])
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
  # if (!inherits(y, "Surv")) {
  #   stop("Response must be a survival object")
  # }
  # if (attr(y, which = "type") == "right") {
  #   left = y[, 1]
  #   right = rep(NA, nrow(y))
  #   icase = which(y[, 2] == 1)
  #   right[icase] = y[icase, 1]
  #   y = survival::Surv(left, right, type = "interval2")
  # }
  # else if (type != "interval") {
  #   stop("\nPlease create the survival object using the option type='interval2'
  #        in the Surv function.\n")
  # }
  # n = nrow(mf)
  # extraArgs <- list(...)
  # if(length(extraArgs)){
  #   controlArgs <- names(formals(phcsh_mpl_control))
  #   m <- pmatch(names(extraArgs), controlArgs, nomatch = 0L)
  #   if (any(m == 0L))
  #   stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m ==
  #   0L]), domain = NA, call. = F)
  # }
}
