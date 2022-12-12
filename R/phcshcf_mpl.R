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
phcshcf_mpl <- function(formula, risk, z, data, control, ...){
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
  if (attr(y, which = "type") == "right") {
    left = y[, 1]
    right = rep(NA, nrow(y))
    icase = which(y[, 2] == 1)
    right[icase] = y[icase, 1]
    y = Surv(left, right, type = "interval2")
  }
  else if (type != "interval") {
    stop("\nPlease create the survival object using the option type='interval2'
         in the Surv function.\n")
  }
  n = nrow(mf)
  extraArgs <- list(...)
  if(length(extraArgs)){
    controlArgs <- names(formals(phcshcf_mpl_control))
    m <- pmatch(names(extraArgs), controlArgs, nomatch = 0L)
    if (any(m == 0L))
      stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m ==
      0L]), domain = NA, call. = F)
  }
  if(missing(control))
    control <- phcshcf_mpl_control(...)
  index = as.vector(row(mf)[,1])
  X = stats::model.matrix(mt, mf)
  X = X[, !apply(X, 2, function(x) all(x == x[1])), drop = FALSE]
  if (ncol(X) == 0) {
    X = matrix(0, n, 1)
    noX = TRUE
  }
  else {
    noX = FALSE
  }
  p = ncol(X)
  xnames <- colnames(X)
  mean_j = apply(X, 2, mean)
  XC = X - rep(mean_j, each = n)
  risk = risk[index]
  z = z[index,]
  Z = z[, !apply(z, 2, function(z) all(z == z[1])), drop = FALSE]
  Z = cbind(rep(1,n),Z)
  colnames(Z) <- c("z0","z1","z2") #HARD CODED, FIX THIS!!
  znames <- colnames(Z)
  o <- ncol(Z)
  mean_o = c(0,apply(Z[,-1], 2, mean))
  ZC = Z - rep(mean_o, each = n)
  data.export = data.frame(y,X,Z,risk)
  data.all = data.frame(y,XC,risk,ZC)
  data.all[,1][,2] = ifelse(data.all[,1][,3]==2,
                     data.all[,1][,1],data.all[,1][,2])
  data.all[,1][,1] = ifelse(data.all[,1][,3]==2, 0, data.all[,1][,1])
  l.risk <- sort(unique(risk))
  n.risk <- length(l.risk)
  tr = data.all[which(data.all[,1][,3]==0),1][,1]
  if(length(tr)==0) tr = NULL
  xr = as.matrix(data.all[which(data.all[,1][,3]==0), xnames])
  zr = as.matrix(data.all[which(data.all[,1][,3]==0), znames])
  n.r = length(tr)
  t.min = min(data.all[,1][,c(1,2)])
  t.max = max(data.all[,1][,c(1,2)])
  b.knots <- c(t.min - 1e-12, t.max)
  knots.perc.limit = control$knots.perc.limit
  gq.points = control$gq.points
  dgr = control$dgr
  basis.intercept = control$basis.intercept
  gq <- statmod::gauss.quad(gq.points, "legendre")
  psif <- function(x, bknots, iknots){splines2::mSpline(x, knots = iknots,
          Boundary = bknots, degree = dgr, intercept = basis.intercept)}
  PSIf <- function(x, bknots, iknots) {splines2::iSpline(x, knots = iknots,
          Boundary = bknots, degree = dgr,intercept = basis.intercept)}
  te <- xe <- tl <- xl <- til <- tir <- xi <- n.e <- n.l <- n.i <- list()
  n.per.risk <- n.iknots <- i.knots <- perc.iknots <- t.all.r <- list()
  ti.rml <- ti.rpl <- ti.rml.gq <- ti.rml.gq.w <- ti.rml.gq.w.l <- list()
  ti.rml.gq.l <- ti.rml.gq.w.psi <- tmid <- list()
  n.basis <- R.mat <- oldbeta <- newbeta <- oldtheta <- newtheta <- list()
  tr.PSI <- te.psi <- oldlambda <- oldgamm <- newgamm <- ze <- zi <- list()
  for(q in 1:n.risk){
    te[[q]] = data.all[which(risk==q & data.all[,1][,3]==1),1][,1]
    xe[[q]] = as.matrix(data.all[which(risk==q & data.all[,1][,3]==1), xnames])
    ze[[q]] <- as.matrix(data.all[which(risk==q & data.all[,1][,3]==1), znames])
    n.e[[q]] = length(te[[q]])
    if(n.e[[q]] == 0) {
      te[[q]] = NA
      xe[[q]] = NA
    }
    til[[q]] = data.all[which(risk==q & (data.all[,1][,3]==3 |
               data.all[,1][,3]==2)),1][,1]
    tir[[q]] = data.all[which(risk==q & (data.all[,1][,3]==3 |
               data.all[,1][,3]==2)),1][,2]
    tmid[[q]] = rowMeans(cbind(til[[q]],tir[[q]]))
    xi[[q]] = as.matrix(data.all[which(risk==q & (data.all[,1][,3]==3 |
              data.all[,1][,3]==2)), xnames])
    zi[[q]] = as.matrix(data.all[which(risk==q & (data.all[,1][,3]==3 |
              data.all[,1][,3]==2)), znames])
    n.i[[q]] = length(til[[q]])
    if(n.i[[q]] == 0){
      til[[q]] = NA
      tir[[q]] = NA
      xi[[q]] = NA
    }
  }
}





































