#' Fit Cause Specific Proportional Hazards Regression Model Via MPL
#'
#' Simultaneously estimate the regression coefficients and the baseline hazard
#' function of proportional hazard Cause specific hazards models for competing
#' risks using maximum penalised likelihood (MPL).
#'
#' @param formula a formula object, with the response on the left of a ~
#' operator, and the terms on the right. The response must be a survival object
#' as returned by the Surv function (in the Surv, use type = "interval" with
#' interval censored data).
#'
#' @param risk test
#' @param data test
#' @param control test
#' @param ... Other arguments. In coxph_mpl, these elements, will be passed to
#' coxph_mpl_control. In print.coxph_mpl, these elements will be passed to the
#' print function.
#' @return A character vector.
#' @export
#'
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
  if (!inherits(y, "Surv")) {
    stop("Response must be a survival object")
  }
  if (attr(y, which = "type") == "right") {
    left = y[, 1]
    right = rep(NA, nrow(y))
    icase = which(y[, 2] == 1)
    right[icase] = y[icase, 1]
    y = survival::Surv(left, right, type = "interval2")
  }
  else if (type != "interval") {
    stop("\nPlease create the survival object using the option type='interval2'
         in the Surv function.\n")
  }
  n = nrow(mf)
  extraArgs <- list(...)
  if(length(extraArgs)){
    controlArgs <- names(formals(phcsh_mpl_control))
    m <- pmatch(names(extraArgs), controlArgs, nomatch = 0L)
    if (any(m == 0L))
    stop(gettextf("Argument(s) %s not matched", names(extraArgs)[m ==
    0L]), domain = NA, call. = F)
  }
  if(missing(control))
    control <- phcsh_mpl_control(...)
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
  mean_j = apply(X, 2, mean)
  XC = X - rep(mean_j, each = n)
  risk = risk[index]
  data.export = data.frame(y,X,risk)
  data.all = data.frame(y,X,risk)
  data.all[,1][,2] = ifelse(data.all[,1][,3]==2,
                            data.all[,1][,1],data.all[,1][,2])
  data.all[,1][,1] = ifelse(data.all[,1][,3]==2, 0, data.all[,1][,1])
  l.risk <- sort(unique(risk))
  n.risk <- length(l.risk)
  tr = data.all[which(data.all[,1][,3]==0),1][,1]
  if(length(tr)==0) tr = NULL
  xr = data.all[which(data.all[,1][,3]==0),-c(1,ncol(data.all))]
  n.r = length(tr)
  t.min = min(data.all[,1][,c(1,2)])
  t.max = max(data.all[,1][,c(1,2)])
  b.knots <- c(t.min - 1e-12, t.max)
  knots.perc.limit = control$knots.perc.limit
  gq.points = control$gq.points
  dgr = control$dgr
  basis.intercept = control$basis.intercept
  gq <- statmod::gauss.quad(gq.points, "legendre")
  psif <- function(x, bknots, iknots){
    splines2::mSpline(x, knots = iknots, Boundary = bknots, degree = dgr,
            intercept = basis.intercept)
  }
  PSIf <- function(x, bknots, iknots)
    splines2::iSpline(x, knots = iknots, Boundary = bknots, degree = dgr,
            intercept = basis.intercept)
  te <- xe <- tl <- xl <- til <- tir <- xi <- n.e <- n.l <- n.i <- list()
  n.per.risk <- n.iknots <- i.knots <- perc.iknots <- t.all.r <- list()
  ti.rml <- ti.rpl <- ti.rml.gq <- ti.rml.gq.w <- ti.rml.gq.w.l <- list()
  ti.rml.gq.l <- ti.rml.gq.w.psi <- tmid <- list()
  n.basis <- R.mat <- oldbeta <- newbeta <- oldtheta <- newtheta <- list()
  tr.PSI <- te.psi <- oldlambda <- list()
  for(q in 1:n.risk){
    te[[q]] = data.all[which(risk==q & data.all[,1][,3]==1),1][,1]
    xe[[q]] = data.all[which(risk==q & data.all[,1][,3]==1),
                       -c(1,ncol(data.all))]
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
    xi[[q]] = data.all[which(risk==q & (data.all[,1][,3]==3 |
              data.all[,1][,3]==2)), -c(1,ncol(data.all))]
    n.i[[q]] = length(til[[q]])
    if(n.i[[q]] == 0){
      til[[q]] = NA
      tir[[q]] = NA
      xi[[q]] = NA
    }
    ti.rml[[q]] = (tir[[q]] - til[[q]]) / 2
    ti.rpl[[q]] = (tir[[q]] + til[[q]]) / 2
    if(control$tmid==TRUE)
      t.all.r[[q]] = c(te[[q]], tmid[[q]])
    else
      t.all.r[[q]] = c(te[[q]], til[[q]], tir[[q]])
    n.per.risk[[q]] = n.e[[q]] + n.i[[q]]
    n.iknots[[q]] = max(floor(n.per.risk[[q]]^(1/3)) - 2,1)
    if(length(control$n.basis)!=0){
      n.basis[[q]] = control$n.basis[[q]]
      n.iknots[[q]] = n.basis[[q]] - 2
    }
    n.basis[[q]] = n.iknots[[q]] + dgr + basis.intercept
    if(n.iknots[[q]]==1)
      perc.iknots[[q]] = 0.5
    else if(n.iknots[[q]]==2)
      perc.iknots[[q]] = seq(knots.perc.limit[1], knots.perc.limit[2],
                         length.out = 4)[2:3]
    else
      perc.iknots[[q]] = seq(knots.perc.limit[1], knots.perc.limit[2],
                             length.out = n.iknots[[q]])
    i.knots[[q]]  = stats::quantile(unique(t.all.r[[q]]),
                    prob = perc.iknots[[q]], names = FALSE, na.rm = TRUE,
                    type = 3)
    if(!is.null(control$iknots.pos)){
      i.knots[[q]] = control$iknots.pos[[q]]
    }
    ti.rml.gq[[q]] = sweep(matrix(ti.rml[[q]], nrow = n.i[[q]],
                     ncol = gq.points), MARGIN = 2, gq$nodes, `*`) +
                     ti.rpl[[q]]
    ti.rml.gq.w[[q]] = sweep(matrix(ti.rml[[q]], nrow = n.i[[q]],
                       ncol = gq.points), MARGIN = 2, gq$weights, `*`)
    ti.rml.gq.l[[q]] = split(ti.rml.gq[[q]], rep(1:gq.points,
                       each = n.i[[q]]))
    ti.rml.gq.w.l[[q]] = split(ti.rml.gq.w[[q]], rep(1:gq.points,
                         each = n.i[[q]]))
    ti.rml.gq.w.psi[[q]] = lapply(ti.rml.gq.l[[q]], function(a) psif(a,
                           b.knots, i.knots[[q]]))
    Rtemp = matrix(0, nrow = n.basis[[q]], ncol = n.basis[[q]])
    xknots = c(rep(t.min, dgr), i.knots[[q]], rep(t.max, dgr))
    for(ii in 1:n.basis[[q]]){
      for(jj in 1:n.basis[[q]]){
        if(jj - ii<dgr){
          kntset = xknots[xknots >= xknots[jj] & xknots <= xknots[ii+dgr]]
          kntsum = 0
          for(kk in 1:(length(kntset) - 1)){
            kntsum = kntsum + splines2::mSpline(kntset[kk], knots =
                     i.knots[[q]], degree = dgr, intercept = basis.intercept,
                     Boundary.knots = b.knots,
                     derivs = 2)[ii]*splines2::mSpline(
                     kntset[kk], knots = i.knots[[q]], degree = dgr,
                     intercept = basis.intercept, Boundary.knots = b.knots,
                     derivs = 2)[jj]*(kntset[kk + 1] - kntset[kk])
          }
          Rtemp[ii, jj] = kntsum
        }
      }
    }
    R.mat[[q]] <- Rtemp
    R.mat[[q]][lower.tri(R.mat[[q]], diag = FALSE)] =
      t(R.mat[[q]])[lower.tri(R.mat[[q]], diag = FALSE)]
    oldbeta[[q]] = rep(0, p)
    newbeta[[q]] = rep(0, p)
    oldtheta[[q]] = rep(1, n.basis[[q]])
    newtheta[[q]] = rep(1, n.basis[[q]])
    tr.PSI[[q]] = if(length(tr)!=0) PSIf(tr, b.knots, i.knots[[q]])
    else NA
    te.psi[[q]] = if(n.e[[q]]!=0) psif(te[[q]], b.knots, i.knots[[q]])
    else NA
    oldlambda[[q]] = control$lambda[q]
  }
  n.e.all = sum(unlist(n.e))
  n.i.all = sum(unlist(n.i))
  te.PSI.qr <- ti.gq.PSI.qr <- ti.gq.psi.qr <- te.psi.qr <- list()
  for(q in 1:n.risk){
    te.PSI.r <- ti.gq.PSI.r <- ti.gq.psi.r <- te.psi.r <- list()
    for(r in 1:n.risk){
      te.PSI.r[[r]] = if(n.e[[q]]!=0) PSIf(te[[q]], b.knots, i.knots[[r]])
      else NA
      ti.gq.PSI.r[[r]] = lapply(ti.rml.gq.l[[q]], function(a) PSIf(a,
                         b.knots, i.knots[[r]]))
      ti.gq.psi.r[[r]] = lapply(ti.rml.gq.l[[q]], function(a) psif(a,
                         b.knots, i.knots[[r]]))
      te.psi.r[[r]] = if(n.e[[q]]!=0) psif(te[[q]], b.knots, i.knots[[r]])
      else NA
    }
    te.PSI.qr[[q]] = te.PSI.r
    ti.gq.PSI.qr[[q]] = ti.gq.PSI.r
    ti.gq.psi.qr[[q]] = ti.gq.psi.r
    te.psi.qr[[q]] = te.psi.r
  }
  updbase <- function(theta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
                      ti.gq.PSI.qr){
    tr.H0.q = te.h0.q = ti.h0.q = ti.h0.q.l = te.H0.qr = ti.gq.H0.qr = list()
    ti.gq.S0r.qr = list()
    for(q in 1:n.risk){
      tr.H0.q[[q]] = tr.PSI[[q]] %*% theta[[q]]
      te.h0.q[[q]] = if(n.e[[q]]!=0) {te.psi[[q]] %*% theta[[q]]}
      else{NA}
      ti.h0.q[[q]] = sapply(ti.rml.gq.w.psi[[q]], function(a) a %*% theta[[q]])
      ti.h0.q.l[[q]] = split(ti.h0.q[[q]], rep(1:gq.points, each = n.i[[q]]))
      te.H0.r <- ti.gq.H0.r <- ti.gq.S0r.r <- list()
      for(r in 1:n.risk){
        te.H0.r[[r]] = if(n.e[[q]]!=0){te.PSI.qr[[q]][[r]] %*% theta[[r]]}
        else{NA}
        ti.gq.H0.r[[r]] = sapply(ti.gq.PSI.qr[[q]][[r]], function(a) a %*%
                          theta[[r]])
        ti.gq.S0r.r[[r]] = if(n.i[[q]]!=0) exp(-ti.gq.H0.r[[r]])
      }
      te.H0.qr[[q]] = te.H0.r
      ti.gq.H0.qr[[q]] = ti.gq.H0.r
      ti.gq.S0r.qr[[q]] = ti.gq.S0r.r
    }
    out = list(tr.H0.q=tr.H0.q, te.h0.q=te.h0.q, ti.h0.q=ti.h0.q,
               ti.h0.q.l=ti.h0.q.l, te.H0.qr = te.H0.qr, ti.gq.H0.qr =
               ti.gq.H0.qr, ti.gq.S0r.qr = ti.gq.S0r.qr)
    out
  }
  updllparms <- function(beta, xr, xe, tr.H0.q, xi, te.h0.q, te.H0.qr,
                         ti.h0.q, ti.gq.S0r.qr, ti.rml.gq.w){
    xr.exb.q  = xe.exb.qr = tr.H.q = xi.exb.qr = te.h.q = te.H.qr = list()
    ti.h.q = ti.Sr.gq.qr = ti.S.gq.q = ti.F.q = list()
    for(q in 1:n.risk){
      xr.exb.q[[q]] = exp(data.matrix(xr) %*% beta[[q]])
      tr.H.q[[q]] = tr.H0.q[[q]] * as.vector(xr.exb.q[[q]])
      xe.exb.r = xi.exb.r = te.H.r = ti.Sr.gq.r = list()
      for(r in 1:n.risk){
        xe.exb.r[[r]] = exp(data.matrix(xe[[q]]) %*% beta[[r]])
        xi.exb.r[[r]] = exp(data.matrix(xi[[q]]) %*% beta[[r]])
        te.H.r[[r]] = te.H0.qr[[q]][[r]] * as.vector(xe.exb.r[[r]])
        ti.Sr.gq.r[[r]] = if(n.i[[q]]!=0)
        ti.gq.S0r.qr[[q]][[r]] ^ as.vector(xi.exb.r[[r]])
      }
      xe.exb.qr[[q]] = xe.exb.r
      xi.exb.qr[[q]] = xi.exb.r
      te.h.q[[q]] = te.h0.q[[q]] * as.vector(xe.exb.qr[[q]][[q]])
      te.H.qr[[q]] = te.H.r
      ti.h.q[[q]] = if(n.i[[q]]!=0)
      ti.h0.q[[q]] * as.vector(xi.exb.qr[[q]][[q]])
      ti.Sr.gq.qr[[q]] = ti.Sr.gq.r
      ti.S.gq.q[[q]] = Reduce("*", ti.Sr.gq.r)
      ti.F.q[[q]] = if(n.i[[q]]!=0) {rowSums(ti.h.q[[q]] * ti.S.gq.q[[q]] *
                    ti.rml.gq.w[[q]])
      }
      else {NA}
    }
    out = list(xr.exb.q = xr.exb.q, xe.exb.qr = xe.exb.qr, tr.H.q = tr.H.q,
               xi.exb.qr = xi.exb.qr, te.h.q = te.h.q, te.H.qr = te.H.qr,
               ti.h.q = ti.h.q, ti.Sr.gq.qr = ti.Sr.gq.qr, ti.S.gq.q = ti.S.gq.q,
               ti.F.q = ti.F.q)
    out
  }
  updscorebeta <- function(ti.h.q, ti.S.gq.q, ti.gq.H0.qr, xi.exb.qr,
                           ti.rml.gq.w){
    ti.A.qr = list()
    for(q in 1:n.risk){
      ti.A.r = list()
      for(r in 1:n.risk){
        ti.A.r[[r]] = if(n.i[[q]]!=0) {rowSums(ti.h.q[[q]] * ti.S.gq.q[[q]] *
        ((q==r) - ti.gq.H0.qr[[q]][[r]]*as.vector(xi.exb.qr[[q]][[r]])) *
        ti.rml.gq.w[[q]])}
        else {NA}
      }
      ti.A.qr[[q]] = ti.A.r
    }
    out = list(ti.A.qr = ti.A.qr)
    out
  }
}






















