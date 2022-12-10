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
  updscoretheta = function(theta, ti.S.gq.q, ti.gq.psi.qr, ti.rml.gq.w.psi,
                  ti.h0.q.l, xi.exb.qr){
    ti.S.gq.q.l = list()
    for(q in 1:n.risk){
    ti.S.gq.q.l[[q]] = if(n.i[[q]]!=0) {split(ti.S.gq.q[[q]],
    rep(1:ncol(ti.S.gq.q[[q]]), each =
    nrow(ti.S.gq.q[[q]])))}
    else NA
    }
    ti.B1.qr = ti.B2.qr = list()
    for(q in 1:n.risk){
      ti.B1.r = ti.B2.r = list()
      for(r in 1:n.risk){
        if(q==r){
        ti.B1.r[[r]] = if(n.i[[q]]!=0) {Reduce("+", mapply(function(a,b,c) a *
        as.vector(b) * as.vector(c), ti.gq.psi.qr[[q]][[r]],
        ti.S.gq.q.l[[q]], ti.rml.gq.w.l[[q]],
        SIMPLIFY = FALSE))}
        else NA
        }
        else{
        ti.B1.r[[r]] = if(n.i[[q]]!=0) {matrix(0, nrow = n.i[[q]],
        ncol = n.basis[[r]])}
        else NA
        }
        ti.B2.r[[r]] = if(n.i[[q]]!=0) {Reduce("+", mapply(function(a,b,c,d) a *
        b * c * d * as.vector(xi.exb.qr[[q]][[q]]), ti.h0.q.l[[q]],
        ti.gq.PSI.qr[[q]][[r]], ti.S.gq.q.l[[q]], ti.rml.gq.w.l[[q]],
        SIMPLIFY = FALSE))}
        else NA
      }
      ti.B1.qr[[q]] = ti.B1.r
      ti.B2.qr[[q]] = ti.B2.r
    }
    out = list(ti.S.gq.q.l = ti.S.gq.q.l, ti.B1.qr = ti.B1.qr,
          ti.B2.qr = ti.B2.qr)
    out
  }
  betascore.X <- betascore.1 <- betahess.X <- hess.AA.x <- list()
  theta.num.1 <- theta.den.1 <- list()
  for(q in 1:n.risk){
    betascore.X[[q]] = Reduce("rbind", list(if(length(tr)!=0) xr,
    if(n.e[[q]]!=0) xe[[q]], Reduce("rbind",
    mapply(function(a,b) if(a!=0) b, n.e, xe,
    SIMPLIFY=FALSE)), Reduce("rbind", mapply(function(a,b)
    if (a!=0) b, n.i, xi, SIMPLIFY = FALSE))))
    betascore.1[[q]] = rep(1, nrow(betascore.X[[q]]))
    betahess.X[[q]] = data.matrix(Reduce("rbind", list(if(length(tr)!=0) xr,
    Reduce("rbind", mapply(function(a,b) if(a!=0) b, n.e, xe,
    SIMPLIFY=FALSE)), Reduce("rbind", mapply(function(a,b)
    if(a!=0) b, n.i, xi, SIMPLIFY = FALSE)))))
    hess.AA.x[[q]] = data.matrix(Reduce("rbind", list(if(length(tr)!=0) xr,
    Reduce("rbind",mapply(function(a,b) if(a!=0) b, n.e, xe,
    SIMPLIFY = FALSE)), Reduce("rbind", mapply(function(a,b) if(a!=0) b, n.i,
     xi, SIMPLIFY = FALSE)))))
    theta.num.1[[q]] = matrix(1, ncol = 1, nrow = n.e[[q]] + sum(unlist(n.i)))
    theta.den.1[[q]] = matrix(1, ncol = 1, nrow = n.r + sum(unlist(n.e)) +
                                sum(unlist(n.i)))
  }
  hess.AA.x.diag = as.matrix(Matrix::bdiag(hess.AA.x))
  updhessian <- function(ti.S.gq.q, ti.h.q, gq.points, n.i, ti.gq.H0.qr,
               xi.exb.qr, ti.rml.gq.w, ti.gq.psi.qr, ti.gq.PSI.qr, tr.H.q,
               te.H.qr, ti.F.q, ti.A.qr, tr.PSI, te.PSI.qr, ti.B1.qr,
               ti.B2.qr, te.psi.qr, theta, xr.exb.q, xe.exb.qr, lambda, R.mat){
    ti.h.q.l = ti.gq.H.qr =  list()
    for(q in 1:n.risk){
      ti.h.q.l[[q]] = if(n.i[[q]]!=0) {split(ti.h.q[[q]], rep((1:gq.points),
      each = n.i[[q]]))}
      else NA
      ti.gq.H.r = list()
      for(r in 1:n.risk){ if(n.i[[q]]!=0){ti.gq.H.r[[r]] =
      ti.gq.H0.qr[[q]][[r]] * as.vector(xi.exb.qr[[q]][[r]])}
      else NA
      }
      ti.gq.H.qr[[q]] = ti.gq.H.r
    }

  ti.AA.qrjtk = ti.AB.qrjtu = ti.BB1.qrutz = ti.BB2.qrutz = list()
    for(q in 1:n.risk){
      ti.AA.rjtk = ti.AB.rjtu = ti.BB1.rutz = ti.BB2.rutz = list()
      for(r in 1:n.risk){
        ti.AA.jtk = ti.AB.jtu = ti.BB1.utz = ti.BB2.utz = list()
        for(t in 1:n.risk){
          if(q==r & r==t & q == t){
          ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
          ti.S.gq.q[[q]] * (1 - 3*ti.gq.H.qr[[q]][[r]] +
          ti.gq.H.qr[[q]][[r]]^2)) * ti.rml.gq.w[[q]])}
          else {NA}
          }
          else if(q==r & r!=t & q!=t){
            ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
            ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[t]] *(ti.gq.H.qr[[q]][[r]] - 1))
            * ti.rml.gq.w[[q]])}
            else {NA}
          }
          else if((q!=r & q==t & r!=t) | (q!=r & q!=t & r==t)){
          ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
          ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[r]] *(ti.gq.H.qr[[q]][[t]] - 1))
          * ti.rml.gq.w[[q]])}
          else {NA}
          }
          else if(q!=r & q!=t & r!=t){
          ti.AA.jtk[[t]] = if(n.i[[q]]!=0) {rowSums((ti.h.q[[q]] *
          ti.S.gq.q[[q]] * ti.gq.H.qr[[q]][[r]] * ti.gq.H.qr[[q]][[t]])
          * ti.rml.gq.w[[q]])}
          else {NA}
          }
          ti.AB.ju = list()
          for(w in 1:gq.points){
            if(q==r & q==t & r==t){
            ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
            (ti.S.gq.q[[q]][,w] * ti.gq.psi.qr[[q]][[t]][[w]] + ti.h.q[[q]][,w]*
            ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]] - 2 * ti.h.q[[q]][,w] *
            ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[t]][[w]] -
            ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
            ti.gq.psi.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
            }
            else {NA}
            }
            else if(q==r & q!=t & r!=t){
            ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
            (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]] - ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
            }
            else {NA}
            }
            else if(q!=r & q==t & r!=t){
            ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
            (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]] - ti.S.gq.q[[q]][,w] *
            ti.gq.H.qr[[q]][[r]][,w] * ti.gq.psi.qr[[q]][[t]][[w]]) *
            ti.rml.gq.w[[q]][,w]
            }
            else {NA}
            }
            else if(q!=r & q!=t & r==t){
            ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
            (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] * ti.gq.H.qr[[q]][[r]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]] - ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
            }
            else {NA}
            }
            else if(q!=r & q!=t & r!=t){
            ti.AB.ju[[w]] = if(n.i[[q]]!=0) {
            (ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
            ti.gq.PSI.qr[[q]][[t]][[w]]) * ti.rml.gq.w[[q]][,w]
            }
            else {NA}
            }
          }
          ti.AB.jtu[[t]] = Reduce("+", ti.AB.ju)
          ti.BB1.uz = ti.BB2.uz = list()
          for(u in 1:n.basis[[r]]){
            ti.BB1.z = ti.BB2.z = list()
            for(z in 1:n.basis[[t]]){
              ti.BB1 = ti.BB2 = list()
              for(w in 1:gq.points){
                if(n.i[[q]]!=0){
                  if(q==r & q==t & r==t){
                    ti.BB1[[w]] = -ti.S.gq.q[[q]][,w] *
                    ti.gq.psi.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z] *
                    ti.rml.gq.w[[q]][,w]
                    ti.BB2[[w]] = (ti.S.gq.q[[q]][,w] *
                    ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.psi.qr[[q]][[t]][[w]][,z] - ti.h.q[[q]][,w]
                    * ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if(q==r & q!=t & r!=t){
                    ti.BB1[[w]] = -ti.S.gq.q[[q]][,w] *
                    ti.gq.psi.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z] *
                    ti.rml.gq.w[[q]][,w]
                    ti.BB2[[w]] = (- ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                    ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if(q!=r & q==t & r!=t){
                    ti.BB1[[w]] = 0
                    ti.BB2[[w]] = (ti.S.gq.q[[q]][,w] *
                    ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.psi.qr[[q]][[t]][[w]][,z] - ti.h.q[[q]][,w]
                    * ti.S.gq.q[[q]][,w] * ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                  else if((q!=r & q!=t & r==t) | (q!=r & q!=t & r!=t)){
                    ti.BB1[[w]] = 0
                    ti.BB2[[w]] = (-ti.h.q[[q]][,w] * ti.S.gq.q[[q]][,w] *
                    ti.gq.PSI.qr[[q]][[r]][[w]][,u] *
                    ti.gq.PSI.qr[[q]][[t]][[w]][,z]) * ti.rml.gq.w[[q]][,w]
                  }
                }
                else{
                  ti.BB1[[w]] = NA
                  ti.BB2[[w]] = NA
                }
              }
              ti.BB1.z[[z]] = Reduce("+", ti.BB1)
              ti.BB2.z[[z]] = Reduce("+", ti.BB2)
            }
            ti.BB1.uz[[u]] = ti.BB1.z
            ti.BB2.uz[[u]] = ti.BB2.z
          }
          ti.BB1.utz[[t]] = ti.BB1.uz
          ti.BB2.utz[[t]] = ti.BB2.uz
        }
        ti.AA.rjtk[[r]] = ti.AA.jtk
        ti.AB.rjtu[[r]] = ti.AB.jtu
        ti.BB1.rutz[[r]] = ti.BB1.utz
        ti.BB2.rutz[[r]] = ti.BB2.utz
      }
      ti.AA.qrjtk[[q]] = ti.AA.rjtk
      ti.AB.qrjtu[[q]] = ti.AB.rjtu
      ti.BB1.qrutz[[q]] = ti.BB1.rutz
      ti.BB2.qrutz[[q]] = ti.BB2.rutz
    }
    te.h0.q = list()
    for(q in 1:n.risk){
      te.h0.q[[q]] = te.psi.qr[[q]][[q]] %*% theta[[q]]
    }
    hess.AA.rt.mat = hess.AB.rt.mat = hess.AB.x = list()
    for(r in 1:n.risk){
      hess.AA.t.mat = hess.AB.t.mat = list()
      for(t in 1:n.risk){
        hess.AA.t.mat[[t]] = Matrix::Diagonal(x = c(if(r==t){if(length(tr)!=0)
          -tr.H.q[[t]]} else{rep(0, n.r)},
          if(r==t){unlist(mapply(function(a,b) if(a!=0)
            -b[[t]], n.e, te.H.qr, SIMPLIFY = FALSE))}
          else{rep(0, sum(unlist(n.e)))},
          unlist(mapply(function(a,b,c,d) if(d!=0)
          {(a*b[[r]][[t]] - c[[r]]*c[[t]]) / (a^2)}, ti.F.q,
          ti.AA.qrjtk, ti.A.qr, n.i, SIMPLIFY = FALSE))))
      }
      hess.AA.rt.mat[[r]] = Reduce("rbind", hess.AA.t.mat)
    }
    hess.AA.mat = Reduce("cbind", hess.AA.rt.mat)
    hess.AA = t(hess.AA.x.diag) %*% hess.AA.mat %*% hess.AA.x.diag
    AA.rjtk = list()
    AA.mat = matrix(0, ncol = p*n.risk, nrow = p*n.risk)
    AA.count1 = 1
    AA.count2 = 1
    for(r in 1:n.risk){
      AA.jtk = list()
      for(j in 1:p){
        AA.tk = list()
        for(t in 1:n.risk){
          AA.k = list()
          for(k in 1:p){
            if(r==t){
              xj.r <- xr[,j]
              xk.r <- xr[,k]
              AA.tr = as.vector(-tr.H.q[[t]])
              xj.e = xj.i = xk.e = xk.i = AA.te = AA.ti = vector()
              for(q in 1:n.risk){
                if(n.e[[q]] != 0) {
                  xj.e = c(xj.e, xe[[q]][, j])
                  xk.e = c(xk.e, xe[[q]][, k])
                  AA.te = c(AA.te, -te.H.qr[[q]][[t]])
                }
                if(n.i[[q]] != 0){
                  xj.i = c(xj.i, xi[[q]][, j])
                  xk.i = c(xk.i, xi[[q]][, k])
                  AA.ti = c(AA.ti, as.vector((ti.F.q[[q]] *
                  ti.AA.qrjtk[[q]][[r]][[t]] - ti.A.qr[[q]][[r]] *
                  ti.A.qr[[q]][[t]]) / ti.F.q[[q]]^2))
                }
              }
            }
          }
        }
      }
    }
  }
}






















