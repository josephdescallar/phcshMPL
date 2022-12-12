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
    i.knots[[q]]  = stats::quantile(unique(t.all.r[[q]]), prob =
                    perc.iknots[[q]], names = FALSE, na.rm = TRUE, type = 3)
    if(!is.null(control$iknots.pos)){
      i.knots[[q]] = control$iknots.pos[[q]]
    }
    ti.rml.gq[[q]] = sweep(matrix(ti.rml[[q]], nrow = n.i[[q]],
                     ncol = gq.points), MARGIN = 2, gq$nodes, `*`) + ti.rpl[[q]]
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
            kntsum = kntsum + splines2::mSpline(kntset[kk],
                     knots = i.knots[[q]], degree = dgr, intercept =
                     basis.intercept, Boundary.knots = b.knots, derivs =
                     2)[ii]*splines2::mSpline(kntset[kk], knots = i.knots[[q]],
                     degree = dgr, intercept = basis.intercept,
                     Boundary.knots = b.knots, derivs = 2)[jj]*(kntset[kk + 1] -
                     kntset[kk])
          }
          Rtemp[ii, jj] = kntsum
        }
      }
    }
    R.mat[[q]] <- Rtemp
    R.mat[[q]][lower.tri(R.mat[[q]], diag = FALSE)] =
      t(R.mat[[q]])[lower.tri(R.mat[[q]], diag = FALSE)]
    oldgamm[[q]] <- rep(0, o)
    newgamm[[q]] <- rep(0, o)
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
  max.outer = control$max.outer
  max.iter = control$max.iter
  # updbase <- function(theta, tr.PSI, te.psi, te.PSI.qr, ti.gq.PSI.qr,
  #                     ti.rml.gq.w.psi){
  #   te.h0.qr <- ti.gq.H0.qr <- ti.gq.S0r.qr <- ti.h0.q <- ti.h0.q.l <- list()
  #   te.h0.q <- tr.H0r.r <- te.H0r.qr <- te.S0r.qr <- list()
  #   for(q in 1:n.risk){
  #     te.h0.r <- ti.gq.H0.r <- ti.gq.S0r.r <- te.H0r.r <- te.S0r.r <- list()
  #     te.h0.q[[q]] <- te.psi[[q]] %*% theta[[q]]
  #     for(r in 1:n.risk){
  #       ti.gq.H0.r[[r]] <- sapply(ti.gq.PSI.qr[[q]][[r]], function(a) a %*%
  #                          theta[[r]])
  #       ti.gq.S0r.r[[r]] <- if(n.i[[q]]!=0) exp(-ti.gq.H0.r[[r]])
  #       te.H0r.r[[r]] <- te.PSI.qr[[q]][[r]] %*% theta[[r]]
  #       te.S0r.r[[r]] <- exp(-te.H0r.r[[r]])
  #     }
  #     tr.H0r.r[[q]] <- tr.PSI[[q]] %*% theta[[q]]
  #     ti.gq.H0.qr[[q]] <-  ti.gq.H0.r
  #     ti.gq.S0r.qr[[q]] <- ti.gq.S0r.r
  #     ti.h0.q[[q]] <- sapply(ti.rml.gq.w.psi[[q]], function(a) a %*% theta[[q]])
  #     ti.h0.q.l[[q]] <- split(ti.h0.q[[q]], rep(1:gq.points, each = n.i[[q]]))
  #     te.H0r.qr[[q]] <- te.H0r.r
  #     te.S0r.qr[[q]] <- te.S0r.r
  #   }
  #   out <- list("te.h0.q"=te.h0.q, "ti.gq.S0r.qr"=ti.gq.S0r.qr,
  #               "ti.h0.q"=ti.h0.q, "ti.h0.q"=ti.h0.q, "tr.H0r.r"=tr.H0r.r,
  #               "te.H0r.qr"=te.H0r.qr, "te.S0r.qr"=te.S0r.qr)
  # }
  # updllparms <- function(gamm, ze, zi, te.h0.q, beta, ti.h0.q, ti.gq.S0r.qr,
  #                        xi, ti.rml.gq.w, tr.H0r.r, xr, zr, te.S0r.qr){
  #   ze.exbzg.r <- ze.pi.r <- zi.ezg.qr <- zi.pi.r <- te.h.qr <- list()
  #   xe.exb.qr <- te.h.q <- xi.exb.qr <- ti.Sr.gq.qr <- ti.h.q <- list()
  #   ti.Sr.gq.qr.pop <- ti.S.gq.q <- ti.F.q <- ti.pi.qr <- xe.exb.q <- list()
  #   xi.exb.q <- zi.ezg.q <- zi.pi.q <- tr.H.r <- xr.exb.r <- tr.Sr.r <- list()
  #   tr.Sr.pop.r <- tr.Hr.r <- tr.Sr.pop.r <- zr.ezg.q <- zr.pi.q <- list()
  #   te.Sr.qr <- te.Sr.pop.qr <- ti.S.gq.q.pop <- list()
  #   for(q in 1:n.risk){
  #     xe.exb.r <- te.h.r <- xi.exb.r <- ti.Sr.gq.r <- ti.Sr.gq.r.pop <- list()
  #     ti.pi.r <- zi.ezg.r <- te.Sr.r <- te.Sr.pop.r <- list()
  #     for(r in 1:n.risk){
  #       ze.exbzg.r[[r]] <- exp(ze[[q]] %*% gamm[[r]])
  #       ze.pi.r[[r]] <- ze.exbzg.r[[r]] / (1 + ze.exbzg.r[[r]])
  #       zi.ezg.r[[r]] <- exp(zi[[q]] %*% gamm[[r]])
  #       zi.pi.r[[r]] <- zi.ezg.r[[r]] / (1 + zi.ezg.r[[r]])
  #       xe.exb.r[[r]] <- exp(data.matrix(xe[[q]]) %*% beta[[r]])
  #       xi.exb.r[[r]] = exp(data.matrix(xi[[q]]) %*% beta[[r]])
  #       ti.pi.r[[r]] <- zi.ezg.r[[r]] / (1 + zi.ezg.r[[r]])
  #       ti.Sr.gq.r[[r]] = if(n.i[[q]]!=0)
  #         ti.gq.S0r.qr[[q]][[r]] ^ as.vector(xi.exb.r[[r]])
  #     }
  #   }
  # }
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
                         ti.h0.q, ti.gq.S0r.qr, ti.rml.gq.w, gamm, ze, zi){
    xr.exb.q  = xe.exb.qr = tr.H.q = xi.exb.qr = te.h.q = te.H.qr = list()
    ti.h.q = ti.Sr.gq.qr = ti.S.gq.q = ti.F.q = pi.ze = pi.zi = list()
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
      pi.ze[[q]] <- exp(data.matrix(ze[[q]]) %*% gamm[[q]]) / (1 +
                     exp(data.matrix(ze[[q]]) %*% gamm[[q]]))
      pi.zi[[q]] <- exp(data.matrix(zi[[q]]) %*% gamm[[q]]) / (1 +
                     exp(data.matrix(zi[[q]]) %*% gamm[[q]]))
    }
    out = list(xr.exb.q = xr.exb.q, xe.exb.qr = xe.exb.qr, tr.H.q = tr.H.q,
               xi.exb.qr = xi.exb.qr, te.h.q = te.h.q, te.H.qr = te.H.qr,
               ti.h.q = ti.h.q, ti.Sr.gq.qr = ti.Sr.gq.qr, ti.S.gq.q = ti.S.gq.q,
               ti.F.q = ti.F.q, pi.ze=pi.ze, pi.zi=pi.zi)
    out
  }
  updscoregamm <- function(gamm, ze, zi){
    ze.C.r <- zi.C.r <- ze.ezg <- zi.ezg <- list()
    for(r in 1:n.risk){
      ze.ezg[[r]] <- exp(data.matrix(ze[[r]]) %*% gamm[[r]])
      zi.ezg[[r]] <- exp(data.matrix(zi[[r]]) %*% gamm[[r]])
      ze.C.r[[r]] <- exp(data.matrix(ze[[r]]) %*% gamm[[r]]) /
      (1 + exp(data.matrix(ze[[r]]) %*% gamm[[r]]))^2
      zi.C.r[[r]] <- exp(data.matrix(zi[[r]]) %*% gamm[[r]]) /
        (1 + exp(data.matrix(zi[[r]]) %*% gamm[[r]]))^2
    }
    out = list(ze.C.r=ze.C.r, zi.C.r=zi.C.r, ze.ezg=ze.ezg, zi.ezg=zi.ezg)
    out
  }
  updscorebeta <- function(ti.h.q, ti.S.gq.q, ti.gq.H0.qr, xi.exb.qr,
                  ti.rml.gq.w){
    ti.A.qr = list()
    for(q in 1:n.risk){
      ti.A.r = list()
      for(r in 1:n.risk){
        ti.A.r[[r]] = if(n.i[[q]]!=0) {rowSums(ti.h.q[[q]] * ti.S.gq.q[[q]] *
                      ((q==r) - ti.gq.H0.qr[[q]][[r]] *
                      as.vector(xi.exb.qr[[q]][[r]])) * ti.rml.gq.w[[q]])}
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
                       b * c * d * as.vector(xi.exb.qr[[q]][[q]]),
                       ti.h0.q.l[[q]], ti.gq.PSI.qr[[q]][[r]], ti.S.gq.q.l[[q]],
                       ti.rml.gq.w.l[[q]], SIMPLIFY = FALSE))}
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
  theta.num.1 <- theta.den.1 <- gammascore.Z <- gammascore.1 <- list()
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
    gammascore.Z[[q]] <- Reduce("rbind", list(if(n.e[[q]]!=0) ze[[q]]),
                         if(n.i[[q]]!=0) zi[[q]])
    gammascore.1[[q]] <- rep(1, nrow(gammascore.Z[[q]]))
  }
  gammascore.mat <- gammascore <- gammahess.mat <- gammahess <- list()
  gammainc <- newgamma <- num.mat <- betascore.mat <- betascore <- list()
  betahess.mat <- betahess <- betainc <- dJ <- num.exb <- list()
  thetascore.num <- list()
  for(outer in 1:max.outer){
    for(iter in 1:max.iter){
      if(control$iter.disp==TRUE)
        cat("Outer iteration", outer, "of", max.outer, ": Inner iteration",
            iter, "of",max.iter, "\r")
      prev.oldgamm <- oldgamm
      prev.oldbeta <- oldbeta
      prev.oldtheta <- oldtheta
      # base <- updbase(oldtheta, tr.PSI, te.psi, te.PSI.qr, ti.gq.PSI.qr,
      #                 ti.rml.gq.w.psi)
      base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
               ti.gq.PSI.qr)
      llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                oldgamm, ze, zi)
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
                     oldlambda, SIMPLIFY = FALSE)
      llik0 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                  -unlist(llparms$pi.zi), na.rm = TRUE)
      for(r in 1:n.risk){
        base <- updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi, te.PSI.qr,
        ti.gq.PSI.qr)
        llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
        base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
        oldgamm, ze, zi)
        gammaparms <- updscoregamm(oldgamm, ze, zi)
        gammascore.mat[[r]] <- diag(c(if(n.e[[r]]!=0) gammaparms$ze.C.r[[r]],
        if(n.i[[r]]!=0) gammaparms$zi.C.r[[r]]))
        gammascore[[r]] <- t(gammascore.Z[[r]]) %*% gammascore.mat[[r]] %*%
          gammascore.1[[r]]
        gammahess.mat[[r]] <- diag(c(if(n.e[[r]]!=0) gammaparms$ze.ezg[[r]]^2 /
        (1 + gammaparms$ze.ezg[[r]]^3),
        if(n.i[[r]]!=0) gammaparms$zi.ezg[[r]]^2 /
        (1 + gammaparms$zi.ezg[[r]]^3)))
        gammahess[[r]] <- solve(t(gammascore.Z[[r]]) %*% gammahess.mat[[r]] %*%
        gammascore.Z[[r]])
        gammainc[[r]] <- gammahess[[r]] %*% gammascore[[r]]
        newgamm[[r]] = oldgamm[[r]] + as.vector(gammainc[[r]])
        llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
        base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
        newgamm, ze, zi)
        llik1 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                    -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                    -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                    -unlist(llparms$pi.zi), na.rm = TRUE)
        ome = 0.6
        while(llik1 <= llik0){
          newgamm[[r]] = oldgamm[[r]] + (ome * as.vector(gammainc[[r]]))
          llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                    base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                    newgamm, ze, zi)
          llik1 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                      -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                      -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                      -unlist(llparms$pi.zi), na.rm = TRUE)
          if(ome >= 1e-2) ome = ome * 0.6
          else if (ome < 1e-2 & ome >= 1e-5) ome = ome * 5e-2
          else if (ome < 1e-5 & ome > 1e-20) ome = ome * 1e-5
          else break
        }
        llik0 = llik1
        oldgamm = newgamm
      }
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
              oldlambda, SIMPLIFY = FALSE)
      llik1 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                  -unlist(llparms$pi.zi), na.rm = TRUE)
      for(r in 1:n.risk){
        base = updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi,
                       te.PSI.qr, ti.gq.PSI.qr)
        llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                  base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                  newgamm, ze, zi)
        betaparms = updscorebeta(llparms$ti.h.q, llparms$ti.S.gq.q,
                    base$ti.gq.H0.qr, llparms$xi.exb.qr, ti.rml.gq.w)
        betascore.mat[[r]] = diag(c(if(length(tr)!=0) -llparms$tr.H.q[[r]],
                             rep(1, n.e[[r]]), unlist(mapply(function(a,b)
                             if(a!=0) -b[[r]], n.e, llparms$te.H.qr, SIMPLIFY =
                             FALSE)), unlist(mapply(function(a,b,c)
                             if(a!=0) b[[r]] / c, n.i, betaparms$ti.A.qr,
                             llparms$ti.F.q, SIMPLIFY = FALSE))))
        betascore[[r]] = t(betascore.X[[r]]) %*% betascore.mat[[r]] %*%
                        betascore.1[[r]]
        betahess.mat[[r]] = diag(c(if(length(tr)!=0) llparms$tr.H.q[[r]],
                            unlist(mapply(function(a,b) if(a!=0) b[[r]], n.e,
                            llparms$te.H.qr, SIMPLIFY = FALSE)),
                            unlist(mapply(function(a,b,c) if(a!=0) (b[[r]] / c) ^ 2,
                            n.i, betaparms$ti.A.qr, llparms$ti.F.q,
                            SIMPLIFY = FALSE))))
        betahess[[r]] = solve(t(betahess.X[[r]]) %*% betahess.mat[[r]] %*%
                        betahess.X[[r]])
        betainc[[r]] = betahess[[r]] %*% betascore[[r]]
        newbeta[[r]] = oldbeta[[r]] + as.vector(betainc[[r]])
        print(oldbeta)
        llparms = updllparms(newbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                  base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                  oldgamm, ze, zi)
        llik2 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                    -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                    -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                    -unlist(llparms$pi.zi), na.rm = TRUE)
        ome = 0.6
        while(llik2 <= llik1){
          newbeta[[r]] = oldbeta[[r]] + (ome * as.vector(betainc[[r]]))
          llparms = updllparms(newbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                    base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                    newgamm, ze, zi)
          llik2 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                    -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                      -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                      -unlist(llparms$pi.zi), na.rm = TRUE)
          if(ome >= 1e-2) ome = ome * 0.6
          else if (ome < 1e-2 & ome >= 1e-5) ome = ome * 5e-2
          else if (ome < 1e-5 & ome > 1e-20) ome = ome * 1e-5
          else break
        }
        llik1 = llik2
        oldbeta = newbeta
      }
      R.pen = mapply(function(a,b,c) (t(a) %*% b %*% a)*c, oldtheta, R.mat,
              oldlambda, SIMPLIFY = FALSE)
      llik2 = sum(-unlist(llparms$tr.H.q), log(unlist(llparms$te.h.q)+1e-12),
                  -unlist(llparms$te.H.qr), log(unlist(llparms$ti.F.q)+1e-12),
                  -sum(unlist(R.pen)), -unlist(llparms$pi.ze),
                  -unlist(llparms$pi.zi), na.rm = TRUE)
      for(r in 1:n.risk){
        base = updbase(oldtheta, tr.PSI, te.psi, ti.rml.gq.w.psi,
               te.PSI.qr, ti.gq.PSI.qr)
        llparms = updllparms(oldbeta, xr, xe, base$tr.H0.q, xi, base$te.h0.q,
                  base$te.H0.qr, base$ti.h0.q, base$ti.gq.S0r.qr, ti.rml.gq.w,
                  oldgamm, ze, zi)
        thetaparms = updscoretheta(oldtheta, llparms$ti.S.gq.q, ti.gq.psi.qr,
                     ti.rml.gq.w.l, base$ti.h0.q.l, llparms$xi.exb.qr)
        num.mat[[r]] = Reduce("rbind", list(if(n.e[[r]]!=0) te.psi[[r]] /
                       as.vector(base$te.h0.q[[r]]), Reduce("rbind",
                       mapply(function(a,b,c) if(a!=0) b[[r]] / as.vector(c),
                       n.i, thetaparms$ti.B1.qr, llparms$ti.F.q,
                       SIMPLIFY = FALSE))))
        dJ[[r]] = 2*(oldtheta[[r]] %*% R.mat[[r]])
        num.exb[[r]] = Reduce("rbind", list(matrix(1, nrow = n.e[[r]],
                       ncol = 1), Reduce("rbind", mapply(function(a,b) if(a!=0)
                       {b[[r]]}, n.i, llparms$xi.exb.qr,
                       SIMPLIFY = FALSE))))
        thetascore.num[[r]] = (as.vector(t(theta.num.1[[r]]) %*%
                              (num.mat[[r]]*as.vector(num.exb[[r]])))) -
                              as.vector(oldlambda[[r]]) * pmin(dJ[[r]], 0) +
                              1e-12
      }
    }
  }
  print(llik0)
  print(llik1)
  print(llik2)
}





































