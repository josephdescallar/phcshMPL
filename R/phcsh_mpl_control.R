phcsh_mpl_control <- function(knots.perc.limit = c(0.05,0.95), gq.points = 15,
                     dgr = 2, basis.intercept = FALSE, max.outer = 100L,
                     max.iter = 10000L, lambda, n.basis = NULL,
                     iter.disp = TRUE, iknots.pos = NULL, aps = TRUE,
                     tmid = TRUE, inner.conv = 1e-8){out = list(
                     knots.perc.limit = knots.perc.limit, gq.points = gq.points,
                     dgr=dgr, basis.intercept = basis.intercept,
                     max.outer = max.outer, max.iter = max.iter,
                     lambda = lambda, n.basis = n.basis, iter.disp = iter.disp,
                     iknots.pos = iknots.pos, aps = aps, tmid = tmid,
                     inner.conv = inner.conv)
  out
}
