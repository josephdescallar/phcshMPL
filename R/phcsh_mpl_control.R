#' Ancillary arguments for controling coxph_mpl fits
#'
#' This is used to set various numeric parameters controling a Cox model fit
#' using phcsh_mpl. Typically it would only be used in a call to phcsh_mpl. Some
#' basic checks are performed on inputs, such that impossible argument values
#' (like a negative number of events per base, for example) are avoided.
#'
#' @param knots.perc.limit test
#' @param gq.points test
#' @param dgr test
#' @param basis.intercept test
#' @param max.outer test
#' @param max.iter test
#' @param lambda test
#' @param n.basis test
#' @param iter.disp test
#' @param iknots.pos test
#' @param aps test
#' @param tmid test
#' @param inner.conv test
#' @export
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
