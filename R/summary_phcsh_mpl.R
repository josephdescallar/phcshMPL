#' @export
summary_phcsh_mpl <- function(object, sand = FALSE){
  if(sand == FALSE){
    col.names = c("Estimate", "Std. Error", "z-value", "Pr(>|z|)",
                "Gradient")
  }
  else if(sand == TRUE){
    col.names = c("Estimate", "Std. Error (sandwich)", "z-value", "Pr(>|z|)",
                "Gradient")
  }
  var.names = dimnames(object$data$X)[[2]]
  risks = list()
  for(r in 1:object$n.risk){
    if(sand == FALSE){
      matxB = cbind(object$"beta"[[r]],object$seB[[r]], object$"beta"[[r]] /
              object$seB[[r]],2*(1-stats::pnorm(abs(object$"beta"[[r]] /
              object$seB[[r]]))), object$beta.gradient[[r]])
      matxT = cbind(object$"theta"[[r]], object$seT[[r]], object$"theta"[[r]] /
              object$seT[[r]],2*(1-stats::pnorm(abs(object$"theta"[[r]] /
              object$seT[[r]]))), object$theta.gradient[[r]])
    }
    else if(sand == TRUE){
      matxB = cbind(object$"beta"[[r]],object$seB.sand[[r]], object$"beta"[[r]]
              / object$seB.sand[[r]],2*(1-stats::pnorm(abs(object$"beta"[[r]] /
              object$seB.sand[[r]]))), object$beta.gradient[[r]])
      matxT = cbind(object$"theta"[[r]], object$seT.sand[[r]],
              object$"theta"[[r]] / object$seT.sand[[r]],2 *
              (1-stats::pnorm(abs(object$"theta"[[r]] / object$seT.sand[[r]]))),
              object$theta.gradient[[r]])
    }
    dimnames(matxB) = list(paste(" ", var.names, sep = ""),col.names)
    colnames(matxT) = col.names
    risks[[r]] = list(Beta = matxB, Theta = matxT)
  }
  out <- list(risks=risks)
  out
}
