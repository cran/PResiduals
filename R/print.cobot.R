#' cobot class print method
#' @param x cobot object
#' @param ... arguments passed to print.default
#' @keywords print
#' @export
#' @method print cobot
#' @author Charles Dupont

print.cobot <- function(x, ...) {
  invisible(print(matrix(c(x$T1, x$T2, x$T3,
                           sqrt(x$varT1), sqrt(x$varT2), sqrt(x$varT3),
                           x$pvalT1, x$pvalT2, x$pvalT3),
                         ncol=3,
                         dimnames=list(c("Gamma(Obs) - Gamma(Exp)", "Correlation of Residuals","Covariance of Residuals"),c("est", "stderr", "p"))), ...))
}
