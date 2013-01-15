#' Goodman-Kruskal's \eqn{\gamma}
#'
#' Computes Goodman-Kruskal's \eqn{\gamma}
#' @param M a matrix
#' @return
#' \item{scon}{concordance}
#' \item{sdis}{disconcordance}
#' \item{gamma}{a real number between -1 and 1. calculated as
#' \eqn{\code{gamma} = \frac{\code{scon}-\code{sdis}}{\code{scon}+\code{sdis}}}{(scon-sdis)/(scon+sdis)}}
#' @references Goodman LA, Kruskal WH (1954) Measures of association
#' for cross classifications, Journal of the American Statistical
#' Association, 49, 732-764.
#' @export
#' @author Chun Li

GKGamma <- function(M) {
  nrow <- nrow(M)
  row <- row(M)
  ncol <- ncol(M)
  col <- col(M)

  calc_scon <- function(i, j, M, row, col) {
    M[i,j] * sum(M[row > i & col > j])
  }
  calc_sdis <- function(i, j, M, row, col) {
    if(j > 1L)
      M[i,j] * sum(M[row > i & col < j])
    else
      0
  }

  temp <- expand.grid(i=seq_len(nrow-1L), j=seq_len(ncol))
  
  scon <- sum(unlist(mapply(calc_scon, temp$i, temp$j,
                            MoreArgs=list(M=M,row=row,col=col),
                            SIMPLIFY=FALSE, USE.NAMES=FALSE)))

  sdis <- sum(unlist(mapply(calc_sdis, temp$i, temp$j,
                            MoreArgs=list(M=M, row=row, col=col),
                            SIMPLIFY=FALSE, USE.NAMES=FALSE)))

  return(list(scon=scon, sdis=sdis, gamma=(scon-sdis)/(scon+sdis)))
}
