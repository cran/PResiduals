\name{GKGamma}
\alias{GKGamma}
\title{Goodman-Kruskal's \eqn{\gamma}}
\usage{
  GKGamma(M)
}
\arguments{
  \item{M}{a matrix}
}
\value{
  \item{scon}{concordance} \item{sdis}{disconcordance}
  \item{gamma}{a real number between -1 and 1. calculated
  as \eqn{\code{gamma} =
  \frac{\code{scon}-\code{sdis}}{\code{scon}+\code{sdis}}}{(scon-sdis)/(scon+sdis)}}
}
\description{
  Computes Goodman-Kruskal's \eqn{\gamma}
}
\author{
  Chun Li
}
\references{
  Goodman LA, Kruskal WH (1954) Measures of association for
  cross classifications, Journal of the American
  Statistical Association, 49, 732-764.
}
