ordinal.scores <- function(mf, method) {
# This function is called from COBOT.scores().
# It calculates all values needed for estimating equations.
  ## mf is the model.frame of the data

  ## Fit proportional odds model and obtain the MLEs of parameters.
  mod <- newpolr(mf, Hess=TRUE, method=method)
  y <- model.response(mf)
  
  ## N: number of subjects; ny: number of y categories
  N = length(y)
  ny = length(mod$lev)

  ## na, nb: number of parameters in alpha and beta
  na = ny - 1

  Terms <- attr(mf, "terms")
  x <- model.matrix(Terms, mf, contrasts)
  xint <- match("(Intercept)", colnames(x), nomatch = 0L)

  nb = ncol(x)

  if(xint > 0L) {
    x <- x[, -xint, drop = FALSE]
    nb <- nb - 1L
  }
  
  npar = na + nb

  alpha = mod$zeta
  beta = -coef(mod)
  eta <- mod$lp

  ## Predicted probabilities p0 and dp0.dtheta are stored for individuals.
  p0 = mod$fitted.values
  dp0.dtheta = array(dim=c(N, ny, npar))
  Y <- matrix(nrow=N,ncol=ny)
  .Ycol <- col(Y)

  edcumpr <- cbind(mod$dcumpr, 0)
  e1dcumpr <- cbind(0,mod$dcumpr)

  for(i in seq_len(na)) {
    dp0.dtheta[,,nb+i] <- (.Ycol == i) * edcumpr - (.Ycol == i + 1)*e1dcumpr
  }

  for(i in seq_len(nb)) {
    dp0.dtheta[,,i] <- -mod$dfitted.values * x[,i]
  }
  
  ## Cumulative probabilities
  Gamma = mod$cumpr
  dcumpr <- mod$dcumpr
  
  ## Scores are stored for individuals.
  dl.dtheta = mod$grad

  d2l.dtheta.dtheta <- mod$Hessian
  
  
  list(mod = mod,
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       p0 = p0,
       dp0.dtheta = dp0.dtheta,
       Gamma = Gamma,
       dcumpr=dcumpr)
}

#' Conditional ordinal by ordinal tests for association.
#'
#' \code{cobot} tests for independence between two ordered categorical
#' variables, \var{X} and \var{Y} conditional on other variables,
#' \var{Z}.  The basic approach involves fitting models of \var{X} on
#' \var{Z} and \var{Y} on \var{Z} and determining whether there is any
#' remaining information between \var{X} and \var{Y}.  This is done by
#' computing one of 3 test statistics.  \code{T1} compares empirical
#' distribution of \var{X} and \var{Y} with the joint fitted
#' distribution of \var{X} and \var{Y} under independence conditional
#' on \var{Z}. \code{T2} computes the correlation between ordinal
#' (probability scale) residuals from both models and tests the null
#' of no residual correlation.  \code{T3} evaluates the
#' concordance--disconcordance of data drawn from the joint fitted
#' distribution of \var{X} and \var{Y} under conditional independence
#' with the empirical distribution. Details are given in \cite{Li C and
#' Shepherd BE, Test of association between two ordinal variables
#' while adjusting for covariates. Journal of the American Statistical
#' Association 2010, 105:612-620}.
#'
#' formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothsis to be
#' tested is \eqn{H_0 : X}{H0 : X} independant of \var{Y} conditional
#' on \var{Z}.
#' 
#' Note that \code{T2} can be thought of as an adjust rank
#' correlation.(\cite{Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480})
#'
#' @param formula an object of class \code{\link{Formula}} (or one
#'   that can be coerced to that class): a symbolic description of the
#'   model to be fitted.  The details of model specification are given
#'   under \sQuote{Details}.
#' @param link The link family to be used for ordinal models of both
#' \var{X} and \var{Y}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}.
#' @param link.x The link function to be used for a model of the first
#' ordered variable. Defaults to value of \code{link}.
#' @param link.y The link function to be used for a model of the
#' second variable. Defaults to value of \code{link}.
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{cobot} is called.
#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' @return object of \samp{cobot} class.
#' @references Li C and Shepherd BE, Test of association between two
#' ordinal variables while adjusting for covariates. Journal of the
#' American Statistical Association 2010, 105:612-620.
#' @references Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480
#' @import Formula
#' @export
#' @seealso \code{\link{Formula}}, \code{\link{as.data.frame}}
#' @author Charles Dupont \email{charles.dupont@@vanderbilt.edu}
#' @author Chun Li \email{chun.li@@vanderbilt.edu}
#' @author Bryan Shepherd \email{bryan.shepherd@@vanderbilt.edu}
#' @include newPolr.R
#' @include diagn.R
#' @include GKGamma.R
#' @include pgumbel.R
#' @examples
#' data(PResidData)
#' cobot(x|y~z, data=PResidData)
cobot <- function(formula, link=c("logit", "probit", "cloglog", "cauchit"),
                  link.x=link,
                  link.y=link,
                  data, subset) {
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.fail
  mf[[1L]] <- as.name("model.frame")


  mx <- my <- mf

  mx[["formula"]] <- Fx
  my[["formula"]] <- Fy

  mx <- eval(mx, parent.frame())
  my <- eval(my, parent.frame())

  score.xz <- ordinal.scores(mx, method=link.x)
  score.yz <- ordinal.scores(my, method=link.y)

  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]

  Terms <- attr(mx, "terms")
  zz <- model.matrix(Terms, mx, contrasts)
  zzint <- match("(Intercept)", colnames(zz), nomatch = 0L)

  

  if(zzint > 0L) {
    zz <- zz[, -zzint, drop = FALSE]
  }

  xx = as.integer(model.response(mx)); yy = as.integer(model.response(my))
  nx = length(table(xx))
  ny = length(table(yy))
  N = length(yy)


#### Asymptotics for T3 = mean(Cprob) - mean(Dprob)
  ##  If gamma.x[0]=0 and gamma.x[nx]=1, then
  ##  low.x = gamma.x[x-1], hi.x = (1-gamma.x[x])
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  low.y = cbind(0, score.yz$Gamma)[cbind(1:N, yy)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]
  hi.y = cbind(1-score.yz$Gamma, 0)[cbind(1:N, yy)]
  Cprob = low.x*low.y + hi.x*hi.y
  Dprob = low.x*hi.y + hi.x*low.y
  mean.Cprob = mean(Cprob)
  mean.Dprob = mean(Dprob)
  T3 = mean.Cprob - mean.Dprob

  dcumpr.xz <- cbind(0, score.xz$dcumpr, 0)
  dcumpr.yz <- cbind(0, score.yz$dcumpr, 0)
  dlowx.dthetax <- cbind((col(score.xz$dcumpr) == (xx - 1L)) * score.xz$dcumpr,
                         dcumpr.xz[cbind(1:N,xx)] * zz)
  dlowy.dthetay <- cbind((col(score.yz$dcumpr) == (yy - 1L)) * score.yz$dcumpr,
                         dcumpr.yz[cbind(1:N,yy)] * zz)

  dhix.dthetax <- -cbind((col(score.xz$dcumpr) == xx) * score.xz$dcumpr,
                         dcumpr.xz[cbind(1:N,xx + 1L)] * zz)
  dhiy.dthetay <- -cbind((col(score.yz$dcumpr) == yy) * score.yz$dcumpr,
                         dcumpr.yz[cbind(1:N,yy + 1L)] * zz)

  dCsum.dthetax <- crossprod(dlowx.dthetax, low.y) +
    crossprod(dhix.dthetax, hi.y)
  dCsum.dthetay <- crossprod(dlowy.dthetay, low.x) +
    crossprod(dhiy.dthetay, hi.x)
  dDsum.dthetax <- crossprod(dlowx.dthetax, hi.y) +
    crossprod(dhix.dthetax, low.y)
  dDsum.dthetay <- crossprod(dlowy.dthetay, hi.x) +
    crossprod(dhiy.dthetay, low.x)

  dT3sum.dtheta = c(dCsum.dthetax-dDsum.dthetax, dCsum.dthetay-dDsum.dthetay)

  ## Estimating equations for (theta, p3)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## p3 is the "true" value of test statistic, and the equation is
  ## p3 - (Ci-Di)
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, T3-(Cprob-Dprob))

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta, -Ntheta] = -dT3sum.dtheta
  A[Ntheta, Ntheta] = N

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT3 = var.theta[Ntheta, Ntheta]

  pvalT3 = 2 * pnorm(-abs(T3)/sqrt(varT3))


#### Asymptotics for T4 = (mean(Cprob)-mean(Dprob))/(mean(Cprob)+mean(Dprob))
  T4 = (mean.Cprob - mean.Dprob)/(mean.Cprob + mean.Dprob)

  ## Estimating equations for (theta, P4)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P4 is a vector of (cc, dd, p4).  Their corresponding equations are:
  ## cc - Ci
  ## dd - Di
  ## p4 - (cc-dd)/(cc+dd)
  ## Then p4 is the "true" value of test statistic.
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.Cprob - Cprob, mean.Dprob - Dprob, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 3
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-3+(1:3), Ntheta-3+(1:3)] = diag(N, 3)
  A[Ntheta-2, 1:(npar.xz+npar.yz)] = -c(dCsum.dthetax, dCsum.dthetay)
  A[Ntheta-1, 1:(npar.xz+npar.yz)] = -c(dDsum.dthetax, dDsum.dthetay)

  revcpd = 1/(mean.Cprob + mean.Dprob)
  dT4.dcpd = (mean.Cprob-mean.Dprob)*(-revcpd^2)
  A[Ntheta, Ntheta-3+(1:2)] = -N * c(revcpd+dT4.dcpd, -revcpd+dT4.dcpd)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT4 = var.theta[Ntheta, Ntheta]

  pvalT4 = 2 * pnorm(-abs(T4)/sqrt(varT4))


#### Asymptotics for T2 = cor(hi.x - low.x, hi.y - low.y)
  xresid = hi.x - low.x
  yresid = hi.y - low.y
  T2 = cor(xresid, yresid)

  xbyyresid = xresid * yresid
  xresid2 = xresid^2
  yresid2 = yresid^2
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  ## T2 also equals numT2 / sqrt(varprod) = numT2 * revsvp
  numT2 = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dT2.dvarprod = numT2 * (-0.5) * revsvp^3

  ## Estimating equations for (theta, P5)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## P5 is a vector (ex, ey, crossxy, ex2, ey2, p5).
  ## Their corresponding equations are:
  ## ex - (hi.x-low.x)[i]
  ## ey - (hi.y-low.y)[i]
  ## crossxy - ((hi.x-low.x)*(hi.y-low.y))[i]
  ## ex2 - (hi.x-low.x)[i]^2
  ## ey2 - (hi.y-low.y)[i]^2
  ## p5 - (crossxy-ex*ey)/sqrt((ex2-ex^2)*(ey2-ey^2))
  ## Then p5 is the "true" value of test statistic
  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta,
    mean.xresid - xresid, mean.yresid - yresid, mean.xbyyresid - xbyyresid,
    mean(xresid2) - xresid2, mean(yresid2) - yresid2, 0)

  ## sandwich variance estimate of var(thetahat)
  Ntheta = npar.xz + npar.yz + 6
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)

  dxresid.dthetax = dhix.dthetax - dlowx.dthetax
  dyresid.dthetay = dhiy.dthetay - dlowy.dthetay
  bigpartial = rbind(c(crossprod(dxresid.dthetax, rep(1, N)), rep(0, npar.yz)),
    c(rep(0, npar.xz), crossprod(dyresid.dthetay, rep(1, N))),
    c(crossprod(dxresid.dthetax, yresid), crossprod(dyresid.dthetay, xresid)),
    c(crossprod(dxresid.dthetax, (2*xresid)), rep(0, npar.yz)),
    c(rep(0, npar.xz), crossprod(dyresid.dthetay, (2*yresid))))
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial

  smallpartial = N *
    c(-mean.yresid * revsvp + dT2.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dT2.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dT2.dvarprod * var.yresid,
      dT2.dvarprod * var.xresid)
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod(SS, SS)
  varT2 = var.theta[Ntheta, Ntheta]

  pvalT2 = 2 * pnorm(-abs(T2)/sqrt(varT2))


#### Asymptotics for T1 = tau - tau0
  ## dtau0/dtheta
  ## P0 is the sum of product predicted probability matrix with dim(nx,ny)
  P0 = crossprod(score.xz$p0, score.yz$p0) / N

  cdtau0 = GKGamma(P0)
  C0 = cdtau0$scon
  D0 = cdtau0$sdis
  ## C0 = sum_{l>j,m>k} {P0[j,k] * P0[l,m]}
  ## D0 = sum_{l>j,m<k} {P0[j,k] * P0[l,m]}
  dC0.dP0 = matrix(,nx,ny)
  dD0.dP0 = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC0.dP0[i,j] = ifelse(i>1 & j>1, sum(P0[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(P0[(i+1):nx, (j+1):ny]), 0)
      dD0.dP0[i,j] = ifelse(i>1 & j<ny, sum(P0[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(P0[(i+1):nx, 1:(j-1)]), 0)
    }

  ## tau0 = (C0-D0)/(C0+D0)
  dtau0.dC0 = 2*D0/(C0+D0)^2
  dtau0.dD0 =-2*C0/(C0+D0)^2

  ## ## P0 is already a matrix
  dP0.dtheta.x = array(0, c(nx, ny, npar.xz))
  for(j in 1:ny) {
    aa = matrix(0, nx, npar.xz)
    for(i in 1:N)
      aa = aa + score.xz$dp0.dtheta[i,,] * score.yz$p0[i,j]
    dP0.dtheta.x[,j,] = aa/N
    ## simpler but mind-tickling version
    #dP0.dtheta.x[,j,] = (score.yz$p0[,j] %*% matrix(score.xz$dp0.dtheta,N))/N
  }

  dP0.dtheta.y = array(0, c(nx, ny, npar.yz))
  for(j in 1:nx) {
    aa = matrix(0, ny, npar.yz)
    for(i in 1:N)
      aa = aa + score.yz$dp0.dtheta[i,,] * score.xz$p0[i,j]
    dP0.dtheta.y[j,,] = aa/N
  }

  ## dC0.dtheta and dD0.dtheta
  dC0.dtheta.x = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dD0.dtheta.x = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.x, nx*ny)
  dC0.dtheta.y = as.numeric(dC0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)
  dD0.dtheta.y = as.numeric(dD0.dP0) %*% matrix(dP0.dtheta.y, nx*ny)

  ## dtau0/dtheta
  dtau0.dtheta.x = dtau0.dC0 * dC0.dtheta.x + dtau0.dD0 * dD0.dtheta.x
  dtau0.dtheta.y = dtau0.dC0 * dC0.dtheta.y + dtau0.dD0 * dD0.dtheta.y


  ## dtau/dPa
  ## tau = (C-D)/(C+D)
  Pa = table(xx, yy) / N
  cdtau = GKGamma(Pa)
  C = cdtau$scon
  D = cdtau$sdis
  dtau.dC = 2*D/(C+D)^2
  dtau.dD =-2*C/(C+D)^2

  ## Pa[nx,ny] is not a parameter, but = 1 - all other Pa parameters.
  ## Thus, d.Pa[nx,ny]/d.Pa[i,j] = -1.
  ## Also, d.sum(Pa[-nx,-ny]).dPa[i,j] = 1 when i<nx and j<ny, and 0 otherwise.
  ##
  ## In C = sum_{l>j,m>k} {Pa[j,k] * Pa[l,m]}, Pa[i,j] appears in
  ##   Pa[i,j] * XX (minus Pa[nx,ny] if i<nx & j<ny), and in
  ##   sum(Pa[-nx,-ny]) * Pa[nx,ny].
  ## So, dC.dPa[i,j] = XX (minus Pa[nx,ny] if i<nx & j<ny)
  ##                  + d.sum(Pa[-nx,-ny]).dPa[i,j] * Pa[nx,ny]
  ##                  - sum(Pa[-nx,-ny])
  ##                 = XX (with Pa[nx,ny] if present) - sum(Pa[-nx,-ny])
  ##
  ## D = sum_{l>j,m<k} {Pa[j,k] * Pa[l,m]} doesn't contain Pa[nx,ny]
  dC.dPa = matrix(,nx,ny)
  dD.dPa = matrix(,nx,ny)
  for(i in 1:nx)
    for(j in 1:ny) {
      dC.dPa[i,j] = ifelse(i>1 & j>1, sum(Pa[1:(i-1), 1:(j-1)]), 0) +
        ifelse(i<nx & j<ny, sum(Pa[(i+1):nx, (j+1):ny]), 0) - sum(Pa[-nx,-ny])
      dD.dPa[i,j] = ifelse(i>1 & j<ny, sum(Pa[1:(i-1), (j+1):ny]), 0) +
        ifelse(i<nx & j>1, sum(Pa[(i+1):nx, 1:(j-1)]), 0)
    }

  dtau.dPa = dtau.dC * dC.dPa + dtau.dD * dD.dPa
  dtau.dPa = dtau.dPa[-length(dtau.dPa)]  ## remove the last value


  ## Estimating equations for (theta, phi)
  ## theta is (theta.xz, theta.yz) and the equations are score functions.
  ## phi is (p_ij) for (X,Y), and the equations are
  ## I{subject in cell (ij)} - p_ij
  phi.Pa = matrix(0, N, nx*ny)
  phi.Pa[cbind(1:N, xx+(yy-1)*nx)] = 1
  phi.Pa = phi.Pa - rep(1,N) %o% as.numeric(Pa)
  phi.Pa = phi.Pa[,-(nx*ny)]

  bigphi = cbind(score.xz$dl.dtheta, score.yz$dl.dtheta, phi.Pa)

  ## sandwich variance estimate of var(thetahat, phihat)
  Ntheta = npar.xz + npar.yz + nx*ny-1
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = score.xz$d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = score.yz$d2l.dtheta.dtheta
  A[-(1:(npar.xz+npar.yz)), -(1:(npar.xz+npar.yz))] = -diag(N, nx*ny-1)

  ## One way of coding:
  ##B = t(bigphi) %*% bigphi
  ##var.theta = solve(A) %*% B %*% solve(A)
  ## Suggested coding for better efficiency and accuracy:
  ##SS = solve(A, t(bigphi))
  ##var.theta = SS %*% t(SS)
  ## Or better yet, no need to explicitly obtain var.theta.  See below.

  ## Test statistic T1 = tau - tau0
  T1 = cdtau$gamma - cdtau0$gamma
  ## dT.dtheta has length nx + ny + nx*ny-1
  dT1.dtheta = c(-dtau0.dtheta.x, -dtau0.dtheta.y, dtau.dPa)

  ## variance of T, using delta method
  ##varT = t(dT.dtheta) %*% var.theta %*% dT.dtheta
  SS = crossprod(dT1.dtheta, solve(A, t(bigphi)))
  varT1 = sum(SS^2)
  pvalT1 = 2 * pnorm(-abs(T1)/sqrt(varT1))

  structure(list(T1=T1, varT1=varT1, pvalT1=pvalT1, T2=T2, varT2=varT2, pvalT2=pvalT2, T3=T3, varT3=varT3, pvalT3=pvalT3), class="cobot")
}
