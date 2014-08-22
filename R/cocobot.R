#### better to add to "pgumbel.R"
qgumbel <- function(p, loc = 0, scale = 1, lower.tail = TRUE)
{
  if(!lower.tail) p <- 1-p 
  q <- -log(-log(p))
  q <- loc + scale * q
  return(q)
}

getCI <- function(ts, var, fisher, ci=0.95){
  if(!fisher){
    lower <- ts - abs(qnorm(0.5*(1-ci)))*sqrt(var)
    upper <- ts + abs(qnorm(0.5*(1-ci)))*sqrt(var)  
  } else {
    ts_f <- log( (1+ts)/(1-ts) )
    var_f <- var*(2/(1-ts^2))^2
    lower_f <- ts_f - abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    upper_f <- ts_f + abs(qnorm(0.5*(1-ci)))*sqrt(var_f)
    lower <- (exp(lower_f)-1)/(1+exp(lower_f))
    upper <- (exp(upper_f)-1)/(1+exp(upper_f))
  }
  return(c(lower, upper))
}

### no need for sandwich package now
lm.scores = function(y, X){
  N = length(y)  
  mod = lm(y~X)
  smod = summary(mod)
  resid = smod$residuals
  ## bread = [1/N sum (- partial phi)]^-1 
  ##       = [- 1/N  d2l. dtheta. dtheta]^-1
  #d2l.dtheta.dtheta = - solve(bread(mod))*N
  d2l.dtheta.dtheta = -crossprod(cbind(1, X))
  #dl.dtheta = estfun(mod)
  dl.dtheta <- resid*cbind(1, X)
  presid = 2*pnorm((y - mod$fitted.values)/smod$sigma) -1
  dresid.dtheta = t(cbind(-1, -X))
  dpresid.dtheta = t(cbind(-2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma,
                           -2*dnorm((y - mod$fitted.values)/smod$sigma)/smod$sigma *
                             X))
  
  f.y<-density(resid)
  fy.ry <- NULL
  presid.k <- NULL
  for (i in 1:length(resid)){
    fy.ry[i] <- f.y$y[which(abs(f.y$x-resid[i])==min(abs(f.y$x-resid[i])))]
    presid.k[i] <- sum(resid<resid[i])/length(resid) - sum(resid>resid[i])/length(resid)
  }
  dpresid.dtheta.k <- t(cbind(-2*fy.ry,
                              -2*fy.ry*X))
  list(mod = mod, 
       dl.dtheta = dl.dtheta,
       d2l.dtheta.dtheta = d2l.dtheta.dtheta,
       resid = resid,
       dresid.dtheta = dresid.dtheta,
       presid = presid,
       presid.k= presid.k,
       dpresid.dtheta = dpresid.dtheta,
       dpresid.dtheta.k = dpresid.dtheta.k)
}


corTS = function(xresid, yresid,
                  xz.dl.dtheta, yz.dl.dtheta,
                  xz.d2l.dtheta.dtheta, yz.d2l.dtheta.dtheta,
                  dxresid.dthetax, dyresid.dthetay,fisher=FALSE){
  
  TS = cor(xresid, yresid)
  
  xresid2 = xresid^2
  yresid2 = yresid^2
  xbyyresid = xresid * yresid
  mean.xresid = mean(xresid)
  mean.yresid = mean(yresid)
  mean.xbyyresid = mean(xbyyresid)
  
  bigphi = cbind(xz.dl.dtheta,
                 yz.dl.dtheta,
                 mean.xresid - xresid,
                 mean.yresid - yresid,
                 mean.xbyyresid - xbyyresid,
                 mean(xresid2)-xresid2,
                 mean(yresid2)-yresid2,
                 0)
  
  npar.xz = dim(xz.dl.dtheta)[2]
  npar.yz = dim(yz.dl.dtheta)[2]
  Ntheta = npar.xz + npar.yz + 6
  N = dim(xz.dl.dtheta)[1]
  
  A = matrix(0,Ntheta,Ntheta)
  A[1:npar.xz, 1:npar.xz] = xz.d2l.dtheta.dtheta
  A[npar.xz+(1:npar.yz), npar.xz+(1:npar.yz)] = yz.d2l.dtheta.dtheta
  A[Ntheta-6+(1:6), Ntheta-6+(1:6)] = diag(N, 6)
  
  bigpartial = rbind(c(dxresid.dthetax %*% rep(1, N), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% rep(1, N)),
                     c(dxresid.dthetax %*% yresid, dyresid.dthetay %*% xresid),
                     c(dxresid.dthetax %*% (2*xresid), rep(0, npar.yz)),
                     c(rep(0, npar.xz), dyresid.dthetay %*% (2*yresid)))
  
  A[Ntheta-6+(1:5), 1:(npar.xz+npar.yz)] = -bigpartial
  
  ## TS also equals numTS / sqrt(varprod) = numTS * revsvp
  numTS = mean.xbyyresid - mean.xresid * mean.yresid
  var.xresid = mean(xresid2) - mean.xresid^2
  var.yresid = mean(yresid2) - mean.yresid^2
  varprod = var.xresid * var.yresid
  revsvp = 1/sqrt(varprod)
  dTS.dvarprod = numTS * (-0.5) * revsvp^3
  
  smallpartial = N *
    c(-mean.yresid * revsvp + dTS.dvarprod * (-2*mean.xresid*var.yresid),
      -mean.xresid * revsvp + dTS.dvarprod * (-2*mean.yresid*var.xresid),
      revsvp,
      dTS.dvarprod * var.yresid,
      dTS.dvarprod * var.xresid)
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  A[Ntheta, Ntheta-6+(1:5)] = -smallpartial
  
  SS = solve(A, t(bigphi))
  var.theta = tcrossprod (SS, SS)
  varTS = var.theta[Ntheta, Ntheta]
  pvalTS = 2 * pnorm( -abs(TS)/sqrt(varTS))
  
  if (fisher==TRUE){
    ####Fisher's transformation
    TS_f <- log( (1+TS)/(1-TS) )
    var.TS_f <- varTS*(2/(1-TS^2))^2
    pvalTS <- 2 * pnorm( -abs(TS_f)/sqrt(var.TS_f))
  }
  
  list(TS=TS,var.TS=varTS, pval.TS=pvalTS)
}


#'  Conditional continuous by ordinal tests for association.
#' 
#' \code{cocobot} tests for independence between an ordered categorical
#' variable, \var{X}, and a continuous variable, \var{Y}, conditional on other variables,
#' \var{Z}.  The basic approach involves fitting an ordinal model of \var{X} on
#' \var{Z}, a linear model of \var{Y} on \var{Z}, and then determining whether there is any
#' residual information between \var{X} and \var{Y}.  This is done by
#' computing residuals for both models, calculating their correlation, and 
#' testing the null of no residual correlation.  This procedure is analogous to test statistic 
#' \code{T2} in \code{cobot}.  Two test statistics (correlations) are currently output.  The first
#' is the correlation between probability-scale residuals (PResid). The second is the correlation between 
#' the observed-minus-expected residual for the continuous outcome model and a latent variable residual
#' for the ordinal model.
#' 
#' Formula is specified as \code{\var{X} | \var{Y} ~ \var{Z}}.
#' This indicates that models of \code{\var{X} ~ \var{Z}} and
#' \code{\var{Y} ~ \var{Z}} will be fit.  The null hypothsis to be
#' tested is \eqn{H_0 : X}{H0 : X} independant of \var{Y} conditional
#' on \var{Z}.  The ordinal variable, \code{\var{X}}, must precede the \code{|} and be a factor variable, and \code{\var{Y}} must be continuous.
#'  
#' @references Li C and Shepherd BE (2012) 
#' A new residual for ordinal outcomes.
#' \emph{Biometrika}. \bold{99}: 473--480.
#' @references Shepherd BE, Li C, Liu Q (submitted)
#' Probability-scale residuals for continuous, discrete, and censored data.
#'

#' @param formula an object of class \code{\link{Formula}} (or one
#' that can be coerced to that class): a symbolic description of the
#' model to be fitted.  The details of model specification are given
#' under \sQuote{Details}.
#' 
#' @param link The link family to be used for the ordinal model of 
#' \var{X} on \var{Z}.  Defaults to \samp{logit}. Other options are
#' \samp{probit}, \samp{cloglog}, and \samp{cauchit}.
#' 
#' @param data an optional data frame, list or environment (or object
#' coercible by \code{\link{as.data.frame}} to a data frame)
#' containing the variables in the model.  If not found in
#' \code{data}, the variables are taken from
#' \code{environment(formula)}, typically the environment from which
#' \code{cocobot} is called.
#' 
#' @param subset an optional vector specifying a subset of
#' observations to be used in the fitting process.
#' 
#' @param na.action action to take when \code{NA} present in data.
#' 
#' @param emp logical indicating whether the residuals from the model of
#' \var{Y} on \var{Z} are computed based on the assumption of normality (\code{FALSE}) 
#' or empirically (\code{TRUE}).
#' 
#' @param fisher logical indicating whether to apply fisher transformation to compute confidence intervals and p-values for the correlation.
#' 
#' @param conf.int numeric specifying confidence interval coverage.
#' 
#' @return object of \samp{cocobot} class.
#' @export
#' @examples
#' data(PResidData)
#' cocobot(y|w ~ z, data=PResidData)
#' @importFrom rms lrm
#' @importFrom sandwich bread estfun

cocobot <- function(formula, data, link=c("logit", "probit", "cloglog", "cauchit"),
                      subset, na.action=getOption('na.action'), 
                      emp=TRUE,fisher=FALSE,conf.int=0.95) {

  # Construct the model frames for x ~ z and y ~ z
  F1 <- Formula(formula)
  Fx <- formula(F1, lhs=1)
  Fy <- formula(F1, lhs=2)
  
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action", 
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- na.action
  mf[[1L]] <- as.name("model.frame")
  
  
  mx <- my <- mf

  # NOTE: we add the opposite variable to each model frame call so that
  # subsetting occurs correctly. Later we strip them off.
  mx[["formula"]] <- Fx
  yName <- all.vars(Fy[[2]])[1]
  mx[[yName]] <- Fy[[2]]

  my[["formula"]] <- Fy
  xName <- all.vars(Fx[[2]])[1]
  my[[xName]] <- Fx[[2]]

  mx <- eval(mx, parent.frame())
  mx[[paste('(',yName,')',sep='')]] <- NULL
  
  my <- eval(my, parent.frame())
  my[[paste('(',xName,')',sep='')]] <- NULL

  data.points <- nrow(mx)

  if (!is.factor(mx[[1]])){
    warning("Coercing ",names(mx)[1]," to factor. Check the ordering of categories.")
    mx[[1]] <- as.factor(mx[[1]])
  }

  if (is.factor(my[[1]])){
    stop(names(my)[1]," cannot be a factor.")
  }

  # Construct the model matrix z
  mxz <- model.matrix(attr(mx,'terms'),mx) 
  zzint <- match("(Intercept)", colnames(mxz), nomatch = 0L)
  if(zzint > 0L) {
    mxz <- mxz[, -zzint, drop = FALSE]
  }

  myz <- model.matrix(attr(my,'terms'),my) 
  zzint <- match("(Intercept)", colnames(myz), nomatch = 0L)
  if(zzint > 0L) {
    myz <- myz[, -zzint, drop = FALSE]
  }
  
  score.xz <- ordinal.scores(mx, mxz, method=link)
  score.yz <- lm.scores(y=model.response(my), X=myz)
  
  npar.xz = dim(score.xz$dl.dtheta)[2]
  npar.yz = dim(score.yz$dl.dtheta)[2]
  
  xx = as.integer(model.response(mx))

  nx = length(table(xx))

  N = length(xx)
  low.x = cbind(0, score.xz$Gamma)[cbind(1:N, xx)]
  hi.x = cbind(1-score.xz$Gamma, 0)[cbind(1:N, xx)]

  xz.presid <- low.x - hi.x
  xz.dpresid.dtheta <- score.xz$dlow.dtheta - score.xz$dhi.dtheta
  

  ## return value
  ans <- list(
          TS=list(),
          fisher=fisher,
          conf.int=conf.int,
          data.points=data.points
         )

  ### presid vs obs-exp
  #ta <-  corTS(xz.presid, score.yz$resid,
  #              score.xz$dl.dtheta, score.yz$dl.dtheta,
  #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
  #              xz.dpresid.dtheta, score.yz$dresid.dtheta, fisher)
  #ans$TS$TA <- 
  #    list( 
  #      ts=ta$TS, var=ta$var.TS, pval=ta$pval.TS,
  #      label='PResid vs. Obs-Exp'
  #    )
  
  if (emp==TRUE){
    ### presid vs presid (emprical)
    tb = corTS(xz.presid, score.yz$presid.k,
                score.xz$dl.dtheta, score.yz$dl.dtheta,
                score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                xz.dpresid.dtheta, score.yz$dpresid.dtheta.k,fisher)
    tb.label = "PResid vs. PResid (empirical)"
  } else {
    ### presid vs presid (use pdf of normal)
    tb = corTS(xz.presid, score.yz$presid,
                score.xz$dl.dtheta, score.yz$dl.dtheta,
                score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                xz.dpresid.dtheta, score.yz$dpresid.dtheta,fisher)
    tb.label = "PResid vs. PResid (assume normality)"
    
  }
  ans$TS$TB <- 
      list( 
        ts=tb$TS, var=tb$var.TS, pval=tb$pval.TS,
        label = tb.label
      )

  rij <- cbind(score.xz$Gamma, 1)[cbind(1:N, xx)]
  rij_1 <- cbind(0,score.xz$Gamma)[cbind(1:N, xx)]
  pij <- rij-rij_1

  G.inverse <- switch(link[1], logit = qlogis, probit = qnorm,
                 cloglog = qgumbel, cauchit = qcauchy)
  xz.latent.resid <- rep(NA, N)
  
  inverse_fail <- FALSE 
  for (i in 1:N){
    tmp <- try(integrate(G.inverse, rij_1[i], rij[i])$value/pij[i],silent=TRUE)
    if (inherits(tmp,'try-error')){
      if (link[1] != 'cauchit')
        warning("Cannot compute latent variable residual.")
      else
        warning("Cannot compute latent variable residual with link function cauchit.")
      inverse_fail <- TRUE
      break
    } else {
      xz.latent.resid[i] <- tmp
    }
  }
  
  if (!inverse_fail){
    ### To compute dlatent.dtheta (need dgamma.dtheta and dp0.dtheta from ordinal scores)
    xz.dlatent.dtheta = dpij.dtheta = matrix(, npar.xz, N)

    drij_1.dtheta <- score.xz$dlow.dtheta
    drij.dtheta <- -score.xz$dhi.dtheta
    for(i in 1:N) {  
      dpij.dtheta[,i] <- score.xz$dp0.dtheta[i, xx[i],]
      
      if (xx[i] == 1) {
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - 0 )
      } else if(xx[i] == nx){
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          0 - G.inverse(rij_1[i])*drij_1.dtheta[,i] )
      } else
        xz.dlatent.dtheta[,i] <- -xz.latent.resid[i]/pij[i]*dpij.dtheta[,i] + 1/pij[i]*(
          G.inverse(rij[i])*drij.dtheta[,i] - G.inverse(rij_1[i])*drij_1.dtheta[,i])
    }
    
    
    ### latent.resid vs obs-exp
    tc <- corTS(xz.latent.resid, score.yz$resid,
                 score.xz$dl.dtheta, score.yz$dl.dtheta,
                 score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
                 xz.dlatent.dtheta, score.yz$dresid.dtheta,fisher)
    ans$TS$TC <- 
        list( 
          ts=tc$TS, var=tc$var.TS, pval=tc$pval.TS,
          label = 'Latent.resid vs. Obs-Exp'
        )
    
    #if (emp==TRUE){
    #  ### latent vs presid (emprical)
    #  td = corTS(xz.latent.resid, score.yz$presid.k,
    #              score.xz$dl.dtheta, score.yz$dl.dtheta,
    #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
    #              xz.dlatent.dtheta, score.yz$dpresid.dtheta.k,fisher)
    #  td.label = "Latent.resid  vs. PResid (empirical)"
    #} else {
    #  ### latent vs presid (use pdf of normal)
    #  td = corTS(xz.latent.resid, score.yz$presid,
    #              score.xz$dl.dtheta, score.yz$dl.dtheta,
    #              score.xz$d2l.dtheta.dtheta, score.yz$d2l.dtheta.dtheta,
    #              xz.dlatent.dtheta, score.yz$dpresid.dtheta,fisher)
    #  td.label = "Latent.resid  vs. PResid (assume normality)"
    #  
    #}
    #ans$TS$TD <- 
    #    list( 
    #      ts=td$TS, var=td$var.TS, pval=td$pval.TS,
    #      label = td.label
    #    )
  }
  
  ans <- structure(ans, class="cocobot")

  # Apply confidence intervals
  for (i in seq_len(length(ans$TS))){
    ts_ci <- getCI(ans$TS[[i]]$ts,ans$TS[[i]]$var,ans$fisher,conf.int)
    ans$TS[[i]]$lower <- ts_ci[1]
    ans$TS[[i]]$upper <- ts_ci[2]
  }

  ans
  
}
