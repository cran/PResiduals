###glm()
#' @export
presid.glm <- function(object, emp=FALSE, ...) {
    ##gaussian.emp = (2 * rank(residuals(object))-1-length(y))/length(y)  #need to figure out how to incorporate this in this function
    switch(object$family$family ,
           poisson = 2 * ppois(object$y, object$fitted.values) - dpois(object$y, object$fitted.values) - 1,
           binomial = object$y - object$fitted.values,
           gaussian = if(emp) (2 * rank(residuals(object)) - 1 - length(object$y)) / length(object$y) else 2 * pnorm((object$y - object$fitted.values)/sqrt(summary(object)$dispersion)) - 1,
           stop("Unhandled family", object$family$family))
}

###lm()
#' @export
presid.lm <- function(object, emp=FALSE, ...) {
    y <- model.response(object$model)
  if(emp) {
      (2 * rank(residuals(object)) - 1 - length(y)) / length(y)
  } else {
      2 * pnorm((y - object$fitted.values)/summary(object)$sigma) - 1
  }
}

###negative binomial
#' @export
presid.negbin <- function(object, ...) {
  pnbinom(object$y-1, mu=object$fitted.values, size=object$theta ) + pnbinom(object$y, mu=object$fitted.values, size=object$theta) -1
}


###polr
#' @export
presid.polr <- function(object, ...) {
    pfun <- switch(object$method,
                   logistic = plogis,
                   probit = pnorm, 
                   loglog = pgumbel,
                   cloglog = pGumbel,
                   cauchit = pcauchy)
    n <- length(object$lp)
    q <- length(object$zeta)
    cumpr <- cbind(0, matrix(pfun(matrix(object$zeta, n, q, byrow = TRUE) - object$lp),, q), 1)
    y <- as.integer(model.response(object$model))
    lo <- cumpr[cbind(seq_len(n), y)]
    hi <- 1 - cumpr[cbind(seq_len(n), y+1L)]
    lo - hi
}

###coxph()
#' @export
presid.coxph <- function(object, ...) {
    time <- object$y[,1]
    delta <- object$y[,2]
    
    1 - exp(residuals(object) - delta) - delta*exp(residuals(object) - delta)
}


###survreg()
#' @export
presid.survreg <- function(object, ...){
    time <- object$y[,1]
    delta <- object$y[,2]
    
    switch(object$dist,
           weibull = {
               prob <- pweibull(exp(time), shape=1/summary(object)$scale,
                                scale=exp(object$linear.predictors),
                                lower.tail=TRUE, log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           exponential = {
               prob <- pexp(time, rate=1/exp(object$linear.predictors),
                            lower.tail=TRUE, log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           gaussian = {
               prob <- pnorm(time, mean=object$linear.predictors,
                             sd=summary(object)$scale, lower.tail=TRUE,
                             log.p=FALSE)
               prob + delta*(prob - 1)
           },
           
           logistic = {
               prob <- plogis(time, location=object$linear.predictors,
                              scale=summary(object)$scale, lower.tail=TRUE,
                              log.p=FALSE)
               prob + delta*(prob - 1)
           },
         
           ## loglogistic = {
           ##     prob <- pllogis(time, shape=summary(object)$scale,
           ##                     scale=exp(object$linear.predictors),
           ##                     lower.tail=TRUE, log.p=FALSE)
           ##     prob + delta*(prob - 1)
           ## },
         
           lognormal = {
               prob <- plnorm(time, meanlog=object$linear.predictors,
                              sdlog=summary(object)$scale, lower.tail=TRUE,
                              log.p=FALSE)
               prob + delta*(prob - 1)
           },
           stop("Unhandled dist", object$dist))
}

#' @export
presid.default <- function(object, ...) {
    stop("Unhandled model type")
}


#' Probability-scale Residual
#'
#' \code{presid} Calculates the probability-scale residual for various model
#' function objects. Currently supported models include \code{\link{glm}}
#' (Poisson, binomial, and gaussian families), \code{\link{lm}} in the
#' \pkg{stats} library, \code{\link{survreg}} (Weibull, exponential, gaussian,
#' logistic, and lognormal distributions) and \code{\link{coxph}} in the
#' \pkg{survival} library, and \code{\link{polr}} and \code{\link{glm.nb}} in
#' the \pkg{MASS} library.
#' 
#' Probability-scale residual is \eqn{P(Y* < y) - P(Y* > y)} where \eqn{y} is the observed
#' outcome and \eqn{Y*} is a random variable from the fitted distribution.
#'
#' @param object The model object for which the probability-scale residual is calculated
#' @param ... Additional arguements passed to methods
#' @return The probability scale residual for the model
#' @references Shepherd BE, Li C, Liu Q.  Probability-scale residuals for continuous,
#' discrete, and censored data.  Submitted.
#' @references Li C and Shepherd BE, A new residual for ordinal
#' outcomes. Biometrika 2012; 99:473-480
#' @author Charles Dupont \email{charles.dupont@@vanderbilt.edu}
#' @author Chun Li \email{chun.li@@vanderbilt.edu}
#' @author Bryan Shepherd \email{bryan.shepherd@@vanderbilt.edu}
#' @importFrom actuar pllogis
#' @importFrom stats plnorm pnorm pexp pweibull plogis pnbinom pcauchy
#' @export
#' @examples
#' library(survival)
#' library(stats)
#' 
#' set.seed(100)
#' n <- 1000
#' x <- rnorm(n)
#' t <- rweibull(n, shape=1/3, scale=exp(x))
#' c <- rexp(n, 1/3)
#' y <- pmin(t, c)
#' d <- ifelse(t<=c, 1, 0)
#'
#' mod.survreg <- survreg(Surv(y, d) ~ x, dist="weibull")
#' summary(presid(mod.survreg))
#' plot(x, presid(mod.survreg))
#' 
#' ##### example for proprotional hazards model
#' n <- 1000
#' x <- rnorm(n)
#' beta0 <- 1
#' beta1 <- 0.5
#' t <- rexp(n, rate = exp(beta0 + beta1*x))
#' c <- rexp(n, rate=1)
#' y <- ifelse(t<c, t, c)
#' delta <- as.integer(t<c)
#' 
#' mod.coxph <- coxph(Surv(y, delta) ~ x)
#' presid <- presid(mod.coxph)
#' plot(x, presid, cex=0.4, col=delta+2)
#'
#' #### example for Negative Binomial regression
#' library(MASS)
#' 
#' n <- 1000
#' beta0 <- 1
#' beta1 <- 0.5
#' x <- runif(n, min=-3, max=3)
#' y <- rnbinom(n, mu=exp(beta0 + beta1*x), size=3)
#' 
#' mod.glm.nb <- glm.nb(y~x)
#' presid <- presid(mod.glm.nb)
#' summary(presid)
#' plot(x, presid, cex=0.4)
#' 
#' ##### example for proportional odds model
#' library(MASS)
#' 
#' n <- 1000
#' x  <- rnorm(n)
#' y  <- numeric(n)
#' alpha = c(-1, 0, 1, 2)
#' beta <- 1
#' py  <-  (1 + exp(- outer(alpha, beta*x, "+"))) ^ (-1)
#' aa = runif(n)
#' for(i in 1:n)
#'   y[i] = sum(aa[i] > py[,i])
#' y <-  as.factor(y)
#' 
#' 
#' mod.polr <- polr(y~x, method="logistic")
#' summary(mod.polr)
#' presid <- presid(mod.polr)
#' summary(presid)
#' plot(x, presid, cex=0.4)
presid <- function(object, ...) {
    UseMethod('presid', object)
}
