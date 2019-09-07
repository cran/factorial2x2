
#' Power of the 1/3-1/3-1/3 procedure
#'
#' Computes the power of the 1/3-1/3-1/3 procedure, that is, the power to
#' detect the overall A effect, the simple A effect, or the simple AB effect.
#'
#' @param  n  total subjects with n/4 subjects in each of the C, A, B, and AB groups
#' @param  hrA  group A to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrB  group B to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrAB  group AB to group C hazard ratio; \code{hrAB} < 1 corresponds to group AB superiority
#' @param  avgprob  event probability averaged across the C, A, B, and AB groups
#' @param  probA_C  event probability averaged across the A and C groups
#' @param  probAB_C  event probability averaged across the AB and C groups
#' @param  crit13  rejection critical value for the overall A, simple A, and simple AB logrank statistics
#' @param  dig number of decimal places to \code{\link{roundDown}} the critical value to
#' @param  cormat12  asymptotic correlation matrix for the overall A and simple A, respectively, simple AB logrank statistics
#' @param  cormat23  asymptotic correlation matrix for the simple A and  simple AB logrank statistics
#' @param  cormat123  asymptotic correlation matrix for the overall A, simple A, and  simple AB logrank statistics
#' @param  niter  number of times we call \code{pmvnorm} to average out its randomness
#' @param  abseps  \code{abseps} setting in the \code{pmvnorm} call
#' @return \item{poweroverA }{power to detect the overall A effect}
#' @return \item{powerA }{power to detect the simple A effect}
#' @return \item{powerAB }{power to detect the simple AB effect}
#' @return \item{power13.13.13 }{power to detect the overall A, simple A,  or simple AB effects, i.e.,
#' power of the 1/3-1/3-1/3 procedure}
#' @import mvtnorm
#' @details   For a 2-by-2 factorial design, this function computes
#' the probability that either the overall A
#' or the simple A or the simple AB logrank statistics
#' reject their null hypotheses at the
#' \code{crit13} critical value.  As described in Leifer, Troendle, et al. (2019),
#' the \code{crit13} = -2.32 critical value
#' corresponds to controlling the famiywise error of the 1/3-1/3-1/3 procedure at the
#' two-sided 0.05 significance level.
#' The critical value -2.32 may be computed using the \code{crit2x2} function.
#' The \code{pmvnorm} function
#' from the \code{mvtnorm} package is used to calculate
#' the powers for rejecting the pairwise and three-way intersections of
#' Since these powers involve bivariate, respectively, trivariate,
#'  normal integration over an unbounded region in R^2, respectively, R^3, \code{pmvnorm}
#' uses a random seed for these computations.  To smooth out the
#' randomness, \code{pmvnorm} is called \code{niter} times and
#' the average value over the \code{niter} calls is taken to be those powers.
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @references Slud, E.V. Analysis of factorial survival experiments. Biometrics. 1994; 50: 25-38.
#' @export power13_13_13
#' @seealso \code{\link{crit2x2}}, \code{lgrkPower}, \code{strLgrkPower}, \code{pmvnorm}
#' @examples
#' # Corresponds to scenario 5 in Table 2 from Leifer, Troendle, et al. (2019).
#' rateC <- 0.0445
#' hrA <- 0.80
#' hrB <- 0.80
#' hrAB <- 0.72
#' mincens <- 4.0
#' maxcens <- 8.4
#' evtprob <- eventProb(rateC, hrA, hrB, hrAB, mincens, maxcens)
#' avgprob <- evtprob$avgprob
#' probAB_C <- evtprob$probAB_C
#' probA_C <- evtprob$probA_C
#' dig <- 2
#' alpha <- 0.05
#' corAa  <- 1/sqrt(2)
#' corAab <- 1/sqrt(2)
#' coraab <- 1/2
#' crit13 <- crit2x2(corAa, corAab, coraab, dig, alpha)$crit13
#' n <- 4600
#' power13_13_13(n, hrA, hrB, hrAB, avgprob, probA_C, probAB_C,
#'   crit13, dig, cormat12 = matrix(c(1, sqrt(0.5), sqrt(0.5), 1), byrow = TRUE,
#'   nrow = 2), cormat23 = matrix(c(1, 0.5, 0.5, 1), byrow = TRUE, nrow = 2),
#'   cormat123 = matrix(c(1, sqrt(0.5), sqrt(0.5), sqrt(0.5), 1, 0.5,
#'   sqrt(0.5), 0.5, 1), byrow=TRUE, nrow = 3), niter = 1, abseps = 1e-03)
#'
#' # $poweroverA
#' # [1] 0.5861992
#'
#' # $powerA
#' # [1] 0.5817954
#'
#' # $powerAB
#' # [1] 0.9071236
#'
#' # $power13.13.13
#' # [1] 0.9302078

power13_13_13 <- function(n, hrA, hrB, hrAB, avgprob, probA_C, probAB_C,
                          crit13, dig,
                          cormat12 = matrix(c(1, sqrt(0.5),
                                              sqrt(0.5), 1), byrow = T, nrow = 2),
                                cormat23 = matrix(c(1, 0.5,
                                              0.5, 1), byrow = T, nrow = 2),
                                cormat123 = matrix(c(1, sqrt(0.5), sqrt(0.5),
                                              sqrt(0.5), 1, 0.5,
                                              sqrt(0.5), 0.5, 1), byrow=T, nrow = 3),
                                niter = 5, abseps = 1e-03)
{
  alpha <- 2 * pnorm(crit13)
  muoverA <- (log(hrA) + 0.5 * log(hrAB/(hrA*hrB)))* sqrt((n/4) * avgprob)
  muA <- log(hrA) * sqrt((n/8) * probA_C)
  muAB <- log(hrAB) * sqrt((n/8) * probAB_C)
  # Compute power for overall A effect
  poweroverA <- strLgrkPower(n, hrA, hrB, hrAB, avgprob, dig, alpha)$power
  # Compute power for simple A effect
  powerA <- lgrkPower(hrA, (n/2) * probA_C, alpha)$power
  # Compute power for simple AB effect
  powerAB <- lgrkPower(hrAB, (n/2) * probAB_C, alpha)$power
  # compute the power that:
  # 12. Both the overall A and simple A effects are detected.
  # 13. Both the overall A and simple AB effects are detected.
  # 23. Both the simple A and simple AB effects are detected.
  # 123. The overall A, simple A, and simple AB effects are all detected.
  # Use pmvnorm to compute the power to detect overall A and simple AB effects.
  # Do this niter times to average out the randomness in pmvnorm.
  # Previous versions of crit2x2 set a random seed here
  # to be used in conjunction with the pmvnorm call.  CRAN
  # suggested that this be omitted.
  # set.seed(rseed)
  powermat <- matrix(rep(0, 4 * niter), nrow = niter)
  for(i in 1:niter){
    powermat[i, 1] <- pmvnorm(lower=-Inf, upper=c(crit13, crit13), mean=c(muoverA, muA),
                          corr=cormat12, sigma=NULL, maxpts = 25000, abseps = abseps, releps = 0)
    powermat[i, 2] <- pmvnorm(lower=-Inf, upper=c(crit13, crit13), mean=c(muoverA, muAB),
                          corr=cormat12, sigma=NULL, maxpts = 25000, abseps = abseps, releps = 0)
    powermat[i, 3] <- pmvnorm(lower=-Inf, upper = c(crit13, crit13), mean = c(muA, muAB),
                          corr=cormat23, sigma=NULL, maxpts = 25000, abseps = abseps, releps = 0)

    powermat[i, 4] <- pmvnorm(lower=-Inf, upper=c(crit13, crit13, crit13), mean=c(muoverA, muA, muAB),
                           corr=cormat123, sigma=NULL, maxpts = 25000, abseps = abseps, releps = 0)
  }
  poweraux <- apply(powermat, 2, mean)
  powerinter12 <- poweraux[1]
  powerinter13 <- poweraux[2]
  powerinter23 <- poweraux[3]
  powerinter123 <- poweraux[4]

  power13.13.13 <- poweroverA + powerA + powerAB -
          (powerinter12 + powerinter13 + powerinter23) +
          powerinter123
  list(poweroverA = poweroverA, powerA = powerA, powerAB = powerAB,
       power13.13.13 = power13.13.13)
}



