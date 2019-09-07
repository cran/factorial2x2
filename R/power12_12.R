#' Power of the 1/2-1/2 procedure
#'
#' Computes the power of the 1/2-1/2 procedure, that is, the power to
#' detect the simple A effect or the simple AB effect.
#'
#' @param  n  total subjects with n/4 subjects in each of the C, A, B, and AB groups
#' @param  hrA  group A to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrAB  group AB to group C hazard ratio; \code{hrAB} < 1 corresponds to group AB superiority
#' @param  probA_C  event probability averaged across the A and C groups
#' @param  probAB_C  event probability averaged across the AB and C groups
#' @param  crit12  logrank statistic critical value for both the simple A and simple AB effects
#' @param  cormat  asymptotic correlation matrix for the simple A and simple AB logrank statistics
#' @param  niter  number of times we call \code{pmvnorm} to average out its randomness
#' @param  abseps  \code{abseps} setting in the \code{pmvnorm} call
#' @return \item{poweroverA }{power to detect the overall A effect}
#' @return \item{powerA }{power to detect the simple A effect}
#' @return \item{powerAB }{power to detect the simple AB effect}
#' @return \item{power12.12 }{power to detect the simple A or simple AB effects, i.e.,
#' power of the 1/2-1/2 procedure}
#' @details   For a 2-by-2 factorial design, this function computes
#' the probability that either the simple A or the simple AB logrank statistics
#' reject their null hypotheses using a \code{crit12} critical value.
#' When the two-sided familywise type I error is 0.05, we may use
#' \code{\link{crit2x2}} to compute \code{crit12} = -2.22 which corresponds
#' to a 0.0264 two-sided significance level.  This is described in
#' Leifer, Troendle, et al. (2019).
#' The \code{pmvnorm} function
#' from the \code{mvtnorm} package is used to calculate
#' the power that both (intersection) the simple A and simple B effects are detected.
#' \code{pmvnorm} uses a random seed in its algorithm.
#' To smooth out the randomness,  \code{pmvnorm} is called
#' \code{niter} times and the average value over the \code{niter} calls is taken to be that power.
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @references Slud, E.V. Analysis of factorial survival experiments. Biometrics. 1994; 50: 25-38.
#' @export power12_12
#'
#' @seealso \code{crit2x2}, \code{lgrkPower}, \code{pmvnorm}
#'
#' @examples
#' # Corresponds to scenario 5 in Table 2 from Leifer, Troendle, et al. (2019).
#' rateC <- 0.0445  # one-year C group event rate
#' hrA <- 0.80
#' hrB <- 0.80
#' hrAB <- 0.72
#' mincens <- 4.0
#' maxcens <- 8.4
#' evtprob <- eventProb(rateC, hrA, hrB, hrAB, mincens, maxcens)
#' probA_C <- evtprob$probA_C
#' probAB_C <- evtprob$probAB_C
#' corAa  <- 1/sqrt(2)
#' corAab <- 1/sqrt(2)
#' coraab <- 1/2
#' dig <- 2
#' alpha <- 0.05
#' crit12 <- crit2x2(corAa, corAab, coraab, dig, alpha)$crit12
#' n <- 4600
#' power12_12(n, hrA, hrAB, probA_C, probAB_C,
#'   crit12, cormat = matrix(c(1,0.5,0.5,1), byrow = TRUE, nrow = 2),
#'   niter = 1, abseps = 1e-03)
#'
#' # $powerA
#' # [1] 0.6203837
#'
#' # $powerAB
#' # [1] 0.9226679
#'
#' # $powerAandAB
#' # [1] 0.6018828
#'
#' # $power12.12
#' # [1] 0.9411688
#'
#'
power12_12 <- function(n, hrA, hrAB, probA_C, probAB_C,
                         crit12,
                         cormat = matrix(c(1,0.5,0.5,1), byrow =TRUE,
                         nrow =2), niter = 5, abseps = 1e-03)
{
  alphaA <- 2 * pnorm(crit12)
  alphaAB <- 2 * pnorm(crit12)
  muA <- log(hrA) * sqrt((n/8) * probA_C)
  muAB <- log(hrAB) * sqrt((n/8) * probAB_C)
  # Use pmvnorm to compute the power to detect
  # both (intersection) the simple A and simple AB effects.
  # Do this niter times to average out the randomness in pmvnorm.
  # Previous versions of crit2x2 set a random seed here
  # to be used in conjunction with the pmvnorm call.  CRAN
  # suggested that this be omitted.
  # set.seed(rseed)
  powervec <- rep(0, niter)
  for(i in 1:niter){
  powervec[i] <- pmvnorm(lower=-Inf, upper= rep(crit12, 2), mean=c(muA, muAB),
                               corr=cormat, sigma=NULL, maxpts = 25000,
                               abseps = abseps, releps = 0)
  }
  powerAandAB <- mean(powervec)
  # compute power to detect the simple A effect
  powerA <- lgrkPower(hrA, (n/2) * probA_C, alpha = alphaA)$power
  # compute power to detect the simple AB effect
  powerAB <- lgrkPower(hrAB, (n/2) * probAB_C, alpha = alphaAB)$power
  # compute power of the 1/2-1/2 procedure to detect either a
  # simple A or simple AB effect
  power12.12 <- powerA + powerAB - powerAandAB
  list(powerA = powerA, powerAB = powerAB,
       powerAandAB = powerAandAB, power12.12 = power12.12)
}


