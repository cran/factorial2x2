#' Power of the 2/3-1/3 procedure
#'
#' Computes the power of the 2/3-1/3 procedure, that is, the power to
#' detect the overall A effect or the simple AB effect.
#'
#' @param  n  total subjects with n/4 subjects in each of the C, A, B, and AB groups
#' @param  hrA  group A to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrB  group B to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrAB  group AB to group C hazard ratio; \code{hrAB} < 1 corresponds to group AB superiority
#' @param  avgprob  event probability averaged across the C, A, B, and AB groups
#' @param  crit23A  rejection critical value for the overall A stratified logrank statistic
#' @param  crit23ab  rejection critical value for the simple AB ordinary logrank statistic
#' @param  dig number of decimal places to which we \code{\link{roundDown}} the critical value
#' for the overall A test as calculated in \code{\link{power23_13}}
#' by \code{\link{strLgrkPower}}
#' @param  probAB_C  event probability averaged across the AB and C groups
#' @param  cormat  asymptotic correlation matrix for the overall A and simple AB logrank statistics
#' @param  niter  number of times we call \code{pmvnorm} to average out its randomness
#' @param  abseps  \code{abseps} setting in the \code{pmvnorm} call
#' @return \item{poweroverA }{power to detect the overall A effect}
#' @return \item{powerAB }{power to detect the simple AB effect}
#' @return \item{poweroverAandAB }{power to detect the overall A and simple AB effects}
#' @return \item{power23.13 }{power to detect the overall A or simple AB effects, i.e., power of the 2/3-1/3 procedure}
#'
#' @details  The 2/3-1/3 procedure uses a  two-sided
#' 2/3 * alpha = 0.033 significance level to test the overall A effect.
#' When  the familywise error is alpha = 0.05, this corresponds to a
#' critical value \code{crit23A} = -2.13.
#' Then \code{\link{crit2x2}} is used to compute a critical value
#' \code{crit23ab} = -2.24 to test the simple AB effect.  This corresponds to
#' a two-sided 0.0251 significance level.  This controls the
#' asymptotic familywise type I error for the two hypothesis tests at the
#' two-sided 0.05 level.  This is because of the \code{1/sqrt(2)} asymptotic
#' correlation between the logrank test statistics for the overall A
#' and simple AB effects (Slud, 1994).  The overall A effect's significance
#' level 2/3 * 0.05 is prespecified and the simple AB effect's significance
#' level 0.0251 is computed using \code{crit2x2}.
#' The \code{pmvnorm} function
#' from the \code{mvtnorm} package is used to calculate
#' the power that both (intersection) the overall A and simple AB effects are detected.
#' Since this involves bivariate normal integration over an unbounded region in R^2, \code{pmvnorm}
#' uses a random seed for this computation.  To smooth out the
#' randomness, \code{pmvnorm} is called \code{niter} times and
#' the average value over the \code{niter} calls is taken to be that power.
#' @author Eric Leifer, James Troendle
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @references Slud, E.V. Analysis of factorial survival experiments. Biometrics. 1994; 50: 25-38.
#' @export power23_13
#' @seealso  \code{crit2x2}, \code{eventProb}, \code{lgrkPower}, \code{strLgrkPower}, \code{pmvnorm}
#' @examples
#'  # Corresponds to scenario 5 in Table 2 from Leifer, Troendle, et al. (2019).
#'  rateC <- 0.0445  # one-year C group event rate
#'  hrA <- 0.80
#'  hrB <- 0.80
#'  hrAB <- 0.72
#'  mincens <- 4.0
#'  maxcens <- 8.4
#'  eventvec <- eventProb(rateC, hrA, hrB, hrAB, mincens, maxcens)
#'  avgprob <- eventvec$avgprob
#'  probAB_C <- 0.5 * (eventvec$probAB + eventvec$probC)
#'  dig <- 2
#'  alpha <- 0.05
#'  corAa  <- 1/sqrt(2)
#'  corAab <- 1/sqrt(2)
#'  coraab <- 1/2
#'  critvals <- crit2x2(corAa, corAab, coraab, dig, alpha)
#'  crit23A <- critvals$crit23A
#'  crit23ab <- critvals$crit23ab
#'  n <- 4600
#'  power23_13(n, hrA, hrB, hrAB, avgprob, probAB_C,
#'             crit23A, crit23ab, dig, cormat =
#'             matrix(c(1, sqrt(0.5), sqrt(0.5), 1), byrow = TRUE,
#'             nrow = 2), niter = 1, abseps = 1e-03)
#' # $poweroverA
#' # [1] 0.6582819
#'
#' # $powerAB
#' # [1] 0.9197286
#'
#' # $poweroverAandAB
#' # [1] 0.6490042
#'
#' # $power23.13
#' # [1] 0.9290062

power23_13 <- function(n, hrA, hrB, hrAB, avgprob, probAB_C,
                crit23A, crit23ab, dig, cormat =
                matrix(c(1, sqrt(0.5), sqrt(0.5), 1), byrow = TRUE,
                nrow = 2), niter = 5, abseps = 1e-03)
{
  alphaoverA <- 2 * pnorm(crit23A)
  alphaAB <- 2 * pnorm(crit23ab)
  # muoverA is the mean of the overall A logrank statistic
  muoverA <- (log(hrA) + 0.5 * log(hrAB/(hrA*hrB)))* sqrt((n/4) * avgprob)
  # muAB is the mean of the simple AB logrank statistic
  muAB <- log(hrAB) * sqrt((n/8) * probAB_C)
  poweroverA <- strLgrkPower(n, hrA, hrB, hrAB,
                             avgprob, dig, alpha = alphaoverA)$power
  powerAB <- lgrkPower(hrAB, (n/2) * probAB_C, alpha = alphaAB)$power
  # compute the power that both (intersection) the overall A and simple AB effects are detected.
  # Use pmvnorm to compute the power to detect overall A and simple AB effects.
  # Do this niter times to average out the randomness in pmvnorm.
  # Previous versions of crit2x2 set a random seed here
  # to be used in conjunction with the pmvnorm call.  CRAN
  # suggested that this be omitted.
  # et.seed(rseed)
  powervec <- rep(0, niter)
  for(i in 1:niter){
    powervec[i] <- pmvnorm(lower=-Inf, upper=c(crit23A, crit23ab), mean=c(muoverA, muAB),
                           corr=cormat, sigma=NULL, maxpts = 25000, abseps = abseps,
                           releps = 0)
  }
  poweroverAandAB <- mean(powervec)
  power23.13 <-  poweroverA + powerAB - poweroverAandAB
  # power2 is ANOTHER, easier way to compute power union
  #power2 <- 1 - pmvnorm(lower= c(crit23A, crit23ab), upper=Inf, mean=c(muoverA, muAB),
  #                      corr=cormat, sigma=NULL, maxpts = 25000, abseps = 0.001,
  #                      releps = 0)

  list(poweroverA = poweroverA, powerAB = powerAB, poweroverAandAB =
         poweroverAandAB, power23.13 = power23.13)
}

