#' Power of the Proportional Allocation 2 procedure
#'
#' Computes the Proportional Allocation 2 procedure's power to
#' detect the overall A effect or the simple AB effect, respectively.
#'
#' @param  n  total subjects with n/4 subjects in each of the C, A, B, and AB groups
#' @param  hrA  group A to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrB  group B to group C hazard ratio; \code{hrA} < 1 corresponds to group A superiority
#' @param  hrAB  group AB to group C hazard ratio; \code{hrAB} < 1 corresponds to group AB superiority
#' @param  avgprob  event probability averaged across the C, A, B, and AB groups
#' @param  critPA2A  rejection critical value for the overall A stratified logrank statistic
#' @param  critPA2ab  rejection critical value for the simple AB ordinary logrank statistic
#' @param  dig number of decimal places to which we \code{\link{roundDown}} the critical value
#' for the overall A test as calculated in \code{\link{powerPA2}}
#' by \code{\link{strLgrkPower}}
#' @param  probAB_C  event probability averaged across the AB and C groups
#' @return \item{powerPA2overallA }{power to detect the overall A effect}
#' @return \item{powerPA2simpleAB }{power to detect the simple AB effect}
#'
#' @details  The Proportional Allocation 2 procedure uses a two-sided
#' 2/3 * alpha significance level to test the overall A effect and the
#' remaining Dunnett-corrected type 1 error to thest the simple AB effect.
#' When  the familywise error is alpha = 0.05, this corresponds to a
#' critical value \code{critPA2A} = -2.13.
#' Then \code{\link{crit2x2}} is used to compute a critical value
#' \code{critPA2ab} = -2.24 to test the simple AB effect.  This corresponds to
#' a two-sided 0.0251 significance level.  This controls the
#' asymptotic familywise type I error for the two hypothesis tests at the
#' two-sided 0.05 level.  This is because of the \code{1/sqrt(2)} asymptotic
#' correlation between the logrank test statistics for the overall A
#' and simple AB effects (Slud, 1994).  The overall A effect's significance
#' level 2/3 * 0.05 is prespecified and the simple AB effect's significance
#' level 0.0251 is computed using \code{crit2x2}.
#' @author Eric Leifer, James Troendle
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2020). Submitted.
#' @references Lin, D-Y., Gong, J., Gallo, P., et al. Simultaneous inference on treatment effects
#' in survival studies with factorial designs. Biometrics. 2016; 72: 1078-1085.
#' @references Slud, E.V. Analysis of factorial survival experiments. Biometrics. 1994; 50: 25-38.
#' @export powerPA2
#' @seealso  \code{crit2x2}, \code{eventProb}, \code{lgrkPower}, \code{strLgrkPower}
#' @examples
#'  # Corresponds to scenario 4 in Table 2 from Leifer, Troendle, et al. (2020).
#'  rateC <- 0.0445  # one-year C group event rate
#'  hrA <- 0.80
#'  hrB <- 0.80
#'  hrAB <- 0.72
#'  mincens <- 4.0
#'  maxcens <- 8.4
#'  evtprob <- eventProb(rateC, hrA, hrB, hrAB, mincens, maxcens)
#'  avgprob <- evtprob$avgprob
#'  probAB_C <- evtprob$probAB_C
#'  dig <- 2
#'  alpha <- 0.05
#'  corAa  <- 1/sqrt(2)
#'  corAab <- 1/sqrt(2)
#'  coraab <- 1/2
#'  critvals <- crit2x2(corAa, corAab, coraab, dig, alpha)
#'  critPA2A <- critvals$critPA2A
#'  critPA2ab <- critvals$critPA2ab
#'  n <- 4600
#'  powerPA2(n, hrA, hrB, hrAB, avgprob, probAB_C,
#'             critPA2A, critPA2ab, dig)
#' # $powerPA2overallA
#' # [1] 0.6582819
#'
#' # $powerPA2simpleAB
#' # [1] 0.9197286


powerPA2 <- function(n, hrA, hrB, hrAB, avgprob, probAB_C,
                critPA2A, critPA2ab, dig)
{
  alphaoverA <- 2 * pnorm(critPA2A)
  alphaAB <- 2 * pnorm(critPA2ab)
  poweroverA <- strLgrkPower(n, hrA, hrB, hrAB,
                             avgprob, dig, alpha = alphaoverA)$power
  powerAB <- lgrkPower(hrAB, (n/2) * probAB_C, alpha = alphaAB)$power
  list(powerPA2overallA = poweroverA, powerPA2simpleAB = powerAB)
}

