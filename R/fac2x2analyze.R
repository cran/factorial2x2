#' Significance testing for the 2/3-1/3, 1/3-1/3-1/3, 1/2-1/2 procedures
#'
#' Performs significance testing for the 2/3-1/3, 1/3-1/3-1/3, 1/2-1/2 procedures.
#' Also reports the hazard ratios, 95\% confidence intervals, p-values,
#' nominal significance levels, and correlations for the overall A,
#' simple A, and simple AB test statistics.
#'
#' @param time follow-up times
#' @param event event indicators (0/1)
#' @param indA treatment A indicators (0/1)
#' @param indB treatment B indicators (0/1)
#' @param covmat covariate matrix, must be non-NULL.  Factor variables MUST use 0/1 dummy variables
#' @param alpha two-sided familywise significance level
#' @param dig number of decimal places to which we \code{\link{roundDown}} the critical value
#' @param niter number of interations passed to \code{crit2x2} function call
#'
#' @return \item{loghrA }{overall A log hazard ratio}
#' @return \item{seA }{standard error of the overall A log hazard ratio}
#' @return \item{ZstatA }{Z-statistic for the overall A log hazard ratio}
#' @return \item{pvalA }{two-sided p-value for the overall hazard ratio}
#' @return \item{hrA }{overall A hazard ratio}
#' @return \item{ciA }{95\% confidence interval for the overall A hazard ratio}
#' @return \item{loghra }{simple A log hazard ratio}
#' @return \item{sea }{standard error of the simple A log hazard ratio}
#' @return \item{Zstata }{Z-statistic for the simple A log hazard ratio}
#' @return \item{pvala }{two-sided p-value for the simple A hazard ratio}
#' @return \item{hra }{simple A hazard ratio}
#' @return \item{cia }{95\% confidence interval for the simple A hazard ratio}
#' @return \item{loghrab }{simple AB log hazard ratio}
#' @return \item{seab}{standard error of the simple AB log hazard ratio}
#' @return \item{Zstatab }{Z-statistic for the simple AB log hazard ratio}
#' @return \item{pvalab }{two-sided p-value for the simple AB hazard ratio}
#' @return \item{hrab }{simple AB hazard ratio}
#' @return \item{ciab}{95\% confidence interval for the simple AB hazard ratio}
#' @return \item{crit23A }{2/3-1/3 procedure's critical value for the overall A statistic}
#' @return \item{sig23A }{2/3-1/3 procedure's p-value rejection criterion for the overall A null hypothesis}
#' @return \item{crit23ab}{2/3-1/3 procedure's critical value for the simple AB statistic}
#' @return \item{sig23ab }{2/3-1/3 procedure's p-value rejection criterion for the simple A null hypothesis}
#' @return \item{result23 }{2/3-1/3 procedure's accept/reject decisions for the overall A and simple A null hypotheses results}
#' @return \item{crit13 }{1/3-1/3-1/3 procedure's critical value for all three test statistics}
#' @return \item{sig13 }{1/3-1/3-1/3 procedure's p-value rejection criterion for the overall A, simple A, and simple AB null hypotheses}
#' @return \item{result13 }{1/3-1/3-1/3 procedure's accept/reject decisions for the overall A, simple A, and simple AB null hypotheses results}
#' @return \item{crit12 }{1/2-1/2 procedure's critical value for the simple A and AB statistics}
#' @return \item{sig12 }{1/2-1/2 procedure's p-value rejection criterion for the simple A and simple AB null hypotheses}
#' @return \item{result12 }{1/2-1/2 procedure's accept/reject decisions for the simple A and simple AB null hypotheses results}
#' @return \item{corAa }{correlation between the overall A and simple A Wald statistics}
#' @return \item{corAab }{correlation between the overall A and simple AB Wald statistics}
#' @return \item{coraab }{correlation between the simple A and simple AB Wald statistics}
#' @author Eric Leifer, James Troendle
#' @references Leifer, E.S., Troendle, J.F., Kolecki, A., Follmann, D.
#' Joint testing of overall and simple effect for the two-by-two factorial design. (2019). Submitted.
#' @export fac2x2analyze
#' @examples
#'  # First load the simulated data variables. The "simdat" file is
#'  # a 100-by-9 matrix which is loaded with the factorial2x2 package.
#'  time <- simdat[, "time"]
#'  event <- simdat[, "event"]
#'  indA <- simdat[, "indA"]
#'  indB <- simdat[, "indB"]
#'  covmat <- simdat[, 6:10]
#'  fac2x2analyze(time, event, indA, indB, covmat, alpha = 0.05, niter = 5)
#' #  $loghrA
#' # [1] 0.05613844
#'
#' # $seA
#' # [1] 0.4531521
#'
#' # $ZstatA
#' # [1] 0.1238843
#'
#' # $pvalA
#' # [1] 0.9014069
#'
#' # $hrA
#' # [1] 1.057744
#'
#' # $ciA
#' # [1] 0.4351608 2.5710556
#'
#' # $loghra
#' # [1] 0.1987329
#'
#' # $sea
#' # [1] 0.6805458
#'
#' # $Zstata
#' # [1] 0.2920198
#'
#' # $pvala
#' # [1] 0.7702714
#'
#' # $hra
#' # [1] 1.219856
#'
#' # $cia
#' # [1] 0.3213781 4.6302116
#'
#' # $loghrab
#' # [1] 0.2864932
#'
#' # $seab
#' # [1] 0.6762458
#'
#' # $Zstatab
#' # [1] 0.4236525
#'
#' # $pvalab
#' # [1] 0.6718193
#'
#' # $hrab
#' # [1] 1.331749
#'
#' # $ciab
#' # [1] 0.3538265 5.0125010
#'
#' # $crit23A
#' # [1] -2.129
#'
#' # $sig23A
#' # [1] 0.03325426
#'
#' # $crit23ab
#' # [1] -2.299
#'
#' # $sig23ab
#' # [1] 0.02150494
#'
#' # $result23
#' # [1] "accept overall A" "accept simple AB"
#'
#' # $crit13
#' # [1] -2.338
#'
#' # $sig13
#' # [1] 0.01938725
#'
#' # $result13
#' # [1] "accept overall A" "accept simple A"  "accept simple AB"
#'
#' # $crit12
#' # [1] -2.216
#'
#' # $sig12
#' # [1] 0.0266915
#'
#' # $result12
#' # [1] "accept simple A"  "accept simple AB"
#'
#' # $corAa
#' # [1] 0.6123399
#'
#' # $corAab
#' # [1] 0.5675396
#'
#' # $coraab
#' # [1] 0.4642737


fac2x2analyze <- function(time, event, indA, indB, covmat, alpha, dig = 2, niter = 5){
  # Performs the Wald significance tests for the 2/3-1/3,
  # 1/3-1/3-1/3, and 1/2-1/2 procedures.  NEED to consider the case
  # when covmat = NULL.  It calls crit2x2 which calls
  # cor2x2.
  # time =  follow-up time
  # event = event indicator: 0=censoring, 1=event
  # indA = treatment A indicator (0 = no, 1 = yes)
  # indB = treatment B indicator (0 = no, 1 = yes)
  # covmat = covariate matrix.  NOTE!! Factor variables have to
  #		use 0/1 dummy variables
  # niter = number of iterations in the crit2x2 call
  # It computes the overall A, simple A, and simple AB p-values for
  # the corresponding Cox model log hazard ratio parameter estimates
  # (where the overall A Cox model is stratified
  # on B) to determine the statistical significance of each effect for
  # each of the three procedures.

  coraux <- cor2x2(time, event, indA, indB, covmat)
  corAa <- coraux$corAa
  corAab <- coraux$corAab
  coraab <- coraux$coraab
  loghrA <- coraux$loghrA
  seA <- coraux$seA
  hrA <- coraux$hrA
  ciA <- coraux$ciA
  pvalA <- coraux$pvalA
  loghra <- coraux$loghra
  sea <- coraux$sea
  hra <- coraux$hra
  cia <- coraux$cia
  pvala <- coraux$pvala
  loghrab <- coraux$loghrab
  seab <- coraux$seab
  hrab <- coraux$hrab
  ciab <- coraux$ciab
  pvalab <- coraux$pvalab

  critaux <- crit2x2(corAa, corAab, coraab, dig = dig,
                     alpha = alpha, niter = niter)
  crit23A <- critaux$crit23A
  sig23A <- critaux$sig23A
  crit23ab <- critaux$crit23ab
  sig23ab <- critaux$sig23ab
  crit13 <- critaux$crit13
  sig13 <- critaux$sig13
  crit12 <- critaux$crit12
  sig12 <- critaux$sig12

  # Make accept/reject statements.
  accA <- "accept overall A"
  rejA <- "reject overall A"
  acca <- "accept simple A"
  reja <- "reject simple A"
  accab <- "accept simple AB"
  rejab <- "reject simple AB"

  # Statement for 2/3-1/3 procedure
  if (pvalA <= sig23A) {
    res23A <- rejA
  } else {
    res23A <- accA
  }
  if (pvalab <= sig23ab) {
    res23ab <- rejab
  } else {
    res23ab <- accab
  }
  result23 <- c(res23A, res23ab)

  # Statement for 1/3-1/3-1/3 procedure
  if (pvalA <= sig13) {
    res13A <- rejA
  } else {
    res13A <- accA
  }
  if (pvala <= sig13) {
    res13a <- reja
  } else {
    res13a <- acca
  }
  if (pvalab <= sig13) {
    res13ab <- rejab
  } else {
    res13ab <- accab
  }
  result13 <- c(res13A, res13a, res13ab)

  # Statement for 1/2-1/2 procedure
  if (pvala <= sig12) {
    res12a <- reja
  } else {
    res12a <- acca
  }
  if (pvalab <= sig12) {
    res12ab <- rejab
  } else {
    res12ab <- accab
  }
  result12 <- c(res12a, res12ab)

  list(loghrA = loghrA, seA = seA, ZstatA = loghrA/seA, pvalA = pvalA, hrA = hrA, ciA = ciA,
       loghra = loghra, sea = sea, Zstata = loghra/sea, pvala = pvala, hra = hra, cia = cia,
       loghrab = loghrab, seab = seab, Zstatab = loghrab/seab, pvalab = pvalab, hrab = hrab, ciab = ciab,
       crit23A = crit23A, sig23A = sig23A, crit23ab = crit23ab, sig23ab = sig23ab, result23 = result23,
       crit13 = crit13, sig13 = sig13, result13 = result13,
       crit12 = crit12, sig12 = sig12, result12 = result12,
       corAa = corAa, corAab = corAab, coraab = coraab
       )
}

