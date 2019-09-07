#' Unstratified (ordinary) logrank power
#'
#' Computes the power for the unstratified (ordinary) logrank statistic
#' for two group comparison.
#'
#' @param hr  hazard ratio
#' @param nevent  expected number of events
#' @param alpha  two-sided significance level
#' @param rprob  randomization probability
#'
#' @details Uses the formula at the bottom of p.317 from Schoenfeld (Biometrika, 1981)
#' where the beta should be 1 - beta.
#' @return  \item{power }{logrank power}
#' @references Schoenfeld, D. The asymptotic properties of nonparametric tests for comparing
#'   survival distributions. Biometrika. 1981; 68: 316-319.
#' @author Eric Leifer, James Troendle
#' @export lgrkPower
#' @examples
#'  hr <- 0.5
#'  nevent <- 98
#'  lgrkPower(hr, nevent, alpha = 0.05,  rprob = 0.5)
#'  # $power
#'  # [1] 0.9293463
#'
lgrkPower <- function(hr, nevent, alpha = 0.05, rprob = 0.5){

  # September 15, 2011
  # This function gives the power for two-armed logrank test.
  # hr = hazard ratio
  # nevent = number of events
  # alpha = two-sided significance level
  # rprob = treatment arm randomization probability
  # with a log hazard ratio of LogHazRat, a two-sided alpha,
  # This formula is taken from the bottom of
  # p.317 of Schoenfeld (Biometrika, 1981) WHERE
  # SCHOENFELD'S BETA SHOULD BE 1 - BETA.

  lhr <- log(hr)
  critval <- qnorm(1 - alpha/2)
  power <- pnorm( sqrt(nevent * rprob * (1 - rprob)) * abs(lhr) - critval)
  list(power = power)
}

