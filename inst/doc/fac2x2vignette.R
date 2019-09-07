## ----setup, include = FALSE----------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

## ---- eval = F, echo = F, results = "hide"-------------------------------
#  # When I use knitr to run this whole file, I need to
#  # "require" all the packages that are needed to run the code
#  # since knitr uses a new environment.  This is in contrast to
#  # running the code chunks separately which use all of the packages
#  # which are listed in the factorial2x2 package's DESCRIPTION file
#  # on the DEPENDS line
#  
#  # require(factorial2x2)
#  # require(survival)
#  # require(mvtnorm)

## ---- eval =   F---------------------------------------------------------
#  time <- simdata[, 'time']    # follow-up time
#  event <- simdata[, 'event']  # event indicator
#  indA <- simdata[, 'indA']    # treatment A indicator
#  indB <- simdata[, 'indB']    # treatment B indicator
#  fac2x2analyze(time, event, indA, indB, covmat, alpha = 0.05, dig = 3, niter = 5)
#  # simdata[, 6:10] corresponds to the baseline covariates which include
#  # a history of cardiovascular disease (yes/no) and four indicator
#  # variables which correspond to which of 5 clinical centers enrolled each of the participants
#  $hrA
#  [1] 0.8895135          # overall A effect hazard ratio (HR)
#  
#  $ciA
#  [1] 0.786823 1.005607  # 95% CI for overall A effect HR
#  
#  $pvalA
#  [1] 0.06139083         # p-value for overall A effect HR
#  
#  $hra
#  [1] 0.8096082          # simple A effect HR
#  
#  $cia
#  [1] 0.6832791 0.9592939 # 95% CI for simple A effect HR
#  
#  $pvala
#  [1] 0.01468184         # p-value for simple A effect HR
#  
#  $hrab
#  [1] 0.7583061          # simple AB effect HR
#  
#  $ciab
#  [1] 0.6389355 0.8999785 # 95% CI fo simple A effect HR
#  
#  $pvalab
#  [1] 0.001545967       # p-value for simple AB effect HR
#  
#  $crit23A
#  [1] -2.129            # critical value for the overall A effect for 2/3-1/3 procedire
#  
#  $sig23A
#  [1] 0.03333333        # significance criterion for overall A effect for 2/3-1/3 procedure
#  
#  $crit23ab
#  [1] -2.233            # critical value for the simple AB effect for 2/3-1/3 procedure
#  
#  $sig23ab
#  [1] 0.0256049         # signficance criterion of simple AB effect for 2/3-1/3 procedure
#  
#  $result23
#  [1] "accept overall A" "reject simple AB"   # hypothesis testing results for 2/3-1/3 proceudre
#  
#  $crit13
#  [1] -2.31            # critical value for all effects for the 1/3-1/3-1/3 procedure
#  
#  $sig13
#  [1] 0.02091404        # significance criterion for all effects for 1/3-1/3-1/3 procedure
#  
#  $result13
#  [1] "accept overall A" "reject simple A"  "reject simple AB"  # hypothesis testing results
#  
#  $crit12
#  [1] -2.217            # critical value for all effects for the 1/2-1/2 procedure
#  
#  $sig12
#  [1] 0.02665078       # significance criterion all effects for 1/2-1/2 procedure
#  
#  $result12
#  [1] "reject simple A"  "reject simple AB"   # hypothesis testing results

## ---- eval =   F---------------------------------------------------------
#  # read the COMBINE data into an R data frame
#  Combine <- read.table("c:\\combine_data.txt", header = T, nrows = 1226, na.strings ="",
#                        stringsAsFactors= T)
#  dim(Combine)
#  [1] 1226    9
#  
#  dimnames(Combine)[[2]]
#  [1] "ID"         "AGE"        "GENDER"     "T0_PDA"     "NALTREXONE"
#  [6] "THERAPY"    "site"       "relapse"    "futime"
#  
#  # create the baseline covariate variables
#  T0_PDA <- Combine[,"T0_PDA"]            # baseline percentage of days abstinent
#  site_1 <- Combine[,"site"] == "site_1"  # research site indicator variables
#  site_2 <- Combine[,"site"] == "site_2"
#  site_3 <- Combine[,"site"] == "site_3"
#  site_4 <- Combine[,"site"] == "site_4"
#  site_5 <- Combine[,"site"] == "site_5"
#  site_6 <- Combine[,"site"] == "site_6"
#  site_7 <- Combine[,"site"] == "site_7"
#  site_8 <- Combine[,"site"] == "site_8"
#  site_9 <- Combine[,"site"] == "site_9"
#  site_10 <- Combine[,"site"] == "site_10"
#  
#  # combine the covariates into a single covariate matrix
#  CombineCovMat <- cbind(T0_PDA, site_1, site_2, site_3, site_4, site_5, site_6,
#                           site_7, site_8, site_9, site_10)
#  
#  # define the other required variables
#  relapse <- Combine[,"relapse"]         # heavy drinking relapse indicator
#  futime <- Combine[,"futime"]           # time to first heavy drinking day or censoring
#  NALTREXONE <- Combine[,"NALTREXONE"]   # received naltrexone indicator
#  THERAPY <- Combine[,"THERAPY"]         # received cognitive behavioral intervention (CBI) indicator
#  
#  # reproduce the COMBINE analysis using fac2x2analyze
#  fac2x2analyze(futime, relapse, NALTREXONE, THERAPY, CombineCovMat, alpha = 0.025, dig = 4)
#  
#  $loghrA
#  [1] -0.0847782              # log hazard rato estimate for the overall effect of naltrexone
#  
#  $seA
#  [1] 0.06854294              # std error of the log HR estimate for the overall effect of naltrexone
#  
#  $ZstatA
#  [1] -1.236863               # Z-statistic for the overall effect of naltrexone
#  
#  
#  $pvalA
#  [1] 0.2161381               # p-value for the overall effect of naltrexone
#  
#  $hrA
#  [1] 0.918716                # hazard ratio estimate for the overall effect of naltrexone
#  
#  $ciA
#  [1] 0.8032234 1.0508149     # corresponding 95% confidence interval
#  
#  $loghra
#  [1] -0.2517618              # log hazard rato estimate for the simple effect of naltrexone
#  
#  $sea
#  [1] 0.09786137              # std error of the log HR estimate for the simple effect of naltrexone
#  
#  $Zstata
#  [1] -2.572637               # Z-statistic for the simple effect of naltrexone
#  
#  $pvala
#  [1] 0.0100927               # p-value for the simple effect of naltrexone
#  
#  $hra
#  [1] 0.7774299               # hazard ratio estimate for the simple effect of naltrexone
#  
#  $cia
#  [1] 0.6417413 0.9418083     # corresponding 95% confidence interval
#  
#  $loghrab
#  [1] -0.09132675             # log hazard ratio estimate for the simple effect of naltrexone and CBI
#  
#  $seab
#  [1] 0.09553005              # std error of the log HR estimate for the simple effect of naltrexone
#                              # and CBI
#  
#  $Zstatab
#  [1] -0.9560003              # Z-statistic for the simple effect of naltrexone and CBI
#  
#  $pvalab
#  [1] 0.3390721               # p-value for the simple effect of naltrexone and CBI
#  
#  $hrab
#  [1] 0.9127194               # hazard ratio estimate for the simple effect of naltrexone and CBI
#  
#  $ciab
#  [1] 0.7568686 1.1006624     # corresponding 95% confidence interval
#  
#  $crit13
#  [1] -2.5811                 # critical value for the three tests in Table 4;
#                              # slightly larger in absolute terms than the
#                              # critical value -2.573 reported on p.1083 of Lin, Gong, et al.
#  $sig13
#  [1] 0.009848605
#  
#  $result13
#  [1] "accept overall A" "accept simple A"  "accept simple AB"
#  
#  $crit12
#  [1] -2.2154
#  
#  $sig12
#  [1] 0.02673262
#  
#  $result12
#  [1] "reject simple A"  "accept simple AB"
#  
#  $corAa
#  [1] 0.6727651
#  
#  $corAab
#  [1] 0.7078624
#  
#  $coraab
#  [1] 0.4691156

## ---- eval =   F---------------------------------------------------------
#    n <- 4600          # total sample size
#    rateC <- 0.0445    # one year event rate in the control group
#    hrA <- 0.80        # simple A effect hazard ratio
#    hrB <- 0.80        # simple B effect hazard ratio
#    hrAB <- 0.72       # simple AB effect hazard ratio
#    mincens <- 4.0     # minimum censoring time in years
#    maxcens <- 8.4     # maximum censoring time in years
#  
#    fac2x2design(n, rateC, hrA, hrB, hrAB, mincens, maxcens, dig = 2, alpha = 0.05)
#    $powerA
#    [1] 0.7182932      # power to detect the overall A effect at the two-sided 0.05 level
#  
#    $power23.13
#    [1] 0.9290271      # power to detect the overall A or simple AB effects using the
#                       # 2/3-1/3 procedure
#  
#    $power13.13.13
#    [1] 0.9302084      # power to detect the overall A, simple A, or simple AB effects using
#                       # the 1/3-1/3-1/3 procedure
#  
#    $power12.12
#    [1] 0.9411688      # power to detect the simple A or simple AB effects using the
#                       # 1/2-1/2 procedure
#  
#    $events            # expected number of events
#    [1] 954.8738
#  
#    $evtprob          # event probabilities for the C, A, B, and AB groups, respectively
#    probC     probA     probB    probAB
#    0.2446365 0.2012540 0.2012540 0.1831806

