#' L-Statistics for binary and time-to-event outcomes (bootstrap variance)
#'
#' @description performs a L-statistic for right-censored and binary data (bootstrap variance estimator).
#'
#' @param time The observed time.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param binary 0=nonresponse, 1=response.
#' @param treat The treatment-group indicator, normally 0=control, 1=intervention.
#' @param tau follow-up. Default NULL denoting the last time in which both groups had patients at risk.
#' @param rho A scalar parameter that controls the type of test (see Weights).
#' @param gam A scalar parameter that controls the type of test (see Weights).
#' @param eta A scalar parameter that controls the type of test (see Weights).
#' @param wb A scalar parameter that controls the type of test (see Weights).
#' @param ws A scalar parameter that controls the type of test (see Weights).
#' @param Boot number of bootstrap samples (default Boot=50).
#'
#' @export
#'
#' @return List: standardized statistic, statistic and variance.
#' @author Marta Bofill Roig

lstats_boots <- function(time, status, binary, treat, tau, rho=0, gam=0, eta=1, wb=0.5, ws=0.5, Boot=50){

  db=cbind.data.frame(time=time, status=status, binary=binary, treat=treat)

  db1=subset(db,db$treat==1)
  db0=subset(db,db$treat==0)

  n1=dim(db1)[1]
  n0=dim(db0)[1]
  n=n0+n1

  # KAPLAN-MEIER ESTIMATORS
  ######################################

  B <- bintest(db$binary, db$treat)
  test_b <- B[1]
  sigma_b <- B[3]

  S <- survtest(db$time, db$status, db$treat, tau, rho, gam, eta)
  test_s <- S[1]
  sigma_s <- S[3]

  # Calculate the statistic
  u_bs = wb*test_b + ws*test_s

  # bootstrap
  # B=50
  ubs_boot <- data.frame(0)
  for(i in 1:Boot){
    daux <- db[sample(1:nrow(db),nrow(db),rep=TRUE),]
    btt <- bintest(daux$binary, daux$treat)[1]
    stt <- survtest(daux$time, daux$status, daux$treat)[1]
    ubs_boot[i] = wb*btt + ws*stt
  }
  stdev <- sd(ubs_boot)
  test_l <- u_bs/stdev

  # return(list=c(LTest=test_l,Statistic=u_bs,sd=stdev))

  output <- data.frame(Parameter=c("(Standardized) L-Test","L-Test", "Standard deviation"),
                       Value=c(test_l, u_bs, stdev))

  output_bin <- data.frame(Parameter=c("Standardized L-Test","Binary Test", "Standard deviation"),
                           Value=c(B[1], B[2], B[3]))

  output_surv <- data.frame(Parameter=c("Standardized Test","Survival Test", "Standard deviation"),
                            Value=c(S[1], S[2], S[3]))

  return(list(LTest=output,Binary_Tests=output_bin,Survival_Tests=output_surv))
}

##################################################################################


