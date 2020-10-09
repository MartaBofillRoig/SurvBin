#' L-Statistics for binary and time-to-event outcomes
#'
#' @description performs a L-statistic for right-censored and binary data.
#'
#' @param time The observed time.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param binary 0=nonresponse, 1=response.
#' @param treat The treatment-group indicator, normally 0=control, 1=intervention.
#' @param tau0 starting follow-up (survival outcome). Default =0.
#' @param tau follow-up. Default NULL denoting the last time in which both groups had patients at risk.
#' @param taub time-point at which the binary endpoint is evaluated. If NULL, the binary endpoint is evaluated at tau.
#' @param rho A scalar parameter that controls the type of test (see Weights).
#' @param gam A scalar parameter that controls the type of test (see Weights).
#' @param eta A scalar parameter that controls the type of test (see Weights).
#' @param wb A scalar parameter that controls the type of test (see Weights).
#' @param ws A scalar parameter that controls the type of test (see Weights).
#' @param var_est indicates the variance estimate to use ('Pooled' or 'Unpooled')
#'
#' @export
#'
#' @return List: standardized statistic, statistic and variance.
#' @author Marta Bofill Roig

lstats <- function(time, status, binary, treat, tau0=0, tau=NULL, taub=NULL, rho=0, gam=0, eta=1, wb=0.5, ws=0.5, var_est="Unpooled"){

  db=cbind.data.frame(time=time, status=status, binary=binary, treat=treat)

  db1=subset(db,db$treat==1)
  db0=subset(db,db$treat==0)

  n1=dim(db1)[1]
  n0=dim(db0)[1]
  n=n0+n1

  # KAPLAN-MEIER ESTIMATORS
  ######################################

  B <- bintest(db$binary, db$treat, var_est)
  test_b <- B[1]
  sigma_b <- B[3]

  S <- survtest(db$time, db$status, db$treat, tau, rho, gam, eta, var_est)
  test_s <- S[1]
  sigma_s <- S[3]

  sigma_sb <- survbinCov(db$time,db$status,db$binary,db$treat, tau0, tau, taub, rho, gam, eta, var_est)

  # Calculate the statistic
  u_bs = wb*test_b + ws*test_s

  # Calculate the estimated variance
  variance =  wb^2 + ws^2 + 2*wb*ws*sigma_sb/(sigma_b*sigma_s)

  stdev <- sqrt(variance)
  test_l <- u_bs/stdev

  # return(list=c(LTest=test_l,Statistic=u_bs,sd=stdev))

  output <- data.frame(Parameter=c("(Standardized) L-Test","L-Test", "Standard deviation"),
                       Value=c(test_l, u_bs, stdev))

  return(output)
}

##################################################################################


