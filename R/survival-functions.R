#' Tests for time-to-event outcomes
#'
#' @description performs a test for right-censored data. It uses the Weighted Kaplan-Meier family of statistics for testing the differences of two survival curves.
#'
#' @param time The observed time.
#' @param status The status indicator, normally 0=alive, 1=dead.
#' @param treat The treatment-group indicator, normally 0=control, 1=intervention.
#' @param tau follow-up. Default NULL denoting the last time in which both groups had patients at risk.
#' @param rho A scalar parameter that controls the type of test (see Weights).
#' @param gam A scalar parameter that controls the type of test (see Weights).
#' @param eta A scalar parameter that controls the type of test (see Weights).
#' @param var_est indicates the variance estimate to use ('Pooled' or 'Unpooled')
#'
#' @export
#'
#' @return List: standardized statistic, statistic and variance.
#' @author Marta Bofill Roig
#'

survtest <- function(time, status, treat, tau=NULL, rho=0, gam=0, eta=1,var_est='Unpooled'){
  # require(zoo) # 'rollmean' function
  # require(survival)

  db=cbind.data.frame(time=time, status=status, treat=treat)

  db1=subset(db,db$treat==1)
  db0=subset(db,db$treat==0)

  n1=dim(db1)[1]
  n0=dim(db0)[1]
  n=n0+n1

  # KAPLAN-MEIER ESTIMATORS
  ######################################

  # Estimation of the KM curve
  km1 <- survfit(Surv(time=time,event=status)~1, data=db1)
  km0 <- survfit(Surv(time=time,event=status)~1, data=db0)

  # Estimation of the pooled KM curve
  kmp <- survfit(Surv(time, status) ~ 1, data = db)

  # Estimation of the KM curve for the time-to-censoring
  censkm1 <- survfit(Surv(time=time,status==0)~1, data = db1)
  censkm0 <- survfit(Surv(time=time,status==0)~1, data = db0)

  # Kaplan-Meier function
  km1_f <- stepfun(km1$time, c(1, km1$surv))
  km0_f <- stepfun(km0$time, c(1, km0$surv))
  kmp_f <- stepfun(kmp$time, c(1, kmp$surv))

  # Censoring Kaplan-Meier function
  censkm1_f <- stepfun(censkm1$time, c(1, censkm1$surv))
  censkm0_f <- stepfun(censkm0$time, c(1, censkm0$surv))

  # TIMEPOINTS FOR INTEGRALS' COMPUTATION
  ###########################################

  # Failure and censoring times
  fail_times <- sort(db$time)

  l1=length(km1$time)
  l0=length(km0$time)

  # Define the last time
  if(is.null(tau)){
    fail_times <- sort(c(db$time))
  }else{
    fail_times <- sort(c(db$time, tau))
    fail_times <- fail_times[fail_times<=tau]
  }
  # Delete duplicates
  fail_times <- fail_times[!duplicated(fail_times)]
  # Define the midpoint between failures times
  midfail_times <- rollmean(c(0, fail_times), 2)
  # number of distinct ordered observed failures times
  l=length(fail_times)

  # WEIGHTED KAPLAN-MEIER STATISTIC
  #######################################################

  # Kaplan-Meier estimates' values at midpoints between event times (failure and censoring times)
  KM0 <- km0_f(fail_times)
  KM1 <- km1_f(fail_times)
  preKM0 <- km0_f(midfail_times)
  preKM1 <- km1_f(midfail_times)
  preKM <- kmp_f(midfail_times)
  preKM0_cens <- censkm0_f(midfail_times)
  preKM1_cens <- censkm1_f(midfail_times)

  KMpooled <- kmp_f(fail_times)
  preKMpooled <- kmp_f(midfail_times)

  # Compute the weight function
  ######################################

  # weight function proposed by Pepe-Fleming
  w <- ifelse(preKM0_cens+preKM1_cens == 0, 0, (n*preKM0_cens*preKM1_cens)/(n0*preKM0_cens + n1*preKM1_cens))
  # weight function proposed by Fleming-Harrington
  f <- ifelse(preKM0+preKM1 == 0, 0, (n*preKM0*preKM1)/(n0*preKM0 + n1*preKM1))
  # Define the weight function as a product of w and f
  weight <- w^(eta)*f^rho*(1-f)^gam

  # Weighted Kaplan-Meier Statistic
  ######################################
  WKM <- sqrt((n0*n1)/n)*sum(weight*(preKM1-preKM0)*diff(c(0, fail_times)), na.rm = TRUE)

  # Variance computation
  ######################################
  if(var_est=='Unpooled'){
    # Kaplan-Meier jumps
    # group 0
    KM0_jumps <- diff(c(1,km0_f(fail_times)))
    # group 1
    KM1_jumps <- diff(c(1,km1_f(fail_times)))

    # Calculate the integral: int_{t,tau}(weight * surv)
    # group 0
    Integral0 <- cumsum(diff(c(0, fail_times)) * weight * km0_f(fail_times))
    Int0 <- (Integral0[l] - Integral0)
    # group 1
    Integral1 <- cumsum(diff(c(0, fail_times)) * weight * km1_f(fail_times))
    Int1 <- (Integral1[l] - Integral1)

    # Calculate the estimated (unpooled) variance
    var_inside0 = ifelse(preKM0*KM0==0, 0, Int0^2*(w)^(-1)*KM0_jumps/(preKM0*KM0))
    var_inside1 = ifelse(preKM1*KM1==0, 0, Int1^2*(w)^(-1)*KM1_jumps/(preKM1*KM1))

    variance = (-n1*(sum(var_inside0,na.rm = TRUE))-
                  n0*(sum(var_inside1,na.rm = TRUE)))/n

  }else{
    KMp_jumpsaux <- c(1,kmp_f(fail_times))
    KMp_jumps <- as.numeric(l)

    for(i in 2:(l+1)){
      KMp_jumps[i-1] = KMp_jumpsaux[i]-KMp_jumpsaux[i-1]
    }

    # Calculate the integral: int_{t,tau}(weight * surv)
    Integral1 <- cumsum(diff(c(0, fail_times)) * weight * KMpooled)
    Int <- (Integral1[l] - Integral1)

    # Calculate the estimated (pooled) variance
    var_inside = ifelse(preKMpooled*KMpooled==0, 0, Int^2*(w)^(-1)*KMp_jumps/(preKMpooled*KMpooled))

    variance =  -sum(var_inside,na.rm = TRUE)
  }

  stdev <- sqrt(variance)

  # Standardized statistic
  WKMtest <- WKM/stdev

  if(is.null(tau)){
    return(list=c(Test=WKMtest,Us=WKM,sd=stdev))
  }else{
    return(list=c(Test=WKMtest,Us=WKM,sd=stdev,tau=round(tau,2)))
  }

}


##################################################################################



