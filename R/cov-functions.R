#' Covariance between binary and time-to-event statistics
#'
#' @description computes the covariance between the binary and time-to-event statistics.
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
#'
#' @export
#'
#' @return Covariance.
#' @author Marta Bofill Roig

survbinCov <- function(time, status, binary, treat, tau0=0, tau=NULL, taub=NULL, rho=0, gam=0, eta=1){
  # require(zoo) # 'rollmean' function
  # require(survival)
  # require(muhaz)

  # x1=data$time; x2=data$status; x3=data$binary; x4=data$treat;
  # tau=3; taub=3
  # rho=0; gam=0; eta=1
  db=cbind.data.frame(time=time, status=status, binary=binary, treat=treat)

  db1=subset(db,db$treat==1)
  db0=subset(db,db$treat==0)

  # Subgroup with binary event =1 (responders)
  dbX=subset(db,db$binary==1)

  db1X=subset(db,db$treat==1 & db$binary==1)
  db0X=subset(db,db$treat==0 & db$binary==1)

  n1=dim(db1)[1]
  n0=dim(db0)[1]
  n=n0+n1

  n1x=dim(db1X)[1]
  n0x=dim(db0X)[1]
  nx=n0x+n1x

  # ESTIMATED PROBABILITY BINARY
  ######################################

  phat_group0 = sum(db0$binary)/n0
  phat_group1 = sum(db1$binary)/n1

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

  # Kaplan-Meier function for responders
  km1x <- survfit(Surv(time=time,event=status)~1, data=db1X)
  km0x <- survfit(Surv(time=time,event=status)~1, data=db0X)

  # Estimation of the pooled KM curve responders
  kmpx <- survfit(Surv(time, status)~1, data=dbX)

  # TIME-TO-EVENT FUNCTIONS
  ######################################

  # Kaplan-Meier functions
  km1_f <- stepfun(km1$time, c(1, km1$surv))
  km0_f <- stepfun(km0$time, c(1, km0$surv))
  kmp_f <- stepfun(kmp$time, c(1, kmp$surv))

  km1x_f <- stepfun(km1x$time, c(1, km1x$surv))
  km0x_f <- stepfun(km0x$time, c(1, km0x$surv))
  kmpx_f <- stepfun(kmpx$time, c(1, kmpx$surv))

  # Censoring Kaplan-Meier function
  censkm1_f <- stepfun(censkm1$time, c(1, censkm1$surv))
  censkm0_f <- stepfun(censkm0$time, c(1, censkm0$surv))

  # Hazard function
  # version with kernels
  fit <- muhaz(db1$time,db1$status)
  hazard1_f <- approxfun(fit$est.grid, fit$haz.est)
  fit <- muhaz(db0$time,db0$status)
  hazard0_f <- approxfun(fit$est.grid, fit$haz.est)

  # RESPONDERS
  # Hazard function for responders
  # version with kernels
  # we define those individuals that failed at time t, without having X as censored

  db1_aux=db1
  db1_aux[which(db1$status==1 & db1$binary==0),]$status=0
  db0_aux=db0
  db0_aux[which(db0$status==1 & db0$binary==0),]$status=0

  fit <- muhaz(db1_aux$time,db1_aux$status)
  hazard1X_f <- approxfun(fit$est.grid, fit$haz.est)
  fit <- muhaz(db0_aux$time,db0_aux$status)
  hazard0X_f <- approxfun(fit$est.grid, fit$haz.est)

  # TIMEPOINTS FOR INTEGRALS' COMPUTATION
  ###########################################

  # Failure and censoring times
  fail_times <- sort(db$time)

  l1=length(km1$time)
  l0=length(km0$time)

  # Define the last time
  if(is.null(tau)){
    fail_times <- sort(c(db$time))
    tau = max(fail_times)
  }else{
    fail_times <- sort(c(db$time, tau))
    fail_times <- fail_times[fail_times<=tau]
  }

  # Define times pre-taub and post-taub
  if(is.null(taub)){
    taub <- tau
    fail_times_pretaub <- fail_times
    fail_times_postaub <- 0
  }else{
    fail_times_pretaub <- fail_times[fail_times<=taub]
    fail_times_postaub <- fail_times[fail_times>taub]
  }

  # Delete duplicates
  fail_times <- fail_times[!duplicated(fail_times)]
  fail_times_pretaub <- fail_times_pretaub[!duplicated(fail_times_pretaub)]
  fail_times_postaub <- fail_times_postaub[!duplicated(fail_times_postaub)]

  # Define the midpoint between failures times
  midfail_times <- rollmean(c(0, fail_times), 2)
  midfail_times_pretaub <- rollmean(c(0, fail_times_pretaub), 2)
  midfail_times_postaub <- rollmean(c(0, fail_times_postaub), 2)

  # Number of distinct ordered observed failures times
  l = length(fail_times)
  l_pretaub = length(fail_times_pretaub)
  l_postaub = length(fail_times_postaub)


  ############################################################################################
  # Compute the covariate integral
  # INTEGRAL PART 1: from tau0 to taub
  ############################################################################################

  # ESTIMATES' VALUES FOR INTEGRALS' COMPUTATION
  ################################################
  # Note that km0_f(midfail_times_pretaub) corresponds to Kaplan-Meier estimates' values at midpoints between event times (failure and censoring times)
  # Note that km0_f(fail_times_pretaub) corresponds to Kaplan-Meier estimates' values at event times (failure and censoring times)
  # Analogously with km1_f(), kmp_f(), censkm0_f(), censkm1_f(), cpn1_f(), cpn0_f(), ...

  # Counting process: dN function

  # version with kernel
  hazard1_values <- hazard1_f(fail_times_pretaub)
  hazard0_values <- hazard0_f(fail_times_pretaub)
  hazard1X_values <- hazard1X_f(fail_times_pretaub)
  hazard0X_values <- hazard0X_f(fail_times_pretaub)

  # Compute the weight function
  ######################################

  # weight function proposed by Pepe-Fleming: w
  w <- ifelse(censkm0_f(midfail_times_pretaub)+censkm1_f(midfail_times_pretaub) == 0,
              0,
              (n*censkm0_f(midfail_times_pretaub)*censkm1_f(midfail_times_pretaub))/(n0*censkm0_f(midfail_times_pretaub)+n1*censkm1_f(midfail_times_pretaub))
  )
  # weight function proposed by Fleming-Harrington: f^rho*(1-f)^gam
  f <- kmp_f(midfail_times_pretaub)
  # Define the weight function as a product of w and f
  weight <- w^(eta)*f^rho*(1-f)^gam

  # Kaplan-Meier jumps
  # group 0
  KM0_jumps <- diff(c(1,km0_f(fail_times_pretaub)))

  # group 1
  KM1_jumps <- diff(c(1,km1_f(fail_times_pretaub)))

  # Calculate the integral: int_{t,tau}(weight * surv)
  # group 0
  Integral0 <- cumsum(diff(c(0, fail_times_pretaub)) * weight * km0_f(fail_times_pretaub))
  Int0 <- (Integral0[l_pretaub] - Integral0)
  # group 1
  Integral1 <- cumsum(diff(c(0, fail_times_pretaub)) * weight * km1_f(fail_times_pretaub))
  Int1 <- (Integral1[l_pretaub] - Integral1)


  # Calculate the integral 1A
  integral1A = -(n1/n)*sum(Int0*hazard0X_values*diff(c(0, fail_times_pretaub)),na.rm = TRUE) - (n0/n)*sum(Int1*hazard1X_values*diff(c(0, fail_times_pretaub)),na.rm = TRUE)

  # Calculate the integral 1B
  integral1B = -(n1/n)*sum(Int0*phat_group0*KM0_jumps/km0_f(fail_times_pretaub),na.rm = TRUE) - (n0/n)*sum(Int1*phat_group1*KM1_jumps/km1_f(fail_times_pretaub),na.rm = TRUE)

  # Calculate the integral part 1
  cov_value_pretaub =  integral1A + integral1B

  ############################################################################################
  # Compute the covariate integral
  # INTEGRAL 2: from taub to tau
  ############################################################################################
  cov_value_posttau = 0
  if(taub<tau){
    # Compute the weight function
    ######################################
    # weight function proposed by Pepe-Fleming
    w <- ifelse(censkm0_f(midfail_times_postaub)+censkm1_f(midfail_times_postaub) == 0,
                0,
                (n*censkm0_f(midfail_times_postaub)*censkm1_f(midfail_times_postaub))/(n0*censkm0_f(midfail_times_postaub) + n1*censkm1_f(midfail_times_postaub)))
    # weight function proposed by Fleming-Harrington
    f <- kmp_f(midfail_times_postaub)
    # Define the weight function as a product of w and f
    weight <- w^(eta)*f^rho*(1-f)^gam

    # Compute the Kaplan-Meier jumps
    ######################################
    # Kaplan-Meier jumps
    # group 0
    KM0_jumps <- diff(c(1,km0_f(fail_times_postaub)))
    # group 1
    KM1_jumps <- diff(c(1,km1_f(fail_times_postaub)))

    # Kaplan-Meier responders
    # group 0
    KM0x_jumps <- diff(c(1,km0x_f(fail_times_postaub)))
    # group 1
    KM1x_jumps <- diff(c(1,km1x_f(fail_times_postaub)))

    # Calculate the integral: int_{t,tau}(weight * surv)
    # group 0
    Integral0 <- cumsum(diff(c(0, fail_times_postaub)) * weight * km0_f(fail_times_postaub))
    Int0 <- (Integral0[l_postaub] - Integral0)
    # group 1
    Integral1 <- cumsum(diff(c(0, fail_times_postaub)) * weight * km1_f(fail_times_postaub))
    Int1 <- (Integral1[l_postaub] - Integral1)

    # Calculate the integral part 1
    sum_part1 = ifelse(km0x_f(fail_times_postaub)*km0_f(fail_times_postaub)*km0_f(midfail_times_postaub) == 0,
                       0,
                       Int0*phat_group0*(KM0x_jumps*km0x_f(midfail_times_postaub)/km0x_f(fail_times_postaub)-
                                           KM0_jumps*km0_f(midfail_times_postaub)/km0_f(fail_times_postaub))/km0_f(midfail_times_postaub)
    )

    sum_part2 = ifelse(km1x_f(fail_times_postaub)*km1_f(fail_times_postaub)*km1_f(midfail_times_postaub) == 0,
                       0,
                       Int1*phat_group1*(KM1x_jumps*km1x_f(midfail_times_postaub)/km1x_f(fail_times_postaub) -
                                           KM1_jumps*km1_f(midfail_times_postaub)/km1_f(fail_times_postaub))/km1_f(midfail_times_postaub)
    )

    cov_value_posttau =  (n1*sum(sum_part1) + n0*sum(sum_part2))/n

  }

  ############################################################################################
  # Compute the covariance
  ############################################################################################

  # Calculate the estimated (pooled) covariance
  cov_value =  cov_value_pretaub + cov_value_posttau

  return(cov_value)
}

#######################################################


