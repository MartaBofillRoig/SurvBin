#######################################################################################

# simsurvbin : simulates binary and survival data; returns the database
# fCS.TEST : simulates binary and survival data; returns the L-statistic
# fCS.TEST_Bonf : simulates binary and survival data; returns the surv and bin tests

####################################################################################### 
##################################################################################

simsurvbin <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton"){

  # CENSORING
  ######################################
  if(censoring=="Exp"){
    TC1 = rexp(n= ss, rate = rate.param)
    TC0 = rexp(n= ss, rate = rate.param)
  }
  if(censoring=="Unif"){
    TC1 = runif(n= ss, min=0,max=rate.param)
    TC0 = runif(n= ss, min=0,max=rate.param)
  }

  # RE-PARAMETRIZATION for time-to-event simulation
  ######################################
  lambda = (1/b.scale)^a.shape
  rho = a.shape
  
  # Copula for binary and time-to-event
  ######################################
  # Select the copula
  if(copula=="clayton"){
    cp <- claytonCopula(param = c(ass.par), dim = 2)
  }
  if(copula=="frank"){
    cp <- frankCopula(param = c(ass.par), dim = 2)
  } 
  
  # Generate the multivariate distribution
  copulaSB <- mvdc(copula = cp,
                   margins = c("unif", "unif"),
                   paramMargins = list(list(0,1), list(0,1)))
  

  # TREATMENT GROUP
  ######################################
  v = rMvdc(n1,copulaSB)

  # without latent variable
  BE1 = ifelse(v[,2]<prob0, 1, 0)

  # time-to-event (survival)
  TE1 = b.scale*(-log(1-v[,1]))^(1/a.shape)
  time1= ifelse(TE1<=TC1, TE1,TC1)
  status1 = ifelse(TE1<=TC1,1,0)
  treat1 = rep(1,ss)

  # CONTROL GROUP
  ######################################
  v = rMvdc(n0,copulaSB)

  # without latent variable
  BE0 = ifelse(v[,2]<prob0, 1, 0)

  # time-to-event (survival)
  TE0 = b.scale*(-log(1-v[,1]))^(1/a.shape)
  time0= ifelse(TE0<=TC0, TE0,TC0)
  status0 = ifelse(TE0<=TC0,1,0)
  treat0 = rep(0,ss)

  # TWO-SAMPLE db
  ######################################
  treat0=as.vector(treat0)
  treat1=as.vector(treat1)
  db = data.frame(binary=c(BE1, BE0), time=c(time1,time0), status=c(status1,status0),treat=c(treat1,treat0))

  return(db)
}


##################################################################################
simsurvbin_H1 <- function(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring="Exp", copula="clayton", H0=FALSE){

  # CENSORING
  ######################################
  if(censoring=="Exp"){
    TC1 = rexp(n= n1, rate = rate.param)
    TC0 = rexp(n= n0, rate = rate.param)
  }
  if(censoring=="Unif"){
    TC1 = runif(n= n1, min=0,max=rate.param)
    TC0 = runif(n= n0, min=0,max=rate.param)
  }

  # RE-PARAMETRIZATION for time-to-event simulation
  ######################################
  lambda=(1/b.scale)^a.shape
  rho=a.shape

  # Copula for binary and time-to-event
  ######################################
  # Select the copula
  if(copula=="clayton"){
    cp <- claytonCopula(param = c(ass.par), dim = 2)
  }
  if(copula=="frank"){
    cp <- frankCopula(param = c(ass.par), dim = 2)
  } 

  # Generate the multivariate distribution
  copulaSB <- mvdc(copula = cp,
                   margins = c("unif", "unif"),
                   paramMargins = list(list(0,1), list(0,1)))

  # CONTROL GROUP
  ######################################
  
  v = rMvdc(n0,copulaSB)
  BE0 = ifelse(v[,2]<p0, 1, 0)
  TE0 = (-log(v[,1])/(lambda))^(1/rho)

  time0= ifelse(TE0<=TC0, TE0, TC0) 
  status0 = ifelse(TE0<=TC0,1,0)
  treat0 = rep(0,n0)

  # TREATMENT GROUP
  ######################################

  v = rMvdc(n1,copulaSB)
  if(H0==TRUE){

    BE1 = ifelse(v[,2]<p0, 1, 0)
    TE1 = (-log(v[,1])/(lambda))^(1/rho)

  }else{

    BE1 <-  ifelse(v[,2]<p1, 1, 0)
    TE1 <- (- log(v[,1])/(lambda*HR))^(1/rho)
  }

  time1 = ifelse(TE1<=TC1,TE1,TC1) 
  status1 = ifelse(TE1<=TC1,1,0)
  treat1 = rep(1,n1)

  # TWO-SAMPLE DATA
  ######################################
  db = data.frame(binary=c(BE1, BE0), time=c(time1,time0), status=c(status1,status0),treat=c(treat1,treat0))


  return(db)
}


##################################################################################

fCS.TEST <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, var_est){

  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring, copula)

  # STATISTICS
  ######################################
  TestBS = lstats(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws, var_est)

  return(TestBS[1])
}

##################################################################################

fCS.TEST_H1 <- function(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, var_est, PH=TRUE, tstar=0){

  # TWO-SAMPLE db
  ######################################
  if(PH==TRUE){
    db = simsurvbin_H1(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, H0=FALSE) 
  }else{
    db = simsurvbin_H1_nonPH(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, tstar) 
  }
  

  # STATISTICS
  ######################################
  TestBS = lstats(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws, var_est)

  return(TestBS[1])
}



##################################################################################

fCS.TEST_Bonf <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta){

  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring, copula)

  # STATISTICS
  ######################################
  B <- bintest(db$binary, db$treat, var_est="Unpooled")
  test_b <- B[1]

  S <- survtest(db$time, db$status, db$treat, tau, rho, gam, eta)
  test_s <- S[1]

  return(c(test_b,test_s))
}


##################################################################################

fCS.TEST_Bonf_H1 <- function(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, PH=TRUE, tstar=0){

  # TWO-SAMPLE db
  ######################################
  if(PH==TRUE){
    db = simsurvbin_H1(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, H0=FALSE) 
  }else{
    db = simsurvbin_H1_nonPH(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, tstar) 
  }
  
  
  # STATISTICS
  ######################################
  B <- bintest(db$binary, db$treat, var_est="Unpooled")
  test_b <- B[1]

  S <- survtest(db$time, db$status, db$treat, tau, rho, gam, eta)
  test_s <- S[1]

  return(c(test_b,test_s))
}


##################################################################################

fCS.TEST_s <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, var_est='Unpooled'){

  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring, copula)

  # STATISTICS
  ######################################
  S <- survtest(db$time, db$status, db$treat, tau, rho, gam, eta, var_est)

  return(S[1])
}


##################################################################################

fCS.TEST_b <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, var_est="Unpooled"){

  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring, copula)

  # STATISTICS
  ######################################
  B <- bintest(db$binary, db$treat, var_est)

  return(B[1])
}
