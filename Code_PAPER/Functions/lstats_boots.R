
#####################################################################################

lstats_boots <- function(x1, x2, x3, x4, tau0=0, tau=NULL, taub=NULL, rho=0, gam=0, eta=1, wb=0.5, ws=0.5, Boot=50){

  db=cbind.data.frame(time=x1, status=x2, binary=x3, treat=x4)

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

  return(list=c(LTest=test_l,Statistic=u_bs,sd=stdev))
}

#####################################################################################

##################################################################################

fCS.TEST_boots <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, Boot){

  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring, copula)

  # STATISTICS
  ######################################
  TestBS = lstats_boots(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws, Boot)

  return(TestBS[1])
}


##################################################################################

fCS.TEST_boots_H1 <- function(a.shape, b.scale, HR, rate.param, p0, p1, ass.par,  n0, n1, censoring="Exp", copula="clayton", tau, taub, rho, gam, eta, wb, ws, Boot, PH=TRUE, tstar=0){

  # TWO-SAMPLE db
  ######################################
  if(PH==TRUE){
    db = simsurvbin_H1(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, H0=FALSE) 
  }else{
    db = simsurvbin_H1_nonPH(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring, copula, tstar) 
  } 

  # STATISTICS
  ######################################
  TestBS = lstats_boots(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws, Boot)

  return(TestBS[1])
}
