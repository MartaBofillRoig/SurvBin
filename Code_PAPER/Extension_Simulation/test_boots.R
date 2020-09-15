
#####################################################################################

lstats_boots <- function(x1, x2, x3, x4, tau0=0, tau=NULL, taub=NULL, rho=0, gam=0, eta=1, w1=0.5, w2=0.5){ 
  
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
  u_bs = w1*test_b + w2*test_s
  
  # bootstrap
  # B=1000
  # ubs_boot <- data.frame(0)
  # for(i in 1:B){
  #   daux <- db[sample(1:nrow(db),nrow(db),rep=TRUE),]  
  #   btt <- bintest(daux$binary, daux$treat)[1] 
  #   stt <- survtest(daux$time, daux$status, daux$treat)[1]
  #   ubs_boot[i] = w1*btt + w2*stt
  # }
  # stdev <- sd(ubs_boot)  
  # test_l <- u_bs/stdev
  # TestBS = lstats(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws, var_est)
  
  cbind.data.frame(time=time, status=status, binary=binary, treat=treat)
  
  data=cbind.data.frame(time=db$time, status=db$status, binary=db$binary, treat=db$treat)
  boot(db, lstats_bnames, R=1000, tau0=0, tau=tau, taub=taub, rho=rho, gam=gam, eta=eta, wb=wb, ws=ws, var_est=var_est)
  
  myBootstrap <- boot(db, lstats, R=1000, tau0=0, tau=tau, taub=taub, rho=rho, gam=gam, eta=eta, wb=wb, ws=ws, var_est=var_est)
  
  return(list=c(LTest=test_l,Statistic=u_bs,sd=stdev)) 
}

#####################################################################################

lstats_bnames <- function(db, tau0=0, tau=NULL, taub=NULL, rho=0, gam=0, eta=1, wb=0.5, ws=0.5, var_est="Unpooled"){
  # time, status, binary, treat
  
  names(db)[1] <- "time"
  names(db)[2] <- "status"
  names(db)[3] <- "binary"
  names(db)[4] <- "treat"
  
  # db=cbind.data.frame(time=time, status=status, binary=binary, treat=treat)
  
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
  
  return(list=c(LTest=test_l,Statistic=u_bs,sd=stdev))
}
##################################################################################

fCS.TEST_boots <- function(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring="Exp", tau, taub, rho, gam, eta, wb, ws){
  
  # TWO-SAMPLE db
  ######################################
  db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring)
  
  # STATISTICS
  ######################################
  TestBS = lstats_boots(db$time,db$status, db$binary, db$treat, tau0=0, tau, taub, rho, gam, eta, wb, ws)
  
  return(TestBS[1])
}
