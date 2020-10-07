#######################################################################################

# simsurvbin_H1_nonPH : simulates binary and survival data under H1 
# and with delayed effects on survival; returns the database 

#######################################################################################


##################################################################################
simsurvbin_H1_nonPH <- function(a.shape, b.scale, HR, rate.param, p0, p1, ass.par, n0, n1, censoring="Exp", tstar){

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
  cp <- claytonCopula(param = c(0.5), dim = 2)
  # cp <- frankCopula(param = c(2), dim = 2)

  # Generate the multivariate distribution
  copulaSB <- mvdc(copula = cp,
                   margins = c("unif", "unif"),
                   paramMargins = list(list(0,1), list(0,1)))
  

  copulaSB_delbefore <- mvdc(copula = cp,
                             margins = c("unif", "unif"),
                             paramMargins = list(list(min=exp(-lambda*tstar^(1/rho)), max=1),
                                                 list(0,1)))
  
  copulaSB_delbeyond <- mvdc(copula = cp,
                             margins = c("unif", "unif"),
                             paramMargins = list(list(min=0, max=exp(-lambda*tstar^(1/rho))),
                                                 list(0,1)))
  
    

  # CONTROL GROUP
  ######################################
  v = rMvdc(n1,copulaSB)
  BE0 = ifelse(v[,2]<p0, 1, 0)
  TE0 = (-log(v[,1])/(lambda))^(1/rho)

  time0= ifelse(TE0<=TC0, TE0, TC0) 
  status0 = ifelse(TE0<=TC0,1,0)
  treat0 = rep(0,n0)

  # TREATMENT GROUP
  ######################################
  
  p = pweibull(q=tstar, shape=a.shape, scale = b.scale)
  
  # before tstar 
  v = rMvdc(n=round(n1*p),copulaSB_delbefore)
  
  BE1_bef <-  ifelse(v[,2]<p1, 1, 0)
  TE1_bef <- (- log(v[,1])/(lambda))^(1/rho)
  
  # beyond tstar  
  v = rMvdc(n=n1-round(n1*p),copulaSB_delbeyond)
  
  BE1_bey = ifelse(v[,2]<p0, 1, 0)
  TE1_bey <- (- log(v[,1])/(lambda*HR))^(1/rho) 
  
  BE1=c(BE1_bef,BE1_bey)
  TE1=c(TE1_bef,TE1_bey)

  time1 = ifelse(TE1<=TC1,TE1,TC1) 
  status1 = ifelse(TE1<=TC1,1,0)
  treat1 = rep(1,n1)

  # TWO-SAMPLE DATA
  ######################################
  db = data.frame(binary=c(BE1, BE0), time=c(time1,time0), status=c(status1,status0),treat=c(treat1,treat0))

  return(db)
}

# 
# # OLD version
# if(H0==TRUE){
#   v = rMvdc(n1,copulaSB)
#   
#   BE1 = ifelse(v[,2]<p0, 1, 0)
#   TE1 = (-log(v[,1])/(lambda))^(1/rho)
#   
# }else{
#   p = pweibull(q=tstar, shape=a.shape, scale = b.scale)
#   
#   # before tstar 
#   v = rMvdc(n=round(n1*p),copulaSB_delbefore)
#   
#   BE1_bef <-  ifelse(v[,2]<p1, 1, 0)
#   TE1_bef <- (- log(v[,1])/(lambda*HR))^(1/rho)
#   
#   # beyond tstar 
#   # v = rMvdc(n=round(n1*(1-p)),copulaSB_delbeyond)
#   v = rMvdc(n=n1-round(n1*p),copulaSB_delbeyond)
#   
#   BE1_bey = ifelse(v[,2]<p0, 1, 0)
#   TE1_bey = (-log(v[,1])/(lambda))^(1/rho)
#   
#   BE1=c(BE1_bef,BE1_bey)
#   TE1=c(TE1_bef,TE1_bey)
# }
