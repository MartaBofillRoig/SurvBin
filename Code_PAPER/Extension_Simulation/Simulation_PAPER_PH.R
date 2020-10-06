
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())
# Working directory
# setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation")
setwd("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation")

#####################################################################################
# PREAMBLE
#####################################################################################

library(copula)
library(pracma) # for using the function 'integral'
library(survival)
library(zoo) # 'rollmean' function
library(muhaz)
library(tidyr)
library(purrr) # 'possibly' function

#####################################################################################

# Functions for the binary and survival setting; for the covariance computation; and for simulating the binary and time-to-event data
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/binary-functions.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/survival-functions.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/cov-functions.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/simulation-functions.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/simulation-functions_H1.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/lstats-functions.R')
# source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/lstats_boots.R')

source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/binary-functions.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/survival-functions.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/cov-functions.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/simulation-functions.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/simulation-functions_H1.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/lstats-functions.R')
source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/lstats_boots.R')

# Parameters
alpha=0.05;
z.alpha <- qnorm(1-alpha,0,1)
z.alphac <- qnorm(1-alpha/2,0,1)

# nsim: number of simulations
nsim=100

# Note:
# alpha=0.05
# > nsim=500
# > sd=sqrt(alpha*(1-alpha)/nsim)
# > z.alpha <- qnorm(1-alpha,0,1)
# > c(alpha-z.alpha*sd,alpha+z.alpha*sd)
# [1] 0.03396795 0.06603205
# > nsim=1000
# > sd=sqrt(alpha*(1-alpha)/nsim)
# > z.alpha <- qnorm(1-alpha,0,1)
# > c(alpha-z.alpha*sd,alpha+z.alpha*sd)
# [1] 0.03866363 0.06133637

# simulation seed
set.seed(1939)

#####################################################################################
# Scenarios H0 FALSE -- case 1
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=c(1)
taub= c(0.5,1)
r=c(3) #unif
p0=c(0.1, 0.3)
d=c(0.1)
HR=c(0.75)
theta=c(0.001,2,3)
n=c(500)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
# omegab=c(0.5)
omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,HR,tau,taub,r,p0,d,theta,n,eta,rho,gamma,omegab)

data$p1=data$d+data$p0
data$omegas=1-data$omegab 
rm(a,b,taub,tau,r,p0,d,HR,omegab,theta,n,eta,rho,gamma) 

#####################################################################################
#####################################################################################

t0=Sys.time()
data$Test_Power_pluginU=0
data$Test_Power_pluginP=0
data$Test_Power_Boots=0
data$Test_Power_Bonf=0
data$Test_Power_S=0
data$Test_Power_B=0

#####################################################################################
# H0 FALSE: Empirical powers
#####################################################################################


# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){

  data$Test_Power_pluginU[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Unpooled')) > z.alpha)/nsim

  data$Test_Power_pluginP[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Pooled')) > z.alpha)/nsim

  data$Test_Power_Boots[i] <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                   rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                   ass.par=data$theta[i],
                                                                   n0=data$n[i]/2, n1=data$n[i]/2,
                                                                   censoring="Unif",
                                                                   tau=data$tau[i],
                                                                   taub=data$taub[i],
                                                                   rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                                   wb=data$omegab[i],
                                                                   ws=data$omegas[i],
                                                                   Boot=50)) > z.alpha)/nsim

  data$Test_Power_Bonf[i] <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                      rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                      ass.par=data$theta[i],
                                                                      n0=data$n[i]/2, n1=data$n[i]/2,
                                                                      censoring="Unif",
                                                                      tau=data$tau[i],
                                                                      taub=data$taub[i],
                                                                      rho=data$rho[i], gam=data$gamma[i],
                                                                      eta=data$eta[i])>z.alphac)>=1)))/nsim

  data$Test_Power_S[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[2]) > z.alpha)/nsim

  data$Test_Power_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim

  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_H1_case1.txt", append=TRUE)
  save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case1.RData")
  # save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case1.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case1.RData")
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")


#####################################################################################
# Scenarios H0 FALSE -- case 2
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=c(1)
taub= c(0.5,1)
r=c(3) #unif
p0=c(0.1, 0.3)
d=c(0)
HR=c(0.75)
theta=c(0.001,2,3)
n=c(500)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
# omegab=c(0.5)
omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,HR,tau,taub,r,p0,d,theta,n,eta,rho,gamma,omegab)

data$p1=data$d+data$p0
data$omegas=1-data$omegab

rm(a,b,taub,tau,r,p0,d,HR,omegab,theta,n,eta,rho,gamma)


#####################################################################################
##################################################################################### 

t0=Sys.time()
data$Test_Power_pluginU=0
data$Test_Power_pluginP=0
data$Test_Power_Boots=0
data$Test_Power_Bonf=0
data$Test_Power_S=0
data$Test_Power_B=0

#####################################################################################
# H0 FALSE: Empirical powers
#####################################################################################


# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){
  
  data$Test_Power_pluginU[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Unpooled')) > z.alpha)/nsim
  
  data$Test_Power_pluginP[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Pooled')) > z.alpha)/nsim
  
  data$Test_Power_Boots[i] <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                   rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                   ass.par=data$theta[i],
                                                                   n0=data$n[i]/2, n1=data$n[i]/2,
                                                                   censoring="Unif",
                                                                   tau=data$tau[i],
                                                                   taub=data$taub[i],
                                                                   rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                                   wb=data$omegab[i],
                                                                   ws=data$omegas[i],
                                                                   Boot=50)) > z.alpha)/nsim
  
  data$Test_Power_Bonf[i] <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                      rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                      ass.par=data$theta[i],
                                                                      n0=data$n[i]/2, n1=data$n[i]/2,
                                                                      censoring="Unif",
                                                                      tau=data$tau[i],
                                                                      taub=data$taub[i],
                                                                      rho=data$rho[i], gam=data$gamma[i],
                                                                      eta=data$eta[i])>z.alphac)>=1)))/nsim
  
  data$Test_Power_S[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[2]) > z.alpha)/nsim
  
  data$Test_Power_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim
  
  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_H1_case2.txt", append=TRUE)
  save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case2.RData")
  # save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case2.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case2.RData")
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")


#####################################################################################
# Scenarios H0 FALSE -- case 3
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=c(1)
taub= c(0.5,1)
r=c(3) #unif
p0=c(0.1, 0.3)
d=c(0.1)
HR=c(1)
theta=c(0.001,2,3)
n=c(500)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
# omegab=c(0.5)
omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,HR,tau,taub,r,p0,d,theta,n,eta,rho,gamma,omegab)

data$p1=data$d+data$p0
data$omegas=1-data$omegab

rm(a,b,taub,tau,r,p0,d,HR,omegab,theta,n,eta,rho,gamma)


#####################################################################################
##################################################################################### 

t0=Sys.time()
data$Test_Power_pluginU=0
data$Test_Power_pluginP=0
data$Test_Power_Boots=0
data$Test_Power_Bonf=0
data$Test_Power_S=0
data$Test_Power_B=0

#####################################################################################
# H0 FALSE: Empirical powers
#####################################################################################


# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){
  
  data$Test_Power_pluginU[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Unpooled')) > z.alpha)/nsim
  
  data$Test_Power_pluginP[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Pooled')) > z.alpha)/nsim
  
  data$Test_Power_Boots[i] <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                   rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                   ass.par=data$theta[i],
                                                                   n0=data$n[i]/2, n1=data$n[i]/2,
                                                                   censoring="Unif",
                                                                   tau=data$tau[i],
                                                                   taub=data$taub[i],
                                                                   rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                                   wb=data$omegab[i],
                                                                   ws=data$omegas[i],
                                                                   Boot=50)) > z.alpha)/nsim
  
  data$Test_Power_Bonf[i] <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                                      rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                                      ass.par=data$theta[i],
                                                                      n0=data$n[i]/2, n1=data$n[i]/2,
                                                                      censoring="Unif",
                                                                      tau=data$tau[i],
                                                                      taub=data$taub[i],
                                                                      rho=data$rho[i], gam=data$gamma[i],
                                                                      eta=data$eta[i])>z.alphac)>=1)))/nsim
  
  data$Test_Power_S[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[2]) > z.alpha)/nsim
  
  data$Test_Power_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim
  
  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_H1_case3.txt", append=TRUE)
  save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case3.RData")
  # save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case3.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1_case3.RData")
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H1.RData")




#####################################################################################
# H0 TRUE: Empirical sign level
#####################################################################################


#####################################################################################
# Scenarios H0 TRUE
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=c(1)
taub= c(0.5,1)
r=c(3) #unif
p0=c(0.1, 0.3)
d=c(0)
HR=c(1)
theta=c(0.001,2,3)
n=c(500)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
omegab=c(0.5)
# omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,HR,tau,taub,r,p0,d,theta,n,eta,rho,gamma,omegab)

data$p1=data$d+data$p0
data$omegas=1-data$omegab

# data = subset(data,data$eta+data$rho+data$gamma>0) #eta=1


rm(a,b,taub,tau,r,p0,d,HR,omegab,theta,n,eta,rho,gamma)
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/database_Scenarios_H0_True.RData")
save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/database_Scenarios_H0_True.RData")


# nsim: number of simulations
nsim=1000

#####################################################################################
# simulation seed
set.seed(1987)

t0=Sys.time()
data$Test_Alpha_pluginU=0
data$Test_Alpha_pluginP=0
data$Test_Alpha_Boots=0
data$Test_Alpha_Bonf=0
data$Test_Alpha_S=0
data$Test_Alpha_B=0

# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){
  
  data$Test_Alpha_pluginU[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Unpooled')) > z.alpha)/nsim
  
  data$Test_Alpha_pluginP[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Pooled')) > z.alpha)/nsim
  
  data$Test_Alpha_Boots[i] <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                                   rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                                   ass.par=data$theta[i],
                                                                   n0=data$n[i]/2, n1=data$n[i]/2,
                                                                   censoring="Unif",
                                                                   tau=data$tau[i],
                                                                   taub=data$taub[i],
                                                                   rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                                   wb=data$omegab[i],
                                                                   ws=data$omegas[i],
                                                                   Boot=50)) > z.alpha)/nsim
  
  data$Test_Alpha_Bonf[i] <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                                      rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                                      ass.par=data$theta[i],
                                                                      n0=data$n[i]/2, n1=data$n[i]/2,
                                                                      censoring="Unif",
                                                                      tau=data$tau[i],
                                                                      taub=data$taub[i],
                                                                      rho=data$rho[i], gam=data$gamma[i],
                                                                      eta=data$eta[i])>z.alphac)>=1)))/nsim
  
  data$Test_Alpha_S[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[2]) > z.alpha)/nsim
  
  data$Test_Alpha_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim
  
  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_H0.txt", append=TRUE)
  save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0.RData")
  # save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H0.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0.RData")
# save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0.RData")

#####################################################################################
# results empirical alpha without correcting small sample size 

source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/cov-functions_wc.R')

data$Test_Alpha_pluginU_c=0
data$Test_Alpha_pluginP_c=0 


set.seed(1983)

# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){
  
  data$Test_Alpha_pluginU_c[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Unpooled')) > z.alpha)/nsim
  
  data$Test_Alpha_pluginP_c[i] <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], HR=1,
                                                               rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                               ass.par=data$theta[i],
                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                               censoring="Unif",
                                                               tau=data$tau[i],
                                                               taub=data$taub[i],
                                                               rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                               wb=data$omegab[i], ws=data$omegas[i],
                                                               var_est='Pooled')) > z.alpha)/nsim  
  
  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_H0.txt", append=TRUE)
  save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0_wc.RData")
  # save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0.RData")
}

save.image("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_H0_wc.RData")


