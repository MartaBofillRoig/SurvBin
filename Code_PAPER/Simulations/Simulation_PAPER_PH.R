
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())
# Working directory 
setwd("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation")

#####################################################################################
# PREAMBLE
#####################################################################################

# install.packages(c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))

library(copula)
library(pracma) # for using the function 'integral'
library(survival)
library(zoo) # 'rollmean' function
library(muhaz)
library(tidyr)
library(purrr) # 'possibly' function
library(foreach)
library(doParallel)

#####################################################################################

# Functions for the binary and survival setting; for the covariance computation; and for simulating the binary and time-to-event data
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/binary-functions.R')
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/survival-functions.R')
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/cov-functions.R')
# source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/simulation-functions.R')
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/simulation-functions_H1.R')
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/lstats-functions.R')
source('C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Functions/lstats_boots.R')

#####################################################################################

# setup parallel backend to use many processors
cores=detectCores()
cl <- makeCluster(cores[1]-1) #not to overload your computer
registerDoParallel(cl) 

#####################################################################################

# Parameters
alpha=0.05;
z.alpha <- qnorm(1-alpha,0,1)
z.alphac <- qnorm(1-alpha/2,0,1)

# nsim: number of simulations 
nsim=1000

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
d=c(0.075)
HR=c(0.75)
n=c(500)

# Assume we consider the following values for the rho spearman
# iRho(clayton.cop,rho=0.001)
# iRho(clayton.cop,rho=0.3)
# iRho(clayton.cop,rho=0.45)
theta=c(0.001,0.51,0.91)
# for the frank copula
# theta=c(0.001,2,3)

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


##################################################################################### 
#####################################################################################

Test_Power_pluginU <- c()
Test_Power_pluginU <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                              
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                           HR=data$HR[i],
                                                                           rate.param=data$r[i], 
                                                                           p0=data$p0[i], p1=data$p1[i],
                                                                           ass.par=data$theta[i],
                                                                           n0=data$n[i]/2, n1=data$n[i]/2,
                                                                           censoring="Unif",
                                                                           tau=data$tau[i],
                                                                           taub=data$taub[i],
                                                                           rho=data$rho[i], gam=data$gamma[i], 
                                                                           eta=data$eta[i],
                                                                           wb=data$omegab[i], ws=data$omegas[i],
                                                                           var_est='Unpooled')) > z.alpha)/nsim
                                tempMatrix 
}

x=t(Test_Power_pluginU)
data$Test_Power_pluginU <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_pluginu.RData")

#####################################################################################

Test_Power_pluginP <- c()
Test_Power_pluginP <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=data$HR[i],
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p1[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Pooled')) > z.alpha)/nsim
                                
                                tempMatrix
                                 
                              }

x=t(Test_Power_pluginP)
data$Test_Power_pluginP <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_plugin.RData")

#####################################################################################

Test_Power_Boots <- c()
Test_Power_Boots <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                   HR=data$HR[i],
                                                                                   rate.param=data$r[i], 
                                                                                   p0=data$p0[i], p1=data$p1[i],
                                                                                   ass.par=data$theta[i],
                                                                                   n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                   censoring="Unif",
                                                                                   tau=data$tau[i],
                                                                                   taub=data$taub[i],
                                                                                   rho=data$rho[i], gam=data$gamma[i], 
                                                                                   eta=data$eta[i],
                                                                                   wb=data$omegab[i],
                                                                                   ws=data$omegas[i],
                                                                                   Boot=50)) > z.alpha)/nsim
                                
                                tempMatrix
                                
                              }

x=t(Test_Power_Boots)
data$Test_Power_Boots <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_boots.RData")

#####################################################################################

Test_Power_Bonf <- c()
Test_Power_Bonf <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                            .combine=cbind,
                            .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                              
                              tempMatrix <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                     HR=data$HR[i],
                                                                                     rate.param=data$r[i], 
                                                                                     p0=data$p0[i], p1=data$p1[i],
                                                                                     ass.par=data$theta[i],
                                                                                     n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                     censoring="Unif",
                                                                                     tau=data$tau[i],
                                                                                     taub=data$taub[i],
                                                                                     rho=data$rho[i], gam=data$gamma[i],
                                                                                     eta=data$eta[i])>z.alphac)>=1)))/nsim
                              
                              tempMatrix
                              
                            }

x=t(Test_Power_Bonf)
data$Test_Power_Bonf <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_bonf.RData")

#####################################################################################

Test_Power_S <- c()
Test_Power_S <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                           .combine=cbind,
                           .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                             
                             tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                               HR=data$HR[i],
                                                                               rate.param=data$r[i], 
                                                                               p0=data$p0[i], p1=data$p1[i],
                                                                               ass.par=data$theta[i],
                                                                               n0=data$n[i]/2, n1=data$n[i]/2,
                                                                               censoring="Unif",
                                                                               tau=data$tau[i],
                                                                               taub=data$taub[i],
                                                                               rho=data$rho[i], gam=data$gamma[i],
                                                                               eta=data$eta[i])[2]) > z.alpha)/nsim
                             
                             tempMatrix
                             
                           }

x=t(Test_Power_S)
data$Test_Power_S <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_surv.RData")

#####################################################################################

Test_Power_B <- c()
Test_Power_B <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=data$HR[i],
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p1[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[1]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Power_B)
data$Test_Power_B <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_bin.RData")

#####################################################################################


t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case1.txt", append=TRUE)
(t1) 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_case1.RData")



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
n=c(500)

# Assume we consider the following values for the rho spearman
# iRho(clayton.cop,rho=0.001)
# iRho(clayton.cop,rho=0.3)
# iRho(clayton.cop,rho=0.45)
theta=c(0.001,0.51,0.91)
# for the frank copula
# theta=c(0.001,2,3)

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

##################################################################################### 
#####################################################################################

Test_Power_pluginU <- c()
Test_Power_pluginU <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=data$HR[i],
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p1[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Unpooled')) > z.alpha)/nsim
                                tempMatrix 
                              }

x=t(Test_Power_pluginU)
data$Test_Power_pluginU <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_pluginu.RData")


#####################################################################################

Test_Power_pluginP <- c()
Test_Power_pluginP <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=data$HR[i],
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p1[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Pooled')) > z.alpha)/nsim
                                
                                tempMatrix
                                
                              }

x=t(Test_Power_pluginP)
data$Test_Power_pluginP <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_pluginp.RData")

#####################################################################################

Test_Power_Boots <- c()
Test_Power_Boots <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                            .combine=cbind,
                            .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                              
                              tempMatrix <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                 HR=data$HR[i],
                                                                                 rate.param=data$r[i], 
                                                                                 p0=data$p0[i], p1=data$p1[i],
                                                                                 ass.par=data$theta[i],
                                                                                 n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                 censoring="Unif",
                                                                                 tau=data$tau[i],
                                                                                 taub=data$taub[i],
                                                                                 rho=data$rho[i], gam=data$gamma[i], 
                                                                                 eta=data$eta[i],
                                                                                 wb=data$omegab[i],
                                                                                 ws=data$omegas[i],
                                                                                 Boot=50)) > z.alpha)/nsim
                              
                              tempMatrix
                              
                            }

x=t(Test_Power_Boots)
data$Test_Power_Boots <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_boot.RData")

#####################################################################################

Test_Power_Bonf <- c()
Test_Power_Bonf <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                           .combine=cbind,
                           .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                             
                             tempMatrix <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                    HR=data$HR[i],
                                                                                    rate.param=data$r[i], 
                                                                                    p0=data$p0[i], p1=data$p1[i],
                                                                                    ass.par=data$theta[i],
                                                                                    n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                    censoring="Unif",
                                                                                    tau=data$tau[i],
                                                                                    taub=data$taub[i],
                                                                                    rho=data$rho[i], gam=data$gamma[i],
                                                                                    eta=data$eta[i])>z.alphac)>=1)))/nsim
                             
                             tempMatrix
                             
                           }

x=t(Test_Power_Bonf)
data$Test_Power_Bonf <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_bonf.RData")

#####################################################################################

Test_Power_S <- c()
Test_Power_S <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=data$HR[i],
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p1[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[2]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Power_S)
data$Test_Power_S <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_surv.RData")


#####################################################################################

Test_Power_B <- c()
Test_Power_B <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=data$HR[i],
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p1[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[1]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Power_B)
data$Test_Power_B <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case2_bin.RData")

#####################################################################################

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case2.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_case2.RData")


#####################################################################################
# Scenarios H0 FALSE -- case 3
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=c(1)
taub= c(0.5,1)
r=c(3) #unif
p0=c(0.1, 0.3)
d=c(0.075)
HR=c(1) 
n=c(500)

# Assume we consider the following values for the rho spearman
# iRho(clayton.cop,rho=0.001)
# iRho(clayton.cop,rho=0.3)
# iRho(clayton.cop,rho=0.45)
theta=c(0.001,0.51,0.91)
# for the frank copula
# theta=c(0.001,2,3)

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

##################################################################################### 
#####################################################################################

Test_Power_pluginU <- c()
Test_Power_pluginU <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=data$HR[i],
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p1[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Unpooled')) > z.alpha)/nsim
                                tempMatrix 
                              }

x=t(Test_Power_pluginU)
data$Test_Power_pluginU <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_pluginu.RData")

#####################################################################################

Test_Power_pluginP <- c()
Test_Power_pluginP <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=data$HR[i],
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p1[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Pooled')) > z.alpha)/nsim
                                
                                tempMatrix
                                
                              }

x=t(Test_Power_pluginP)
data$Test_Power_pluginP <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_pluginp.RData")

#####################################################################################

Test_Power_Boots <- c()
Test_Power_Boots <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                            .combine=cbind,
                            .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                              
                              tempMatrix <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                 HR=data$HR[i],
                                                                                 rate.param=data$r[i], 
                                                                                 p0=data$p0[i], p1=data$p1[i],
                                                                                 ass.par=data$theta[i],
                                                                                 n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                 censoring="Unif",
                                                                                 tau=data$tau[i],
                                                                                 taub=data$taub[i],
                                                                                 rho=data$rho[i], gam=data$gamma[i], 
                                                                                 eta=data$eta[i],
                                                                                 wb=data$omegab[i],
                                                                                 ws=data$omegas[i],
                                                                                 Boot=50)) > z.alpha)/nsim
                              
                              tempMatrix
                              
                            }

x=t(Test_Power_Boots)
data$Test_Power_Boots <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_boot.RData")

#####################################################################################

Test_Power_Bonf <- c()
Test_Power_Bonf <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                           .combine=cbind,
                           .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                             
                             tempMatrix <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                    HR=data$HR[i],
                                                                                    rate.param=data$r[i], 
                                                                                    p0=data$p0[i], p1=data$p1[i],
                                                                                    ass.par=data$theta[i],
                                                                                    n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                    censoring="Unif",
                                                                                    tau=data$tau[i],
                                                                                    taub=data$taub[i],
                                                                                    rho=data$rho[i], gam=data$gamma[i],
                                                                                    eta=data$eta[i])>z.alphac)>=1)))/nsim
                             
                             tempMatrix
                             
                           }

x=t(Test_Power_Bonf)
data$Test_Power_Bonf <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_bonf.RData")

#####################################################################################

Test_Power_S <- c()
Test_Power_S <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=data$HR[i],
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p1[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[2]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Power_S)
data$Test_Power_S <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_surv.RData")

#####################################################################################

Test_Power_B <- c()
Test_Power_B <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=data$HR[i],
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p1[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[1]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Power_B)
data$Test_Power_B <- x

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1case3_bin.RData")

#####################################################################################

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H1_case3.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H1_case3.RData")



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
n=c(500)

# Assume we consider the following values for the rho spearman
# iRho(clayton.cop,rho=0.001)
# iRho(clayton.cop,rho=0.3)
# iRho(clayton.cop,rho=0.45)
theta=c(0.001,0.51,0.91)
# for the frank copula
# theta=c(0.001,2,3)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
omegab=c(0.5)
# omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,HR,tau,taub,r,p0,d,theta,n,eta,rho,gamma,omegab)

data$p1=data$d+data$p0
data$omegas=1-data$omegab


rm(a,b,taub,tau,r,p0,d,HR,omegab,theta,n,eta,rho,gamma)
# save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/database_Scenarios_H0_True.RData")

#####################################################################################

# simulation seed
set.seed(1987)
# nsim: number of simulations
nsim=100000

#####################################################################################
##################################################################################### 

t0=Sys.time()  

##################################################################################### 
#####################################################################################

Test_Alpha_pluginU <- c()
Test_Alpha_pluginU <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=1,
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p0[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Unpooled')) > z.alpha)/nsim
                                tempMatrix 
                              }

x=t(Test_Alpha_pluginU)
data$Test_Alpha_pluginU <- x 


save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphaplugu.RData")

#####################################################################################

Test_Alpha_pluginP <- c()
Test_Alpha_pluginP <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                              .combine=cbind,
                              .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                                
                                tempMatrix <- sum(replicate(nsim,fCS.TEST_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                             HR=1,
                                                                             rate.param=data$r[i], 
                                                                             p0=data$p0[i], p1=data$p0[i],
                                                                             ass.par=data$theta[i],
                                                                             n0=data$n[i]/2, n1=data$n[i]/2,
                                                                             censoring="Unif",
                                                                             tau=data$tau[i],
                                                                             taub=data$taub[i],
                                                                             rho=data$rho[i], gam=data$gamma[i], 
                                                                             eta=data$eta[i],
                                                                             wb=data$omegab[i], ws=data$omegas[i],
                                                                             var_est='Pooled')) > z.alpha)/nsim
                                
                                tempMatrix
                                
                              }

x=t(Test_Alpha_pluginP)
data$Test_Alpha_pluginP <- x 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphaplugp.RData")

#####################################################################################

Test_Alpha_Boots <- c()
Test_Alpha_Boots <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                            .combine=cbind,
                            .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                              
                              tempMatrix <- sum(replicate(nsim,fCS.TEST_boots_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                                 HR=1,
                                                                                 rate.param=data$r[i], 
                                                                                 p0=data$p0[i], p1=data$p0[i],
                                                                                 ass.par=data$theta[i],
                                                                                 n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                 censoring="Unif",
                                                                                 tau=data$tau[i],
                                                                                 taub=data$taub[i],
                                                                                 rho=data$rho[i], gam=data$gamma[i], 
                                                                                 eta=data$eta[i],
                                                                                 wb=data$omegab[i],
                                                                                 ws=data$omegas[i],
                                                                                 Boot=50)) > z.alpha)/nsim
                              
                              tempMatrix
                              
                            }

x=t(Test_Alpha_Boots)
data$Test_Alpha_Boots <- x 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphaboot.RData")


#####################################################################################

Test_Alpha_Bonf <- c()
Test_Alpha_Bonf <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                           .combine=cbind,
                           .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                             
                             tempMatrix <- sum(replicate(nsim,(sum(fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i],
                                                                                    HR=1,
                                                                                    rate.param=data$r[i], 
                                                                                    p0=data$p0[i], p1=data$p0[i],
                                                                                    ass.par=data$theta[i],
                                                                                    n0=data$n[i]/2, n1=data$n[i]/2,
                                                                                    censoring="Unif",
                                                                                    tau=data$tau[i],
                                                                                    taub=data$taub[i],
                                                                                    rho=data$rho[i], gam=data$gamma[i],
                                                                                    eta=data$eta[i])>z.alphac)>=1)))/nsim
                             
                             tempMatrix
                             
                           }

x=t(Test_Alpha_Bonf)
data$Test_Alpha_Bonf <- x 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphabonf.RData")

#####################################################################################

Test_Alpha_S <- c()
Test_Alpha_S <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=1,
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p0[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[2]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Alpha_S)
data$Test_Alpha_S <- x 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphas.RData")

#####################################################################################

Test_Alpha_B <- c()
Test_Alpha_B <- foreach(i= 1:dim(data)[1], #i= 1:2, # just for testing
                        .combine=cbind,
                        .packages=c("copula","pracma","survival","zoo","muhaz","tidyr","purrr"))%dopar% { 
                          
                          tempMatrix <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], 
                                                                            HR=1,
                                                                            rate.param=data$r[i], 
                                                                            p0=data$p0[i], p1=data$p0[i],
                                                                            ass.par=data$theta[i],
                                                                            n0=data$n[i]/2, n1=data$n[i]/2,
                                                                            censoring="Unif",
                                                                            tau=data$tau[i],
                                                                            taub=data$taub[i],
                                                                            rho=data$rho[i], gam=data$gamma[i],
                                                                            eta=data$eta[i])[1]) > z.alpha)/nsim
                          
                          tempMatrix
                          
                        }

x=t(Test_Alpha_B)
data$Test_Alpha_B <- x 

save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0_alphab.RData")

#####################################################################################

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_H0.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/Marta/Nextcloud/R_code/Code_PAPER/Rev_Simulation/results/RESULTS_PAPER_H0.RData")

##################################################################################### 
#####################################################################################


# stop cluster
stopCluster(cl)


