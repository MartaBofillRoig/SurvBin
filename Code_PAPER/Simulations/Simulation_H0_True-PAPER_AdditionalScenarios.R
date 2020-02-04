
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Simulating bivariate binary and survival data 
# Marta Bofill and Guadalupe Gómez
################################################################

# Functions for the binary and survival setting; for the covariance computation; and for simulating the binary and time-to-event data
source('C:/Users/Marta.Bofill/Desktop/Code_PAPER/Functions/binary-functions.R')
source('C:/Users/Marta.Bofill/Desktop/Code_PAPER/Functions/survival-functions.R') 
source('C:/Users/Marta.Bofill/Desktop/Code_PAPER/Functions/cov-functions.R')
source('C:/Users/Marta.Bofill/Desktop/Code_PAPER/Functions/simulation-functions.R') 
source('C:/Users/Marta.Bofill/Desktop/Code_PAPER/Functions/lstats-functions.R') 

# Scenarios 
load("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations/Results/database_Final_Scenarios_H0_True.RData")
data=subset(data,data$theta==5)
# Working directory
setwd("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations") 

##################################################################################### 
# PREAMBLE 
##################################################################################### 
library(copula) 
library(pracma) # for using the function 'integral'
library(survival)  
require(zoo) # 'rollmean' function
require(muhaz)

# Parameters
alpha=0.05;  
z.alpha <- qnorm(1-alpha,0,1)  
z.alphac <- qnorm(1-alpha/2,0,1)  

# nsim: number of simulations

nsim=1000 
# nsim=500

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
set.seed(2023)

t0=Sys.time()
data$Test_Reject_Unpooled=0
data$Test_Reject_Pooled=0

# for(i in 1:1){ # just for testing
for(i in 1:dim(data)[1]){
  data$Test_Reject_Unpooled[i] <- sum(replicate(nsim,
                                       fCS.TEST(a.shape=data$a[i], b.scale=data$b[i], rate.param=data$r[i],
                                                prob0=data$p0[i],
                                                ass.par=3,
                                                ss=data$n[i],
                                                censoring="Unif",
                                                tau=data$tau[i], taub=data$taub[i], rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                wb=data$omegab[i], ws=data$omegas[i], var_est='Unpooled')) > z.alpha,na.rm = T)/nsim
 
  
  data$Test_Reject_Pooled[i] <- sum(replicate(nsim,
                                              fCS.TEST(a.shape=data$a[i], b.scale=data$b[i], rate.param=data$r[i], 
                                                       prob0=data$p0[i], 
                                                       ass.par=3, 
                                                       ss=data$n[i], 
                                                       censoring="Unif",
                                                       tau=data$tau[i], taub=data$taub[i], rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i], 
                                                       wb=data$omegab[i], ws=data$omegas[i], var_est='Pooled')) > z.alpha,na.rm = T)/nsim
  
  
  t1=Sys.time()-t0
  cat(i, "\t", data$Test_Reject_Unpooled[i], "\t", data$Test_Reject_Pooled[i], "\t", t1, "\n", file="LOG_PAPER_Results_add.txt", append=TRUE)
  save.image("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations/Results/RESULTS_PAPER_Add.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_PAPER_Results_add.txt", append=TRUE)
(t1)

rm(i) 

save.image("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations/Results/RESULTS_PAPER_Add.RData")
