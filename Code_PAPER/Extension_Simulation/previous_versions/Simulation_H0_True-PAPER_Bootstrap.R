
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())
# Working directory
setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation")

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
# Scenarios
#####################################################################################

a=c(0.5,1,2)
b=c(1)
tau=1
taub= c(0.5,1)
r=c(1,3) #unif
p0=c(0.2, 0.4)
theta=c(0.001,2,3)
n=c(100)

# test
eta=c(1)
rho=c(0,1)
gamma=c(0,1)
omegab=c(0.5)
# omegab=c(0.25,0.5,0.75)

data = tidyr::expand_grid(a,b,tau,taub,r,p0,theta,n,eta,rho,gamma,omegab)

data$omegas=1-data$omegab

# data = subset(data,data$eta+data$rho+data$gamma>0) #eta=1


rm(a,b,taub,tau,r,p0,omegab,theta,n,eta,rho,gamma)
save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/database_Scenarios_H0_True.RData")

#####################################################################################

# Functions for the binary and survival setting; for the covariance computation; and for simulating the binary and time-to-event data
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/binary-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/survival-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/cov-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/simulation-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/lstats-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/lstats_boots.R')

#####################################################################################
#####################################################################################
# Parameters
alpha=0.05;
z.alpha <- qnorm(1-alpha,0,1)
z.alphac <- qnorm(1-alpha/2,0,1)

# nsim: number of simulations
nsim=1000

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
set.seed(1994)

t0=Sys.time()
# data$Test_Reject_Unpooled=0
# data$Test_Reject_bonferroni=0
data$Test_Reject_boots=0

# for(i in 1:5){ # just for testing
for(i in 1:dim(data)[1]){

  data$Test_Reject_boots[i] <- sum(replicate(nsim,
                                             fCS.TEST_boots(a.shape=data$a[i], b.scale=data$b[i], rate.param=data$r[i],
                                                         prob0=data$p0[i],
                                                         ass.par=data$theta[i],
                                                         ss=data$n[i],
                                                         censoring="Unif",
                                                         tau=data$tau[i], taub=data$taub[i], rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                                         wb=data$omegab[i], ws=data$omegas[i])) > z.alpha,na.rm = T)/nsim




  t1=Sys.time()-t0
  cat(i, "\t", t1, "\n", file="LOG_Results_boots.txt", append=TRUE)
  save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_boots.RData")
}

t1=Sys.time()-t0
cat(t1, "\n", file="LOG_Results_boots.txt", append=TRUE)
(t1)

rm(i)
save.image("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_boots.RData")


#########
library(boot)

a.shape=data$a[i]; b.scale=data$b[i]; rate.param=data$r[i];
prob0=data$p0[i];
ass.par=data$theta[i];
ss=data$n[i];
censoring="Unif";
tau=data$tau[i]; taub=data$taub[i]; rho=data$rho[i]; gam=data$gamma[i]; eta=data$eta[i];
wb=data$omegab[i]; ws=data$omegas[i]; var_est='Unpooled'



db = simsurvbin(a.shape, b.scale, rate.param, prob0, ass.par, ss, censoring)

##
i=1
nsim=1000
test = replicate(nsim,
                 fCS.TEST(a.shape=data$a[i], b.scale=data$b[i], rate.param=data$r[i],
                          prob0=data$p0[i],
                          ass.par=data$theta[i],
                          ss=data$n[i],
                          censoring="Unif",
                          tau=data$tau[i], taub=data$taub[i], rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                          wb=data$omegab[i], ws=data$omegas[i], var_est='Unpooled'))

tboots = replicate(nsim,
                   fCS.TEST_boots(a.shape=data$a[i], b.scale=data$b[i], rate.param=data$r[i],
                                  prob0=data$p0[i],
                                  ass.par=data$theta[i],
                                  ss=data$n[i],
                                  censoring="Unif",
                                  tau=data$tau[i], taub=data$taub[i], rho=data$rho[i], gam=data$gamma[i], eta=data$eta[i],
                                  wb=data$omegab[i], ws=data$omegas[i]))


plot(density(test),ylim=c(0,0.5), lwd = 2)
lines(density(tboots),col=3)
x=rnorm(n=10000)
lines(density(x), col=2)
