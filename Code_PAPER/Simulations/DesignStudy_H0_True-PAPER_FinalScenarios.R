
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Defining scenarios H0 true  
# Marta Bofill and Guadalupe Gómez
################################################################

setwd("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations")

#####################################################################################
#####################################################################################
# DESIGN OF SIMULATION STUDY
# Design
#####################################################################################

# install.packages("rms");install.packages("copula");install.packages("survival");install.packages("plyr")
library(survival)
library(rms) 
library(copula)
library(plyr) 

#####################################################################################
# Defining the database 
# X0 =  'a': shape parameter weibull
# X0.1 = 'b':  scale parameter weibull
# X0.2 = 'tau': total follow-up
# X0.3 = 'r': parameter censoring distribution
# X0.4 = 'p0':   probability of observing the binary endpoint 
# X0.5 = 'copula': copula chosen for the binary-survival joint distribution
# X0.6 = 'theta': association parameter copula
# X0.7 = 'n': sample size per group 

#####################################################################################
# SURVIVAL 
a=c(0.5,1,2) 
b=c(1)
tau=1
taub= c(0.5,1) 
r=c(1,3) #unif 
p0=c(0.2, 0.4)  
# copula=c(1)  # fgmCopula(), 1: frankCopula(), 2: claytonCopula()
theta=c(0.001,2,5)  
n=c(1000) 
# test
eta=c(0,1)
rho=c(0,1)
gamma=c(0,1)
omegab=0.5

# Notes:
# > rho(copula=frankCopula(param=0.001,dim=2))
# [1] 0.0001666667
# > rho(copula=frankCopula(param=2,dim=2))
# [1] 0.3168122
# > rho(copula=frankCopula(param=5,dim=2))
# [1] 0.6434871

#####################################################################################


data <- data.frame(0,0,0,0,0,0,0,0,0,0,0,0,0)
data <- rename(data,c('c(0)'='Scenario',
                      X0 = 'a', 
                      X0.1 = 'b',
                      X0.2 = 'tau',
                      X0.3 = 'r',
                      X0.4 = 'p0',
                      X0.5 = 'taub',
                      X0.6 = 'theta',
                      X0.7 = 'n',
                      X0.8 = 'omegab',
                      X0.9 = 'omegas',
                      X0.10 = 'eta',
                      X0.11 = 'rho',
                      X0.12 = 'gamma'))



# Creating the scenarios
i=1;j=1;k=1;l=1;m=1;u=1;v=1;w=1;
it=1;  

for(i in 1:length(a)){ 
  for(k in 1:length(b)){
    for(l in 1:length(taub)){
      for(m in 1:length(p0)){
        for(u in 1:length(theta)){
          for(v in 1:length(r)){
            for(j1 in 1:length(eta)){
              for(j2 in 1:length(rho)){
                for(j3 in 1:length(gamma)){
                  if(eta[j1]+rho[j2]+gamma[j3]>=1){
                    for(w in 1:length(n)){ 
                      data[it,]<- c(a[i],b[k],tau,r[v],p0[m],taub[l],theta[u],n[w],omegab,1-omegab,eta[j1],rho[j2],gamma[j3])
                      it=it+1;
                    }
                  }
                }
              }
            }  
          }
        }
      }
    }
  } 
}


rm(i,j,k,l,m,u,v,w,it,j1,j2,j3)
rm(a,b,taub,tau,r,p0,omegab,theta,n,eta,rho,gamma)

##################################################################################### 
save.image("C:/Users/Marta.Bofill/Desktop/Code_PAPER/Simulations/Results/database_Final_Scenarios_H0_True.RData")

##################################################################################### 
