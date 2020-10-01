
# TEST - toy examples


#########
# BINARY STATISTIC

data$Test_Alpha_B=0
data$Test_Power_B=0

# for(i in 1:3){ # just for testing
for(i in 1:dim(data)[1]){
  
  data$Test_Alpha_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p0[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim
  
  data$Test_Power_B[i] <- sum(replicate(nsim,fCS.TEST_Bonf_H1(a.shape=data$a[i], b.scale=data$b[i], HR=data$HR[i],
                                                              rate.param=data$r[i], p0=data$p0[i], p1=data$p1[i],
                                                              ass.par=data$theta[i],
                                                              n0=data$n[i]/2, n1=data$n[i]/2,
                                                              censoring="Unif",
                                                              tau=data$tau[i],
                                                              taub=data$taub[i],
                                                              rho=data$rho[i], gam=data$gamma[i],
                                                              eta=data$eta[i])[1]) > z.alpha)/nsim
  print(i)
}

summary(data$Test_Alpha_B)
summary(data$Test_Power_B)

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