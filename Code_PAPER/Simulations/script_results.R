
################################################################
# SIMULATIONS RESULTS
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())

library(ggplot2)
library(gridExtra)
library(ggpubr)


load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Simulations/Results/RESULTS_PAPER_Unpooled.RData")
data_complete=data
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Simulations/Results/RESULTS_PAPER_Pooled.RData")
data_pooled=data
data_complete$Test_Reject_Pooled = data_pooled$Test_Reject_Pooled
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Simulations/Results/RESULTS_PAPER_Add.RData")
data2=data
data2$theta= rep(3,dim(data2)[1])
data= rbind(data_complete,data2)
data= subset(data,data$theta<5)

dim(data)
head(data)
summary(data)



data$cases = 0
for(i in 1:dim(data)[1]){
  if(data$a[i]==0.5 & data$taub[i]==0.5) data$cases[i]= 1
  if(data$a[i]==1 & data$taub[i]==0.5) data$cases[i]= 2
  if(data$a[i]==2 & data$taub[i]==0.5) data$cases[i]= 3

  if(data$a[i]==0.5 & data$taub[i]==1) data$cases[i]= 4
  if(data$a[i]==1 & data$taub[i]==1) data$cases[i]= 5
  if(data$a[i]==2 & data$taub[i]==1) data$cases[i]= 6
}
data$cases = as.factor(data$cases)
summary(data$cases)

##

data$cases_p = 0
for(i in 1:dim(data)[1]){
  if(data$p0[i]==0.2 & data$taub[i]==0.5) data$cases_p[i]= 1
  if(data$p0[i]==0.4 & data$taub[i]==0.5) data$cases_p[i]= 2

  if(data$p0[i]==0.2 & data$taub[i]==1) data$cases_p[i]= 3
  if(data$p0[i]==0.4 & data$taub[i]==1) data$cases_p[i]= 4
}
data$cases_p = as.factor(data$cases_p)
summary(data$cases_p)

##

data$cases_tau = 0
for(i in 1:dim(data)[1]){
  if(data$theta[i]==0.001 & data$taub[i]==0.5) data$cases_tau[i]= 1
  if(data$theta[i]==2 & data$taub[i]==0.5) data$cases_tau[i]= 2
  if(data$theta[i]==3 & data$taub[i]==0.5) data$cases_tau[i]= 3

  if(data$theta[i]==0.001 & data$taub[i]==1) data$cases_tau[i]= 4
  if(data$theta[i]==2 & data$taub[i]==1) data$cases_tau[i]= 5
  if(data$theta[i]==3 & data$taub[i]==1) data$cases_tau[i]= 6
}
data$cases_tau = as.factor(data$cases_tau)
summary(data$cases_tau)

##

data$i= 1:dim(data)[1]

data$diff_alpha = data$Test_Reject_Unpooled - data$Test_Reject_Pooled

############

windows(height = 7, width = 7)
fig_b <-  ggplot(data, aes(x=i, y=Test_Reject_Unpooled, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))

fig_b1 <-  ggplot(subset(data,data$eta==0), aes(x=i, y=Test_Reject_Unpooled, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))
fig_b2 <-  ggplot(subset(data,data$eta==1), aes(x=i, y=Test_Reject_Unpooled, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))


ggarrange(fig_b, fig_b1,fig_b2, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

