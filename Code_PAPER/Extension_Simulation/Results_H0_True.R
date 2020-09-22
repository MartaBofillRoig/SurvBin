
# RESULTS SIMULATIONS UNDER H0

library(ggplot2)
library(gridExtra)
library(ggpubr)

load("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_boots.RData")
Test_Reject_boots = data$Test_Reject_boots

load("C:/Users/Marta/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Extension_Simulation/results/RESULTS_PAPER_bonf.RData")
data$Test_Reject_boots = Test_Reject_boots

summary(data)

#

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

#


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

#

windows(height = 7, width = 7)
fig_b <-  ggplot(data, aes(x=i, y=Test_Reject_boots, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))

fig_b1 <-  ggplot(subset(data,data$eta==0), aes(x=i, y=Test_Reject_boots, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))
fig_b2 <-  ggplot(subset(data,data$eta==1), aes(x=i, y=Test_Reject_boots, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Significance level (Unpooled variance)", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))


ggarrange(fig_b, fig_b1,fig_b2, ncol=3, nrow=1, common.legend = TRUE, legend="bottom")

#

windows(height = 7, width = 7)
fig_boots <-  ggplot(data, aes(x=i, y=Test_Reject_boots, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Test_Reject_boots", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))


fig_bonf <-  ggplot(data, aes(x=i, y=Test_Reject_bonferroni, color=as.factor(cases_tau)))  +
  geom_point(size=2) + labs(y = "Test_Reject_bonferroni", x="Scenarios", color=expression(tau[b]))  + coord_cartesian(ylim = c(0, 0.1))   + scale_color_manual(name= expression(paste("(",theta, "," ,tau, ")")),labels = c("(0.001,0.5)","(2,0.5)","(3,0.5)","(0.001,1)","(2,1)","(3,1)"), values=c(1,2,3,4,5,6))



fig_boots_box <-  ggplot(data, aes(x=i, y=Test_Reject_boots, color=as.factor(cases_tau)))  +
  geom_boxplot()
fig_bonf_box <-  ggplot(data, aes(x=i, y=Test_Reject_bonferroni, color=as.factor(cases_tau)))  +
  geom_boxplot()

ggarrange(fig_boots,fig_bonf,fig_boots_box,fig_bonf_box, ncol=2, nrow=2, common.legend = TRUE, legend="bottom")
