
################################################################
# Case study - Statistics Binary and Survival outcomes
# Marta Bofill and Guadalupe G�mez
################################################################

rm(list = ls())

# setwd("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/CaseStudy")
# load("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/CaseStudy/DigitizeIt/Dataset_Survival.RData")

setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/CaseStudy")
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/CaseStudy/DigitizeIt/Dataset_Survival.RData")

set.seed(2020)

n0=dim(pbo_IPD)[1]
n1=dim(trt_IPD)[1]

prob1 = 5.7/100
prob0 = 1.5/100

v0 = runif(n=n0)
v1 = runif(n=n1)

BE0 = ifelse(v0<prob0, 1, 0)
BE1 = ifelse(v1<prob1, 1, 0)

trt_data = data.frame(treat=trt_IPD$arm, time=trt_IPD$time, status=trt_IPD$status, binary = BE1)
tbo_data = data.frame(treat=pbo_IPD$arm, time=pbo_IPD$time, status=pbo_IPD$status, binary = BE0)

data = rbind(trt_data, tbo_data)
head(data)

######################################

sum(subset(data, data$treat==1)$binary)/n1
sum(subset(data, data$treat==0)$binary)/n0

######################################
# PLOTS KAPLAN-MEIER ESTIMATORS
######################################
library(rms)
windows()
fit.rms <- npsurv(survival::Surv(time, status) ~ arm, data=recon_IPD)
survplot(fit  = fit.rms,
         conf = c("none","bands","bars")[2],
         xlab = "Time in years", ylab = "Overall Survival",
         xlim=c(0,5),
         # label.curves = TRUE,                     # label curves directly
         label.curves = list(keys = "lines"),  # legend instead of direct label
         levels.only  = FALSE,                    # show only levels, no label
         abbrev.label = FALSE,                    # if label used, abbreviate
         ## fun = function(x) {1 - x},            # Cumulative probability plot
         loglog   = FALSE,                        # log(-log Survival) plot
         logt     = FALSE,                        # log time
         time.inc = 1,                          # time increment
         dots     = TRUE,                         # dot grid
         n.risk   = TRUE,                          # show number at risk
         cex.n.risk=0.7,
         col=c("coral1","chartreuse4"),lwd=2,lty=1,legend.pos = "topright"
)

library(grDevices)
savePlot(filename="Rplot",
         type=c("png"),
         device=dev.cur())

######################################
# Functions for the binary and survival setting; for the covariance computation; and for simulating the binary and time-to-event data
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/binary-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/survival-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/cov-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/lstats-functions.R')
source('C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Functions/lstats_boots.R')

# source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/binary-functions.R')
# source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/survival-functions.R')
# source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/cov-functions.R')
# source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/lstats-functions.R')
# source('C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Functions/lstats_boots.R')

require(zoo)
require(survival)
require(muhaz)

sigma_sb <- survbinCov(time=data$time, status=data$status, binary=data$binary, treat=data$treat, tau0=0, tau=4, taub=0.5, rho=0, gam=1, eta=1, var_est = "Pooled")
sigma_sb

z_sb = lstats(time=data$time, status=data$status, binary=data$binary, treat=data$treat, tau0=0, tau=4, taub=0.5, rho=0, gam=1, eta=1, wb=0.25, ws=0.75, var_est = "Pooled")
z_sb
# z_sb[1]

z_sb = lstats(time=data$time, status=data$status, binary=data$binary, treat=data$treat, tau0=0, tau=4, taub=0.5, rho=0, gam=1, eta=1, wb=0.25, ws=0.75, var_est = "Unpooled")
z_sb
# z_sb[1]

z_sb = lstats_boots(data$time, data$status, data$binary, data$treat, tau=4, rho=0, gam=1, eta=1, wb=0.25, ws=0.75,
                    Boot = 100)
z_sb


######################################

# install.packages("survRM2")
library(survRM2)

rmst2(time=data$time, status=data$status, arm=data$treat, tau=4)


######################################
# PLOTS HR Over time
######################################
windows()
fit.rms <- npsurv(Surv(time=time,event=status)~treat, data)
hzh <- hazard.ratio.plot(data$treat,Surv(time=data$time,event=data$status),  e=20,  xlim=c(0,30), legendloc='ll',antilog=T, xlab = "Time") #, ylim=c(0,1.5))

se <- hzh$se[1,]
lhr <- hzh$log.hazard.ratio[1,]
t <- hzh$time
hz <- exp(lhr)
LI <- exp(lhr - 2*se)
LS <- exp(lhr + 2*se)
data2 <- data.frame(t,hz,LI,LS)

library(ggplot2)
# ggplot(data2,aes(x=t,y=hz)) + geom_line()
ggplot(data2,aes(x=t,y=hz)) + geom_smooth(colour="#000099") + xlab('Time in years') + ylab("Hazard Ratio")  +theme(axis.text=element_text(size=12), axis.title=element_text(size=12,face="bold")) #+ ylim(0, 2)

# library("survival")
res.cox <- coxph(Surv(time=time,event=status)~treat, data)
res.cox
test.ph <- cox.zph(res.cox)
test.ph


