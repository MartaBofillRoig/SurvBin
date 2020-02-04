
# install.packages("dplyr")
# install.packages("magrittr")
# install.packages("survminer")
# install.packages("devtools")

library(dplyr)
library(magrittr)
library(survival)
library(survminer)
library(devtools)
library(reconstructKM)
install_github("ryanrsun/reconstructKM")

# read in the CSVs provided by digitizeIt and save them as rda
setwd("C:/Users/Marta.Bofill/Desktop/Code_PAPER/CaseStudy/DigitizeIt") 

# GROUP 1 (Lpi plus gp100)
km_trt_clicks <- read.csv('Data_group1.csv',header = T,dec = ",",sep = ";") %>%
  as.data.frame(.) %>% 
  set_colnames(c('time', 'survival')) 
km_trt_clicks[1, ] <- c(0, 1)
km_trt_clicks$time=km_trt_clicks$time/12
save(km_trt_clicks, file="Group_lpi.rda")

# GROUP 0 (gp100)
km_pbo_clicks <- read.csv('Data_group0.csv',header = T,dec = ",",sep = ";") %>%
  as.data.frame(.) %>%
  set_colnames(c('time', 'survival')) 
km_pbo_clicks[1, ] <- c(0, 1)

km_pbo_clicks$time=km_pbo_clicks$time/12
save(km_pbo_clicks, file="Group_gp100.rda")

# time=c(0,4,8,12,16,20,24,28,32,36,40,44,48,52,56)
# NAR=c(137,106,79,56,38,30,24,18,13,13,8,5,2,1,0)
# NAR=c(136,93,58,32,23,17,16,7,5,5,3,1,0,0,0)
# NAR=c(403,297,223,163,115,81,54,42,33,24,17,7,6,4,0)

# monthsyears = c(0,4,8,12,16,20,24,28,32,36,40,44,48,52,56) 
monthsyears = c(0,4,8,12,16,20,24,28,32,36,40,44,48,52,56)/12

trt_NAR <- data.frame(time=monthsyears, 
                      NAR=c(403,297,223,163,115,81,54,42,33,24,17,7,6,4,0))
pbo_NAR <- data.frame(time=monthsyears, 
                      NAR=c(136,93,58,32,23,17,16,7,5,5,3,1,0,0,0)) 
save(trt_NAR, file="trt_NAR.rda")
save(pbo_NAR, file="pbo_NAR.rda") 

# augment  
trt_aug <- format_raw_tabs(raw_NAR=trt_NAR,
                           raw_surv=km_trt_clicks) 
pbo_aug <- format_raw_tabs(raw_NAR=pbo_NAR,
                           raw_surv=km_pbo_clicks) 
# reconstruct KM only
trt_recon <- KM_reconstruct(aug_NAR=trt_aug$aug_NAR, aug_surv=trt_aug$aug_surv)
pbo_recon <- KM_reconstruct(aug_NAR=pbo_aug$aug_NAR, aug_surv=pbo_aug$aug_surv) 

# put the treatment and control arms into one dataset
trt_IPD <- data.frame(arm=1, time=trt_recon$IPD_time, status=trt_recon$IPD_event)
pbo_IPD <- data.frame(arm=0, time=pbo_recon$IPD_time, status=pbo_recon$IPD_event)
recon_IPD <- rbind(trt_IPD, pbo_IPD)

save.image("C:/Users/Marta.Bofill/Desktop/Code_PAPER/CaseStudy/DigitizeIt/Dataset_Survival.RData")

# plot
recon_KM_fit <- survival::survfit(survival::Surv(time, status) ~ arm, data=recon_IPD)
recon_KM <- survminer::ggsurvplot(recon_KM_fit, data = recon_IPD, risk.table = TRUE, 
                                  palette=c('black', 'red'),
                                  legend=c(0.86,0.9), legend.title='',legend.labs=c('Gp100-alone', 'Ipilimumab-plus-gp100'),
                                  title='Reconstructed',
                                  ylab='Survival Probability (%)', xlab='Time (years)',
                                  tables.y.text=TRUE,
                                  tables.y.text.col=FALSE, risk.table.title='Number at Risk', break.time.by=1,
                                  censor=TRUE, font.x=22, font.y=18, font.tickslab=16, font.legend=22, 
                                  font.subtitle=20, font.caption=20, risk.table.fontsize=7,
                                  tables.theme = survminer::theme_survminer(font.main = 22, font.y=22,
                                                                            font.x=22, font.tickslab=16))
windows(width = 7, height = 7)
recon_KM        

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

