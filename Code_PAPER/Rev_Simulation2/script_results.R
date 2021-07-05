
################################################################
# SIMULATIONS Statistics Binary and Survival outcomes
# Results
# Marta Bofill and Guadalupe Gómez
################################################################

rm(list = ls())

# Install and load packages
# install.packages("ggplot2")
# install.packages("gridExtra")
library("ggplot2")
library("gridExtra")

# Working directory
setwd("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2")
# setwd("C:/Users/Marta/Nextcloud/Gitkraken/SurvBin/Code_PAPER/Extension_Simulation")

################################################################
# Unified version results
# Under H1
################################################################

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H1_case1.RData")
data$cases = 1
data$tstar = 0

# names(data)
# names(data)[names(data)=="Test_Power_pluginU.V1"] <- 'Test_Power_pluginU'
# names(data)[names(data)=="Test_Power_pluginP.V1"] <- 'Test_Power_pluginP'
# names(data)[names(data)=="Test_Power_Boots.V1"] <- 'Test_Power_Boots'
# names(data)[names(data)==" Test_Power_Bonf.V1"] <- ' Test_Power_Bonf'
# names(data)[names(data)=="Test_Power_S.V1"] <- 'Test_Power_S'
# names(data)[names(data)==" Test_Power_B.V1"] <- ' Test_Power_B'

data_1 = data
dim(data_1)
summary(data_1)

#
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H1_case2.RData")
data$cases = 2
data$tstar = 0
data_2 = data
dim(data_2)
summary(data_2)

#
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H1_case3.RData")
data$cases = 3
data$tstar = 0
data_3 = data
dim(data_3)
summary(data_3)

#
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H1_nPH.RData")
data$cases = 4
data_4 = data
dim(data_4)
summary(data_4)

################################################################
################################################################
#
# Unified dataset for the paper
#
################################################################
################################################################

data_H1 = rbind(data_1,data_2,data_3,data_4)
data_H1$cases = as.factor(data_H1$cases)
summary(data_H1)

# General summary powers
summary(data_H1[,17:22])

# Summary per cases
summary(data_H1[data_H1$cases==1,17:20])
summary(data_H1[data_H1$cases==2,17:22])
summary(data_H1[data_H1$cases==3,17:22])
summary(data_H1[data_H1$cases==4,17:20])


# taub
summary(data_H1[data_H1$cases==1 & data_H1$taub==0.5,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$taub==1,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$taub==0.5,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$taub==1,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$taub==0.5,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$taub==1,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$taub==0.5,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$taub==1,17:20])

# Theta
summary(data_H1[data_H1$cases==1 & data_H1$theta==0.001,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$theta==0.510,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$theta==0.910,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$theta==0.001,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$theta==0.510,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$theta==0.910,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$theta==0.001,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$theta==0.510,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$theta==0.910,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$theta==0.001,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$theta==0.510,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$theta==0.910,17:20])

# Rho and Gamma
summary(data_H1[data_H1$cases==1 & data_H1$gamma==1 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$gamma==0 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$gamma==1 & data_H1$rho==0,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$gamma==0 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$gamma==1 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$gamma==0 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$gamma==1 & data_H1$rho==0,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$gamma==0 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$gamma==1 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$gamma==0 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$gamma==1 & data_H1$rho==0,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$gamma==0 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$gamma==1 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$gamma==0 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$gamma==1 & data_H1$rho==0,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$gamma==0 & data_H1$rho==0,17:20])

# Gamma (late-effects)
summary(data_H1[data_H1$cases==1 & data_H1$gamma==1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$gamma==0,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$gamma==1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$gamma==0,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$gamma==1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$gamma==0,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$gamma==1,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$gamma==0,17:20])

# Rho (censoring)
summary(data_H1[data_H1$cases==1 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$rho==0,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$rho==1,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$rho==0,17:20])

# p0 (more/less effect binary)
summary(data_H1[data_H1$cases==1 & data_H1$p0==0.1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$p0==0.3,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$p0==0.1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$p0==0.3,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$p0==0.1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$p0==0.3,17:20])

summary(data_H1[data_H1$cases==4 & data_H1$p0==0.1,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$p0==0.3,17:20])


# a
summary(data_H1[data_H1$cases==1 & data_H1$a==0.5,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$a==1,17:20])
summary(data_H1[data_H1$cases==1 & data_H1$a==2,17:20])

summary(data_H1[data_H1$cases==2 & data_H1$a==0.5,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$a==1,17:20])
summary(data_H1[data_H1$cases==2 & data_H1$a==2,17:20])

summary(data_H1[data_H1$cases==3 & data_H1$a==0.5,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$a==1,17:20])
summary(data_H1[data_H1$cases==3 & data_H1$a==2,17:20])

# summary(data_H1[data_H1$cases==4 & data_H1$a==0.5,17:20])
summary(data_H1[data_H1$cases==4 & data_H1$a==1,17:20])
# summary(data_H1[data_H1$cases==4 & data_H1$a==2,17:20])

################################################################
# Plot per case
################################################################

colors <- c("0.25" = "#FD6467", "0.75" = "#018F00", "0.5" = "#0001CE", "Individual tests" = "#333C45")

# case1
power_data <- data.frame(power=c(data_H1[data_H1$cases==1,]$Test_Power_pluginU,
                                 data_H1[data_H1$cases==1,]$Test_Power_pluginP,
                                 data_H1[data_H1$cases==1,]$Test_Power_Boots,
                                 data_H1[data_H1$cases==1,]$Test_Power_Bonf,
                                 data_H1[data_H1$cases==1,]$Test_Power_B,
                                 data_H1[data_H1$cases==1,]$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1[data_H1$cases==1,]$Test_Power_pluginU)),
       rep("Pooled",length(data_H1[data_H1$cases==1,]$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1[data_H1$cases==1,]$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1[data_H1$cases==1,]$Test_Power_Bonf)),
       rep("BE",length(data_H1[data_H1$cases==1,]$Test_Power_B)),
       rep("SE",length(data_H1[data_H1$cases==1,]$Test_Power_S))
),
omegab=c(data_H1[data_H1$cases==1,]$omegab,
         data_H1[data_H1$cases==1,]$omegab,
         data_H1[data_H1$cases==1,]$omegab,
         rep(0.5,length(data_H1[data_H1$cases==1,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==1,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==1,]$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_case1 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot() + scale_color_manual(values = colors) + ggtitle("Case 1") +  theme(legend.position = "none")

########

# Create plot with legend
ggp1_legend <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot() + scale_color_manual(values = colors)+theme(legend.position = "bottom")  + labs(color='Weight (wb)')  #+ labs(color='Weight (expression(omega)b)')

# Create user-defined function, which extracts legends from ggplots
extract_legend <- function(my_ggp) {
  step1 <- ggplot_gtable(ggplot_build(my_ggp))
  step2 <- which(sapply(step1$grobs, function(x) x$name) == "guide-box")
  step3 <- step1$grobs[[step2]]
  return(step3)
}

# Apply user-defined function to extract legend
shared_legend <- extract_legend(ggp1_legend)

########

# case2
power_data <- data.frame(power=c(data_H1[data_H1$cases==2,]$Test_Power_pluginU,
                                 data_H1[data_H1$cases==2,]$Test_Power_pluginP,
                                 data_H1[data_H1$cases==2,]$Test_Power_Boots,
                                 data_H1[data_H1$cases==2,]$Test_Power_Bonf,
                                 data_H1[data_H1$cases==2,]$Test_Power_B,
                                 data_H1[data_H1$cases==2,]$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1[data_H1$cases==2,]$Test_Power_pluginU)),
       rep("Pooled",length(data_H1[data_H1$cases==2,]$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1[data_H1$cases==2,]$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1[data_H1$cases==2,]$Test_Power_Bonf)),
       rep("BE",length(data_H1[data_H1$cases==2,]$Test_Power_B)),
       rep("SE",length(data_H1[data_H1$cases==2,]$Test_Power_S))
),
omegab=c(data_H1[data_H1$cases==2,]$omegab,
         data_H1[data_H1$cases==2,]$omegab,
         data_H1[data_H1$cases==2,]$omegab,
         rep(0.5,length(data_H1[data_H1$cases==2,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==2,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==2,]$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_case2 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot() + scale_color_manual(values = colors) + ggtitle("Case 2") +  theme(legend.position = "none")

# case3
power_data <- data.frame(power=c(data_H1[data_H1$cases==3,]$Test_Power_pluginU,
                                 data_H1[data_H1$cases==3,]$Test_Power_pluginP,
                                 data_H1[data_H1$cases==3,]$Test_Power_Boots,
                                 data_H1[data_H1$cases==3,]$Test_Power_Bonf,
                                 data_H1[data_H1$cases==3,]$Test_Power_B,
                                 data_H1[data_H1$cases==3,]$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1[data_H1$cases==3,]$Test_Power_pluginU)),
       rep("Pooled",length(data_H1[data_H1$cases==3,]$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1[data_H1$cases==3,]$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1[data_H1$cases==3,]$Test_Power_Bonf)),
       rep("BE",length(data_H1[data_H1$cases==3,]$Test_Power_B)),
       rep("SE",length(data_H1[data_H1$cases==3,]$Test_Power_S))
),
omegab=c(data_H1[data_H1$cases==3,]$omegab,
         data_H1[data_H1$cases==3,]$omegab,
         data_H1[data_H1$cases==3,]$omegab,
         rep(0.5,length(data_H1[data_H1$cases==3,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==3,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==3,]$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_case3 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot() + scale_color_manual(values = colors) + ggtitle("Case 3") +  theme(legend.position = "none")


# case4
power_data <- data.frame(power=c(data_H1[data_H1$cases==4,]$Test_Power_pluginU,
                                 data_H1[data_H1$cases==4,]$Test_Power_pluginP,
                                 data_H1[data_H1$cases==4,]$Test_Power_Boots,
                                 data_H1[data_H1$cases==4,]$Test_Power_Bonf,
                                 data_H1[data_H1$cases==4,]$Test_Power_B,
                                 data_H1[data_H1$cases==4,]$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1[data_H1$cases==4,]$Test_Power_pluginU)),
       rep("Pooled",length(data_H1[data_H1$cases==4,]$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1[data_H1$cases==4,]$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1[data_H1$cases==4,]$Test_Power_Bonf)),
       rep("BE",length(data_H1[data_H1$cases==4,]$Test_Power_B)),
       rep("SE",length(data_H1[data_H1$cases==4,]$Test_Power_S))
),
omegab=c(data_H1[data_H1$cases==4,]$omegab,
         data_H1[data_H1$cases==4,]$omegab,
         data_H1[data_H1$cases==4,]$omegab,
         rep(0.5,length(data_H1[data_H1$cases==4,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==4,]$omegab)),
         rep("Individual tests",length(data_H1[data_H1$cases==4,]$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_case4 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot() + scale_color_manual(values = colors) + ggtitle("Case 4") +   theme(legend.position = "none")

windows()
# Draw plots with shared legend
grid.arrange(arrangeGrob(plot_case1, plot_case2, plot_case3,plot_case4,  ncol = 2),
             shared_legend, nrow = 2, heights = c(10, 1), top = "Empirical Powers")


windows()
grid.arrange(arrangeGrob(plot_case1, plot_case4,  ncol = 2),
             shared_legend, nrow = 2, heights = c(10, 1), top = "Empirical Powers (Effect on both endpoints)")

grid.arrange(arrangeGrob(plot_case2, plot_case3, ncol = 2),
             shared_legend, nrow = 2, heights = c(10, 1), top = "Empirical Powers (Effect only on one endpoint)")


################################################################
################################################################

################################################################
# Unified version results
# Under H0
################################################################

rm(list = ls())

# Note:
alpha=0.05
nsim=100000
sd=sqrt(alpha*(1-alpha)/nsim)
z.alpha <- qnorm(1-alpha,0,1)
c(alpha-z.alpha*sd,alpha+z.alpha*sd)

#
rm(list = ls())

# Load plugin (pooled unpooled)
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0_alphaplugp.RData")
data_H0_1=data

# load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0_alphaboot.RData")
# load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0_alphabonf.RData")
# load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0_alphas.RData")
# load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0_alphab.RData")

# Load boots, bonf, univ
load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_H0.RData")
data_H0_2=data

data = cbind(data_H0_1,data_H0_2[,17:20])
summary(data)

# General summary powers
summary(data[,17:22])

# data=cbind(data,data_aux[,17:20])

# taub
summary(data[data$taub==0.5,17:20])
summary(data[data$taub==1,17:20])

summary(data[data$taub==1 & data$p0==0.1,17:20])
summary(data[data$taub==1 & data$p0==0.3,17:20])


# Theta
summary(data[data$theta==0.001,17:20])
summary(data[data$theta==0.510,17:20])
summary(data[data$theta==0.910,17:20])

# Rho and Gamma
summary(data[data$gamma==1 & data$rho==1,17:20])
summary(data[data$gamma==0 & data$rho==1,17:20])
summary(data[data$gamma==1 & data$rho==0,17:20])
summary(data[data$gamma==0 & data$rho==0,17:20])

# Gamma (late-effects)
summary(data[data$gamma==1,17:20])
summary(data[data$gamma==0,17:20])

# Rho (censoring)
summary(data[data$rho==1,17:20])
summary(data[data$rho==0,17:20])

# p0 (more/less effect binary)
summary(data[data$p0==0.1,17:22])
summary(data[data$p0==0.3,17:22])


# a
summary(data[data$a==0.5,17:20])
summary(data[data$a==1,17:20])
summary(data[data$a==2,17:20])


################################################################
# ADDITIONAL RESULTS
################################################################

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_add.RData")

data_H1=data

# add1

colors <- c("0.25" = "#FD6467", "0.75" = "#018F00", "0.5" = "#0001CE", "Individual tests" = "#333C45")
power_data <- data.frame(power=c(data_H1$Test_Power_pluginU,
                                 data_H1$Test_Power_pluginP,
                                 data_H1$Test_Power_Boots,
                                 data_H1$Test_Power_Bonf,
                                 data_H1$Test_Power_B,
                                 data_H1$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1$Test_Power_pluginU)),
       rep("Pooled",length(data_H1$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1$Test_Power_Bonf)),
       rep("BE",length(data_H1$Test_Power_B)),
       rep("SE",length(data_H1$Test_Power_S))
),
omegab=c(data_H1$omegab,
         data_H1$omegab,
         data_H1$omegab,
         rep(0.5,length(data_H1$omegab)),
         rep("Individual tests",length(data_H1$omegab)),
         rep("Individual tests",length(data_H1$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_add1 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot()  + scale_color_manual(values = colors) + ggtitle("Small Effect Binary Endpoint") +theme(legend.position = "bottom")  + labs(color='Weight (wb)')
# +   theme(legend.position = "none")

windows()
plot_add1


################################################################

load("C:/Users/mbofi/Dropbox/C5/Scripts/GitKraken/survivalbinary/Code_PAPER/Rev_Simulation2/results/RESULTS_PAPER_add2.RData")

data_H1=data

# add2

colors <- c("0.25" = "#FD6467", "0.75" = "#018F00", "0.5" = "#0001CE", "Individual tests" = "#333C45")
power_data <- data.frame(power=c(data_H1$Test_Power_pluginU,
                                 data_H1$Test_Power_pluginP,
                                 data_H1$Test_Power_Boots,
                                 data_H1$Test_Power_Bonf,
                                 data_H1$Test_Power_B,
                                 data_H1$Test_Power_S
),
Test=c(rep("Unpooled",length(data_H1$Test_Power_pluginU)),
       rep("Pooled",length(data_H1$Test_Power_pluginP)),
       rep("Bootstrap",length(data_H1$Test_Power_Boots)),
       rep("Bonferroni",length(data_H1$Test_Power_Bonf)),
       rep("BE",length(data_H1$Test_Power_B)),
       rep("SE",length(data_H1$Test_Power_S))
),
omegab=c(data_H1$omegab,
         data_H1$omegab,
         data_H1$omegab,
         rep(0.5,length(data_H1$omegab)),
         rep("Individual tests",length(data_H1$omegab)),
         rep("Individual tests",length(data_H1$omegab))
)
)

power_data$Test <- factor(power_data$Test,
                          levels = c('Bootstrap', 'Unpooled', 'Pooled', 'Bonferroni', 'BE','SE'),
                          ordered = TRUE)

plot_add2 <- ggplot(power_data, aes(x=Test, y=power, color=as.factor(omegab))) + geom_boxplot()  + scale_color_manual(values = colors) + ggtitle("Small Effect Survival Endpoint") +theme(legend.position = "bottom")  + labs(color='Weight (wb)')
# +   theme(legend.position = "none")

windows()
plot_add2
