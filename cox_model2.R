DIR=getwd()
setwd("E:\\Rcode\\data")

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
library(data.table)
library(plyr)

clin_data_all <- read.csv("clinical_merged.csv",header = TRUE)
clin_data_all$days_to_death <- as.character(clin_data_all$days_to_death)
clin_data_all$days_to_death[which(clin_data_all$days_to_death=="--")] <- 0 
clin_data_all$days_to_last_follow_up <- as.character(clin_data_all$days_to_last_follow_up)
clin_data_all$days_to_last_follow_up[which(clin_data_all$days_to_last_follow_up=="--")] <- 0 
clin_data_all$OS <- as.numeric(clin_data_all$days_to_death) + as.numeric(clin_data_all$days_to_last_follow_up)

dat <- clin_data_all[clin_data_all$OS>30,]
summary(clin_data_all)
my.surv <- Surv(dat$OS,dat$vital_status=='dead')

kmfit1 <- survfit(my.surv~dat$gender,data=dat)
summary(kmfit1)
plot(kmfit1,col =rainbow(2),main='Overall Survival Sex ',xlab='Days',ylab='Percent Survival')
legend("topright", legend=c(levels(dat$gender)), col=rainbow(2), lwd=2)

ggsurvplot(kmfit1,conf.int =F, pval = T, fun = "event",
           ggtheme = theme_bw())
ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE,
           ggtheme = theme_bw())

dat_s <- dat[which(dat$tumor_stage!="not reported"),]
dat_s$tumor_stage <- as.character(dat_s$tumor_stage)
dat_s$tumor_stage <- as.factor(dat_s$tumor_stage)
kmfit2 <- survfit(Surv(dat_s$OS,dat_s$vital_status=='dead')~dat_s$tumor_stage,data=dat_s)
summary(kmfit2)
plot(kmfit2,col =rainbow(12),main='Overall Survival Sex ',xlab='Days',ylab='Percent Survival')
legend("topright", legend=c(levels(dat_s$tumor_stage)), col=rainbow(12), lwd=2)

ggsurvplot(kmfit2,conf.int =F, pval = T,fun = "event",
           ggtheme = theme_bw())
ggsurvplot(kmfit2,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE,
           ggtheme = theme_bw())
