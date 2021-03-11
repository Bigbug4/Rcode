# DIR=getwd()
# setwd("E:\\Rcode\\data")
# 
# library(survival)
# library(ggplot2)
# library(ggpubr)
# library(survminer)
# 
# os_result.matrix <- read.table("OS_result.matrix1.txt",blank.lines.skip=F,header = T)
# #os_result.matrix <- read.csv("OS_dat_expr.csv",header = T)
# 
# attach(os_result.matrix)

# ggplot(os_result.matrix,aes(x=Status,y=PMEPA1))+geom_boxplot()
# 
# p <- ggboxplot(os_result.matrix, x="Status", y="PMEPA1", color = "Status",
#                palette = "jco", add = "jitter")
# p+stat_compare_means(method = "t.test") 

PMEPA1_group <- ifelse(PMEPA1 > median(PMEPA1),'high','low')
PMEPA1_group <- as.factor(PMEPA1_group)
table(PMEPA1_group)

my.surv <- Surv(Time,Status)
#my.surv <- Surv(OS,vital_status=="dead")
kmfit1 <- survfit(my.surv~PMEPA1_group,data=os_result.matrix)
# summary(kmfit1)
# plot(kmfit1,col =rainbow(2),main='Overall Survival PMEPA1 ',xlab='Days',ylab='Percent Survival')
# legend("topright", legend=c(levels(PMEPA1_group)), col=rainbow(2), lwd=2)

ggsurvplot(kmfit1,conf.int =F, pval = T,
           ggtheme = theme_bw())

# ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE, fun = "event",
#            ggtheme = theme_bw())

# str(os_result.matrix,  no.list = T, vec.len = 2)
# m <- coxph(my.surv ~ PMEPA1,data = os_result.matrix)
# ggsurvplot(survfit(m, data =  os_result.matrix), palette = "#2E9FDF",
#            ggtheme = theme_minimal())
# 
# beta <- coef(m)
# se <- sqrt(diag(vcov(m)))
# HR <- exp(beta)
# HRse <- HR * se
# p_value <-  summary(m)$sctest[3]
# 
# diff <- survdiff(my.surv ~ PMEPA1,data = os_result.matrix)
# # pvalue <- 1-pchisq(diff$chisq,df=1)
# pvalue <- 1 - pchisq(diff$chisq, length(data.survdiff$n) - 1)

# detach(os_result.matrix)
# save.image("Surv.RData")
# load("Surv.RData")
