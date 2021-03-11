DIR=getwd()
setwd("E:\\Rcode\\data")

library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)

expr_data<-read.table("expr_all_stomach.txt")
names(expr_data) <- gsub("\\.","-",names(expr_data))
expr_data_log <- log10(expr_data)
write.table(expr_data_log,"expr_all_stomach_log.txt",sep = "\t")
clin_data<-read.table("clinical_matrix.txt",header = TRUE)

# os_matrix<-read.table("OS_result.matrix.txt",blank.lines.skip=F,header = T)
# time <- as.numeric(os_matrix[,1])
# status <- os_matrix[,2]
dat <- clin_data[clin_data$OS>30,]
write.csv(dat,"OS_dat.csv",row.names = FALSE)
table(dat$vital_status)
attach(dat)

ggplot(dat,       
       aes(x = OS, group = vital_status,colour = vital_status,           
           fill = vital_status
       )) + geom_density(alpha = 0.5)

## ?��?KM????????
my.surv <- Surv(OS,vital_status=='dead')
## The status indicator, normally 0=alive, 1=dead. 
## Other choices are TRUE/FALSE (TRUE = death) or 1/2 (2=death). 
kmfit <- survfit(my.surv~1)
summary(kmfit)
plot(kmfit,main='Overall Survival',xlab='Days',ylab='Percent Survival')
detach(dat)

dat_expr <- cbind(dat,t(expr_data[,dat$submitter_id]))
# dat_expr <- na.omit(dat_expr)
# dat_expr <- dat_expr[!is.na(dat_expr$vital_status),]
dat_expr$vital_status <- as.character(dat_expr$vital_status)
write.csv(dat_expr,"OS_dat_expr.csv",row.names = FALSE)
attach(dat_expr)
ggplot(dat_expr,aes(x=vital_status,y=GPX3))+geom_boxplot()

p <- ggboxplot(dat_expr, x="vital_status", y="GPX3", color = "vital_status",
               palette = "jco", add = "jitter")
p+stat_compare_means(method = "t.test") 

GPX3_group <- ifelse(GPX3 > median(GPX3),'high','low')
GPX3_group <- as.factor(GPX3_group)
table(GPX3_group)

kmfit1 <- survfit(my.surv~GPX3_group,data=dat_expr)
summary(kmfit1)
plot(kmfit1,col =rainbow(2),main='Overall Survival GPX3 ',xlab='Days',ylab='Percent Survival')
legend("topright", legend=c(levels(GPX3_group)), col=rainbow(2), lwd=2)

ggsurvplot(kmfit1,conf.int =F, pval = T,
           ggtheme = theme_bw())
ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE, fun = "event",
           ggtheme = theme_bw())

str(dat_expr,  no.list = T, vec.len = 2)
m <- coxph(my.surv ~ GPX3,data = dat_expr)
ggsurvplot(survfit(m, data =  dat_expr), palette = "#2E9FDF",
           ggtheme = theme_minimal())

beta <- coef(m)
se <- sqrt(diag(vcov(m)))
HR <- exp(beta)
HRse <- HR * se
p_value <-  summary(m)$sctest[3]

diff <- survdiff(my.surv ~ GPX3,data = dat_expr)
# pvalue <- 1-pchisq(diff$chisq,df=1)
pvalue <- 1 - pchisq(diff$chisq, length(data.survdiff$n) - 1)

detach(dat_expr)
save.image("Surv.RData")
load("Surv.RData")
