DIR=getwd()
setwd("E:\\Rcode\\data")

# library(rms)
library(survival)
library(ggplot2)
library(ggpubr)
library(survminer)
library(data.table)
library(plyr)

dat <- read.csv("OS_dat.csv")
attach(dat)
my.surv <- Surv(OS,vital_status=='dead')
detach(dat)
# set.seed(2)

expr_data<-read.csv("CPM_S.csv",row.names=1)
names(expr_data) <- gsub("\\.","-",names(expr_data))
# CPM_sample <- sample(expr_data)[,1:136]
# write.csv(CPM_sample,"CPM_sample.csv")
# write.table(CPM_sample,"CPM_sample.txt",sep = "\t")
expr_data_log <- read.csv("logCPM_S.csv",row.names=1)
names(expr_data_log) <- gsub("\\.","-",names(expr_data_log))
# logCPM_sample <- sample(expr_data_log)[,1:136]
# write.csv(logCPM_sample,"logCPM_sample.csv")
# write.table(logCPM_sample,"logCPM_sample.txt",sep = "\t" )
dat_expr <- expr_data[,dat$submitter_id]
dat_expr.matrix <- cbind(dat,t(dat_expr))
row.names(dat_expr.matrix) <- c(1:271)
# attach(dat_expr.matrix)

# fit.KM=survfit(my.surv~1)
# fit.KM
# plot(fit.KM)

fun1 <- function(values1){
  group <- ifelse(values1>median(values1),'high','low')
  kmfit <- survfit(my.surv~group)
  # plot(kmfit,col =rainbow(2),main=paste0('Overall Survival',values1),xlab='Days',ylab='Percent Survival')
  # legend("topright", legend=c(levels(group)), col=rainbow(2), lwd=2)
  data.survdiff=survdiff(my.surv~group)
  p.val = 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
}

log_rank_p1 <- apply(dat_expr, 1, fun1)
x1 <- names(log_rank_p1[log_rank_p1<0.05])

fun2 <- function(values1){
  m <- coxph(my.surv ~ values1,data = dat_expr_t)
  # ggsurvplot(survfit(m, data =  dat_expr), palette = "#2E9FDF",
  #            ggtheme = theme_minimal())
  
  beta <- coef(m)
  se <- sqrt(diag(vcov(m)))
  HR <- exp(beta)
  HRse <- HR * se
  p_value <-  summary(m)$sctest[3]
  data.survdiff=survdiff(my.surv~values1)
  p.val = 1-pchisq(data.survdiff$chisq,length(data.survdiff$n) - 1)
  tmp <- rbind(coef = beta, se = se,  HR = HR, HRse = HRse,
                z = beta/se, p_value = p_value, p.val = p.val)
  return(tmp)
}

dat_expr_t <- as.data.frame(t(dat_expr))
attach(dat_expr_t)
cox_results <- apply(dat_expr_t, 2, fun2)
log_rank_p2 <- cox_results[6,]
x2 <- names(log_rank_p2[log_rank_p2<0.05])
row.names(cox_results) <- c('coef', 'se',  'HR', 'HRse','z', 'p_value', 'p.val')
cox_results <- as.data.frame(t(cox_results))

DT::datatable(cox_results ,
              extensions = 'FixedColumns',
              options = list(
                #dom = 't',
                scrollX = TRUE,
                fixedColumns = TRUE
              ))

x <- row.names(dat_expr)
x1 <- paste(x,"+")
x1 <- as.factor(x1)
x1
m <- coxph(my.surv ~   MKI67 +   PMEPA1 +  SOX9 + BGN +  IFI6 + GPRC5A +  SOX4 +  TIMP1   
           + FKBP10 +  CDC25B +  ARFGEF3 + PLOD3 +  MET +  SLC12A7 + TYMP + BOP1   
           +  KPNA2 +   VAV2 +  HELZ2 + RUNX1 + LMNB2 +   LPCAT1 +  SH3KBP1 + THEM6
           + RCC2 +    CAD +  DNMT1 +  SLC1A5 +  TMEM63A + CHD7 +    ENTPD6 +  IRAK1 
           + DKC1 +    PLXNA3 +  GTPBP4 +  NCAPD2 +  MFHAS1 +  PARP14 +  PDCD11 +  MTHFD2
           + MFSD12 +  ATP5PF +  RAP1A +   UBL3 +    PINK1 +  MXI1 + GSN +  CITED2
           + METTL7A + KAT2B + PER1 + KLF4 + SLC25A4 + GPX3 +  FHL1 ,
           data = dat_expr)
m

cox_result <- read.csv("muliti_cox_result.csv",row.names=1)
cox_result_r <-cox_result[which(cox_result$p < 0.05),]
beta <- cox_result_r$coef
x2 <- as.character(row.names(cox_result_r))
x2

cox_result.matrix <- as.data.frame(t(dat_expr))
cox_result.matrix <- cox_result.matrix[which(names(cox_result.matrix)%in%x2)]

#risk score=∑Expi* βi
k <- as.numeric(dim(cox_result.matrix)[1])
Risk_score <-  data.table()
progress.bar <- create_progress_bar("text")
progress.bar$init(k)

for (n in 1:k){
  temp <- as.data.frame(cox_result.matrix[n,])
  for (i in 1:7){
    temp[i] <- temp[i]*beta[i]
  }
  temp$sum <- sum(temp)
  Risk_score <- rbind(Risk_score, temp)
  progress.bar$step()
}

cox_result.matrix$risk_score <- Risk_score$sum
attach(cox_result.matrix)
summary(risk_score)
cox_result.matrix$risk_group<- ifelse(risk_score > median(risk_score),'high','low')
cox_result.matrix$risk_group<- as.factor(cox_result.matrix$risk_group)
attach(cox_result.matrix)
table(risk_group)
# kmfit
kmfit1 <- survfit(my.surv~risk_group,data=cox_result.matrix)
summary(kmfit1)
plot(kmfit1,col =rainbow(2),main='Overall Survival',xlab='Days',ylab='Percent Survival')
legend("topright", legend=c(levels(risk_group)), col=rainbow(2), lwd=2)

ggsurvplot(kmfit1,conf.int =F, pval = T, fun = "event",
           ggtheme = theme_bw())
ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE, 
           ggtheme = theme_bw())

detach(cox_result.matrix)
save.image("cox.RData")
