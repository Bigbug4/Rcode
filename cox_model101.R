DIR=getwd()
setwd("E:\\Rcode\\data")

library(rms)
library(survival)
library(survminer)
library(survMisc)
library(data.table)
library(plyr)

os_result<-read.table("OS_result101.txt",blank.lines.skip=F,header = T)
os_result.matrix <- read.table("OS_result101.matrix.txt",blank.lines.skip=F,header = T)
# os_result<- na.omit(os_result)

os_result_m <- os_result[which(os_result$P_value < 0.05),]
x <- as.character(os_result_m$Terms)

sig_gene <- read.csv("sig_gene_101.csv")
sig_gene <- sig_gene[which(sig_gene$GeneName%in%x),]

attach(os_result.matrix)
my.surv <- Surv(Time,Status)

os_result_c <- os_result[which(os_result$P_value < 0.01),]
x1 <- as.character(os_result_c$Terms)
x1
m <- coxph(my.surv ~  LMNB2 + BGN + IRAK1 + MFSD12 + SLC12A7 + COL5A2 + DNMT1 + SLC1A5 + ENTPD6 + SLC1A5 + TMC6 + NFKB2 + CLN6 + FKBP10 + COL1A1 + MCM7 + HELZ2 + SOX4 + CAD + SULF1 + NCAPD2 + RGS2 + WDR34 + TIMP1,
           data = os_result.matrix)
m

mr <- coxph(my.surv ~  LMNB2 + BGN + MFSD12 + SLC12A7,
            data = os_result.matrix)
mr
os_result_r <- os_result_c[c(1,2,4,6),]
x2 <- as.character(os_result_r$Terms)
detach(os_result.matrix)

kmfit <- survfit(m, data =  os_result.matrix)
summary(kmfit)

ggsurvplot(kmfit, palette= '#2E9FDF',
           ggtheme = theme_minimal())

plot(kmfit,col = "#2E9FDF",main='Overall Survival ',xlab='Days',ylab='Percent Survival')

beta <- coef(m)[c(1,2,4,6)]
se <- sqrt(diag(vcov(m)))[c(1,2,4,6)]
HR <- exp(beta)
HRse <- HR * se
p_value <-  summary(m)$sctest[3]

#risk score=??Expi* ??i
k <- as.numeric(dim(os_result.matrix)[1])
Risk_score <-  data.table()
progress.bar <- create_progress_bar("text")
progress.bar$init(k)

for (n in 1:k){
  temp <- as.data.frame(os_result.matrix[n,c(4,5,7,9)])
  for (i in 1:4){
    temp[i] <- temp[i]*beta[i]
  }
  temp$sum <- sum(temp)
  Risk_score <- rbind(Risk_score, temp)
  progress.bar$step()
}

os_result.matrix$risk_score <- Risk_score$sum
attach(os_result.matrix)
summary(risk_score)
os_result.matrix$risk_group<- ifelse(risk_score > median(risk_score),'high','low')
os_result.matrix$risk_group<- as.factor(os_result.matrix$risk_group)
attach(os_result.matrix)
table(risk_group)
# kmfit
kmfit1 <- survfit(my.surv~risk_group,data=os_result.matrix)
summary(kmfit1)
plot(kmfit1,col =rainbow(2),main='Overall Survival',xlab='Days',ylab='Percent Survival')
legend("topright", legend=c(levels(risk_group)), col=rainbow(2), lwd=2)

ggsurvplot(kmfit1,conf.int =F, pval = T, fun = "event",
           ggtheme = theme_bw())
ggsurvplot(kmfit1,conf.int =F, pval = T,risk.table =T, ncensor.plot = TRUE, 
           ggtheme = theme_bw())

# expr_data_choosed <- read.csv("CPM_S.csv",row.names=1)
# names(expr_data_choosed) <- gsub("\\.","-",names(expr_data_choosed))
# expr_data_choosed <- expr_data_choosed[which(row.names(expr_data_choosed)%in%x2),]
# a <- expr_data_choosed[,which(names(expr_data_choosed)%in%os_result.matrix$Terms)]
# group <- os_result.matrix[which(unique(os_result.matrix$Terms)%in%names(a)),60]
# group <- as.factor(group)
# a <- data.frame(t(a),group)
# a <- a[order(a$group),]
# table(a$group)
# expr_data <- t(a[,1:4])
# pheatmap::pheatmap(expr_data,scale = "row", labels_col = c("high","low"),cluster_cols = FALSE,fontsize_row = 10, fontsize_col = 12)

#predict
# k <- ten(m)
# p <- predict(k)

detach(os_result.matrix)
save.image("cox101.RData")
