DIR=getwd()
setwd("E:\\Rcode\\data")

library(rms)
library(survival)
library(survminer)
library(survMisc)
library(data.table)
library(plyr)

os_result<-read.table("OS_result1.txt",blank.lines.skip=F,header = T)
os_result.matrix <- read.table("OS_result.matrix1.txt",blank.lines.skip=F,header = T)
# os_result<- na.omit(os_result)

os_result_m <- os_result[which(os_result$P_value < 0.05),]
x <- as.character(os_result_m$Terms)

sig_gene <- read.csv("sig_gene_v3.csv")
sig_gene <- sig_gene[which(sig_gene$geneName%in%x),]

attach(os_result.matrix)
my.surv <- Surv(Time,Status)

os_result_c <- os_result[which(os_result$P_value < 0.05),]
x1 <- as.character(os_result_c$Terms)
x1
m <- coxph(my.surv ~  LMNB2 + BGN + IRAK1 + MFSD12 + FKBP10 + SOX4 +  SLC12A7 + SLC1A5 + TIMP1 + ENTPD6 + GPX3 + HELZ2 + PMEPA1 + DNMT1,
           data = os_result.matrix)
m <- coxph(my.surv ~  LMNB2 + BGN + IRAK1 + MFSD12 + FKBP10 + SOX4 +  SLC12A7 + SLC1A5 + TIMP1,
           data = os_result.matrix)
m
summary(m)

m1 <- coxph(my.surv ~  LMNB2 + BGN +  MFSD12 + SOX4,  data = os_result.matrix)
m1

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
m2 <- coxph(my.surv ~  risk_score,  data = os_result.matrix)
m2
summary(m2)
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


Sys.setlocale('LC_ALL','C')
library(survivalROC)

mayo=os_result.matrix[,c(2,3,4,5,7,9,59)]
nobs <- NROW(mayo)
cutoff <- 365*5

Mayo5= survivalROC(Stime=mayo$Time,##生存时间
                     status=mayo$Status,## 终止事件    
                     marker = mayo$risk_score, ## marker value    
                     predict.time = cutoff,## 预测时间截点
                     span = 0.05*nobs^(-0.20))##span,NNE法的namda
str(Mayo5)## list结构

plot(Mayo5$FP, Mayo5$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo5$AUC,3)), ##连接
     ylab="TP",
     main="Mayoscore, Method = NNE \n  Year = 5")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

Mayo5.1= survivalROC(Stime=mayo$Time,##生存时间
                   status=mayo$Status,## 终止事件    
                   marker = mayo$risk_score, ## marker value    
                   predict.time = cutoff,## 预测时间截点
                   method="KM")
str(Mayo5.1)## list结构

plot(Mayo5.1$FP, Mayo5.1$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo5.1$AUC,3)), ##连接
     ylab="TP",
     main="Mayoscore, Method = KM \n  Year = 5")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

Mayo3= survivalROC(Stime=mayo$Time,##生存时间
                     status=mayo$Status,## 终止事件    
                     marker = mayo$risk_score, ## marker value    
                     predict.time = 365*3,## 预测时间截点
                     method="KM")##span,NNE法的namda
str(Mayo3)## list结构

plot(Mayo3$FP, Mayo3$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo3$AUC,3)), ##连接
     ylab="TP",
     main="Mayoscore, Method = KM \n  Year = 5")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

df=data.frame(cbind(Mayo5$FP,Mayo5$TP))
ggplot(data = df) + geom_path(aes(X1,X2,colour = "#0F0"),size = 1,linetype=1) +
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1, color="grey", linetype="dashed") + theme_minimal()

df1=data.frame(cbind(Mayo5.1$FP,Mayo5.1$TP))
ggplot(data = df1) + geom_path(aes(X1,X2,colour = "#0F0"),size = 1,linetype=1) +
  labs(title = "ROC curve") + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1, color="grey", linetype="dashed") + theme_minimal()

Mayo5.2= survivalROC(Stime=mayo$Time,##生存时间
                     status=mayo$Status,## 终止事件    
                     marker = mayo$LMNB2, ## marker value    
                     predict.time = cutoff,## 预测时间截点
                     span = 0.25*nobs^(-0.20))
str(Mayo5.2)## list结构

plot(Mayo5.2$FP, Mayo5.2$TP, ## x=FP,y=TP
     type="l",col="red", ##线条设置
     xlim=c(0,1), ylim=c(0,1),   
     xlab=paste( "FP", "\n", "AUC = ",round(Mayo5.2$AUC,3)), ##连接
     ylab="TP",
     main="Mayoscore, Method = KM \n  Year = 5")## \n换行符
abline(0,1,col="gray",lty=2)##线条颜色

detach(os_result.matrix)
save.image("cox1.RData")
load("cox1.RData")
