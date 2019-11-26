
DIR=getwd()
setwd("E:\\Rcode\\data")

data<-read.csv("expr_all_log.csv")

library(ROSE)

data$class <- as.factor(data$class)
under <- ovun.sample(class ~ ., data, method = "under", N=50,seed = 3)$data
table(under$class)

library(data.table)
library(plyr)
library(e1071)

#交叉验证
k=5
under$id <- sample(1:k, nrow(under), replace = TRUE)
list <- 1:k
# 每次迭代的预测用数据框，测试用数据框
prediction <- data.table()
testsetCopy <- data.table()
acc <- c()
tpr <-c()
spc <- c()
# svm_pred <- data.table()
# svm_y <- data.table()

# 写一个进度条，用来了解CV的进???
progress.bar <- create_progress_bar("text")
progress.bar$init(k)

# k层的函数
for(i in 1:k){
  # 删除i的行，创建训练集
  # 选i的行，创建验证集
  trainingset <- subset(under, id %in% list[-i])[,-57]
  testset <- subset(under, id %in% c(i))[,-57]
  
  #运行一个svm模型
  mymodel <- svm(trainingset$class ~ ., data = trainingset,probability=T)
  summary(mymodel)
  svm.pred <- predict(mymodel, testset[,-56])
  svm.pred_p <- predict(mymodel, testset[,-56],probability=T)
  svm.pred.p <- as.numeric(attr(svm.pred_p,"probabilities")[,1])
  # n <- ifelse(svm.pred>0,1,0) 
  # n <- ifelse(svm.pred==z,1,0)
  # sum(n)
  tab <- table(svm.pred,testset[,56])
  res <- as.numeric(tab)
  acc[i] <- (res[1] + res[4])/sum(res)
  spc[i] <- res[4]/(res[3] + res[4])
  tpr[i] <- res[1]/(res[1] + res[2])
  # 将迭代出的预测结果添加到预测数据框的末尾
  
  svmtest_y <- as.character(testset[,56])
  
  svmtest_y[which(grepl("yes",svmtest_y))] <- 1
  svmtest_y[which(grepl("no",svmtest_y))] <- 0
  
  temp <- as.data.frame(svm.pred.p)
  prediction <- rbind(prediction, temp)
  testsetCopy <- rbind(testsetCopy, as.data.frame(svmtest_y))
  
  # svm_pred <- cbind(svm_pred, temp)
  # svm_y <- cbind(svm_y,as.data.frame(svmtest_y))
  
  progress.bar$step()
}

# 将预测和实际值放在一???
result <- cbind(prediction, testsetCopy[, 1])
names(result) <- c("Predicted", "Actual")
attach(result)

library(ROCR)
pred <- prediction(Predicted, Actual) 
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))

plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

abline(0,1)

library(ggplot2)
df<- data.frame(x = attributes(perf)$x.values[[1]],y = attributes(perf)$y.values[[1]])
save(df,file = "df_SVM.Rdata")

ggplot(data = df) + geom_line(aes(x,y),colour = "#2E9FDF",size = 1) +   
  labs(title = paste("ROC curve (", "AUC = ",auc,")")) + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1, color="grey", linetype="dashed") + theme_minimal()


# a <- c('a','b','c','d','e')
# svm_pred <- as.data.frame(svm_pred)
# names(svm_pred) <- a
# 
# svm_y <- as.data.frame(svm_y)
# names(svm_y) <- a

# ROC
# library(pROC)
# roc1 <- plot.roc(as.numeric(svm_pred$a),as.numeric(svm_y$a), main="Statistical comparison", col="1")
# par(new=TRUE)
# roc2 <- plot.roc(as.numeric(svm_pred$b),as.numeric(svm_y$b), col="2")
# par(new=TRUE)
# roc3 <- plot.roc(as.numeric(svm_pred$c),as.numeric(svm_y$c), col="3")
# par(new=TRUE)
# roc4 <- plot.roc(as.numeric(svm_pred$d),as.numeric(svm_y$d), col="4")
# par(new=TRUE)
# roc5 <- plot.roc(as.numeric(svm_pred$e),as.numeric(svm_y$e), col="5")
# legend("bottomright", legend=c("1", "2","3","4","5"), col=c("1", "2","3","4","5"), lwd=2)
# 

ACC <- mean(acc)
SPC <- mean(spc)
TPR <- mean(tpr)
# AUC <- sum(auc)/k
save.image("SVM_CV5.RData")
#load("SVM_CV5.RData")
