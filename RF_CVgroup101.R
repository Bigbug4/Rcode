
DIR=getwd()
setwd("E:\\Rcode\\data")

data_z <- read.csv("expr_gene_all101.csv",row.names = 1)
data_zscore <- scale(data_z[,-102])
data_zscore <- as.data.frame(data_zscore)
data_zscore$class <- data_z$class
write.csv(data_zscore,"expr_all_zscore101.csv")

data<-read.csv("expr_all_zscore101.csv",row.names = 1)
# data<-read.csv("expr_all_log.csv")

library(ROSE)

data$class <- as.factor(data$class)
under <- ovun.sample(class ~ ., data, method = "under", N=50,seed = 3)$data
table(under$class)

library(data.table)
library(plyr)

library(randomForest)

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
# rf_pred <- data.table()
# rf_y <- data.table()

# 写一个进度条，用来了解CV的进???
progress.bar <- create_progress_bar("text")
progress.bar$init(k)

# k层的函数
for(i in 1:k){
  # 删除i的行，创建训练集
  # 选i的行，创建验证集
  trainingset <- subset(under, id %in% list[-i])[,-103]
  testset <- subset(under, id %in% c(i))[,-103]
  
  #运行一个随机森林模???
  mymodel <- randomForest(trainingset$class ~ ., data = trainingset,importance = TRUE,norm.vote=TRUE)
  #去掉回应列class
  rf.pred <- predict(mymodel, testset[,-102])
  rf.pred.p = predict(mymodel,testset[,-102],type = 'prob')
  tab <- table(rf.pred,testset[,102])
  res <- as.numeric(tab)
  acc[i] <- (res[1] + res[4])/sum(res)
  spc[i] <- res[4]/(res[3] + res[4])
  tpr[i] <- res[1]/(res[1] + res[2])
  # 将迭代出的预测结果添加到预测数据框的末尾
  rftest_rf.pred <- rf.pred.p[,1]
  rftest_y <- as.character(testset[,102])
  
  rftest_y[which(grepl("yes",rftest_y))] <- 1
  rftest_y[which(grepl("no",rftest_y))] <- 0
  
  temp <- as.data.frame(rf.pred.p[,1])
  prediction <- rbind(prediction, temp)
  testsetCopy <- rbind(testsetCopy, as.data.frame(rftest_y))
  
  # rf_pred <- cbind(rf_pred, temp)
  # rf_y <- cbind(rf_y,as.data.frame(rftest_y))
  
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
save(df,file = "df_RF.Rdata")

ggplot(data = df) + geom_line(aes(x,y),colour = "#2E9FDF",size = 1) +   
  labs(title = paste("ROC curve (", "AUC = ",auc,")")) + 
  xlab("1 - Specificity") + 
  ylab("Sensitivity") + 
  theme(plot.title = element_text(size = 15)) +
  geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), size = 1, color="grey", linetype="dashed") + theme_minimal()

# a <- c('a','b','c','d','e')
# rf_pred <- as.data.frame(rf_pred)
# names(rf_pred) <- a
# rf_y <- as.data.frame(rf_y)
# names(rf_y) <- a

#ROC
# library(pROC)
# roc1 <- plot.roc(as.numeric(rf_pred$a),as.numeric(rf_y$a), main="Statistical comparison", col="1")
# par(new=TRUE)
# roc2 <- plot.roc(as.numeric(rf_pred$b),as.numeric(rf_y$b), col="2")
# par(new=TRUE)
# roc3 <- plot.roc(as.numeric(rf_pred$c),as.numeric(rf_y$c), col="3")
# par(new=TRUE)
# roc4 <- plot.roc(as.numeric(rf_pred$d),as.numeric(rf_y$d), col="4")
# roc5 <- lines.roc(as.numeric(rf_pred$e),as.numeric(rf_y$e), col="5")
# legend("bottomright", legend=c("1", "2","3","4","5"), col=c("1", "2","3","4","5"), lwd=2)


ACC <- mean(acc)
SPC <- mean(spc)
TPR <- mean(tpr)
# AUC <- mean(auc)
# save.image("RF_CV5.RData")
#load("RF_CV5.RData")
