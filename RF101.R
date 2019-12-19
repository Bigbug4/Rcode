
DIR=getwd()
setwd("E:\\Rcode\\data")

#randomForest
trainset <- read.csv("expr_gene_choosed101.csv",row.names = 1)
testset <- read.csv("expr_gene_all101.csv",row.names = 1)

trainset_zscore <- scale(trainset[,-102])
trainset_zscore <- as.data.frame(trainset_zscore)
trainset_zscore$class <- trainset$class

testset_zscore <- scale(testset[,-102])
testset_zscore <- as.data.frame(testset_zscore)
testset_zscore$class <- testset$class

trainset <- trainset_zscore 
testset <- testset_zscore

library(randomForest)

set.seed(1245)
rf.model <- randomForest(class ~ .,data = trainset,importance = TRUE,norm.vote=TRUE)
rf.model
plot(rf.model,type="l",main="l")

x <- subset(testset, select = -class)
y <- testset$class
rf.pred = predict(rf.model,x)
rf.pred.p = predict(rf.model,x,type = 'prob')
tab <- table(rf.pred,y)
y <- as.character(y)
y[which(grepl("yes",y))] <- 1
y[which(grepl("no",y))] <- 0
y <- as.numeric(y)
df2 <- data.frame(y,rf.pred.p[,2])
names(df2) <- c("y","pred") 
write.csv(df2,file ="rftest101.csv",row.names = FALSE)

res <- as.numeric(tab)
acc <- (res[1] + res[4])/sum(res)
TPR <- res[4]/(res[3] + res[4])
SPC <- res[1]/(res[1] + res[2])

# value <- importance(rf.model,2)
# varImpPlot(rf.model)
# 
# margins.rf = margin(rf.model,trainset)
# plot(margins.rf)
# hist(margins.rf,main = "Margines of Random Forest for stomach dataset")
# boxplot(margins.rf ~ trainset$class,main = "Margines of Random Forest for stomach dataset by class")

#ROC
rftest <- read.csv("rftest101.csv")

rftest_rf.pred <- as.numeric(rftest$pred)
rftest_y <- as.numeric(rftest$y)

library(ROCR)
pred <- prediction(rftest_rf.pred, rftest_y) 
perf <- performance(pred,"tpr","fpr")
auc <- performance(pred,'auc')
auc = unlist(slot(auc,"y.values"))

plot(perf,
     xlim=c(0,1), ylim=c(0,1),col='red', 
     main=paste("ROC curve (", "AUC = ",auc,")"),
     lwd = 2, cex.main=1.3, cex.lab=1.2, cex.axis=1.2, font=1.2)

abline(0,1)

save.image("RF101.Rdata")
