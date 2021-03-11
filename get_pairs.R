data<-sum_expr_data[1:5,1:5]
data1<-data[,1]
data1<-as.data.frame(data1)
rownames(data1)<-rownames(data)
data1<-t(data1)
res<-cor(data1,data1,method = "pearson")
res<-as.data.frame(res)
write.csv(res,"res.csv")
res1=read.csv("res.csv",header = TRUE,row.names = 1)
library(data.table)
tmp=data.table()
pair=data.table()
name<-rownames(res1)
for(i in name){
  for(j in name)
    if(res1[i,j]>0.5 && res1[i,j]<1){
      tmp<-rbind(tmp,res1[i,j])
      pair<-rbind(pair,paste(i,paste("-",j)))
    }
    }
tmp
pair
ans=cbind(tmp,pair)
