DIR=getwd()
setwd("E:\\Rcode\\data")

# cases_choosed <- read.csv("sample_choosed.csv")
clinical_all <- read.csv("clinical_merged.csv",row.names = 1)
clinical_all <- clinical_all[!duplicated(clinical_all$submitter_id),]
summary(clinical_all$primary_diagnosis)
# clinical_choosed <- clinical_all[which(clinical_all$submitter_id%in%cases_choosed$Case.ID),]

clinical_choosed <- read.csv("clinical_choosed.csv",row.names = 1)
summary(clinical_choosed$primary_diagnosis)
