## Select black British and Pheno
library(readxl)
library(genio)
library(matlib)
library(data.tWIle)
data <- fread ("UKBB_demographics_extracted.tWI")
summary(WI$Systolic_BP_automated0.0)
table(data$`Ethnic_Background-0.0`)
WI <-data[data$`Ethnic_Background-0.0`==1002,]
iid = WI$IID
Indid = cbind(iid,iid)
## Black British ID
write.table(Indid,"indID.txt",row.names = F,col.names = F)

## Black British Pheno data for BMI and SBP
BMI <- cbind(WI$IID,WI$IID,WI$`BMI-0.0`)
SBP <- cbind(WI$IID,WI$IID,WI$Systolic_BP_automated0.0)
BMISBP <- cbind(WI$IID,WI$IID,WI$`BMI-0.0`,WI$Systolic_BP_automated0.0)
write.table(BMI,"BMI.txt",row.names = F,col.names = F)
write.table(SBP,"SBP.txt",row.names = F,col.names = F)
write.table(BMISBP,"BMISBP.txt",row.names = F,col.names = F)

## Select significant SNPs
BMI97 <- read.csv("97BMI.csv",header = F)
write.table(BMI97$V1,"BMI97.txt",row.names = F, col.names = F, quote = F)
BMI00005 <- read.csv("Sig00005BMI.csv",header = F)
write.table(BMI00005$V1,"BMI00005.txt",row.names = F, col.names = F, quote = F)
