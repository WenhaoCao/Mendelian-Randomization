setwd("D:/MR_Analysis/Newdata/UKB")

library(MendelianRandomization)

ReadGRMBin=function(prefix, AllN=F, size=4){
  sum_i=function(i){
    return(sum(1:i))
  }
  BinFileName=paste(prefix,".grm.bin",sep="")
  NFileName=paste(prefix,".grm.N.bin",sep="")
  IDFileName=paste(prefix,".grm.id",sep="")
  id = read.table(IDFileName)
  n=dim(id)[1]
  BinFile=file(BinFileName, "rb");
  grm=readBin(BinFile, n=n*(n+1)/2, what=numeric(0), size=size)
  NFile=file(NFileName, "rb");
  if(AllN==T){
    N=readBin(NFile, n=n*(n+1)/2, what=numeric(0), size=size)
  }
  else N=readBin(NFile, n=1, what=numeric(0), size=size)
  i=sapply(1:n, sum_i)
  return(list(diag=grm[i], off=grm[-i], id=id, N=N))
}

## the univariate MR method
myGRM = ReadGRMBin(prefix="allkeepCPDLD0.05")
BMISBP <- read.table("BMISBP.txt", header = F)
data <- data.frame(ID = BMISBP[,1], BMI = BMISBP[,3], SBP = BMISBP[,4])

newdata <- data.frame(id=myGRM$id[,1], BMI=rep(0,nrow(myGRM$id)),SBP =rep(0,nrow(myGRM$id))) 
A=snpproduct=myGRM$off
for (i in 1:nrow(myGRM$id)) {
  newdata$BMI[i]= data$BMI[data$ID==myGRM$id[i,1]]
  newdata$SBP[i]=data$SBP[data$ID==myGRM$id[i,1]]
}
newdata$BMI <- scale(newdata$BMI)
newdata$SBP <- scale(newdata$SBP)

n=length(newdata$id)
k=1
number=n*(n-1)/2

colno=(k+1)*k
X=newdata[,2:3]

pheno=matrix(NA,nrow=number,3)
l=1

XX= c()
T1=X[,1]%*%t(X[,1])
for (i in 2:n) {
  a = T1[i,1:(i-1)]
  XX=c(XX,a)
}
pheno[,1]=XX

T2=X[,2]%*%t(X[,1])
XY= c()
for (i in 2:n) {
  a = T2[i,1:(i-1)]
  XY=c(XY,a)
}
pheno[,2]=XY

A=snpproduct
mydata = data.frame(XX=pheno[,1], XY=pheno[,2] ,A=snpproduct)
mydata = mydata[!is.na(mydata$XY),]
mydata = mydata[!is.na(mydata$XX),]


theta1=
  cov(mydata$XY,mydata$A)/cov(mydata$XX,mydata$A)

newMRSE <- function(XX,XY,A){
  eta = mean(XY*A)
  delta = mean(XX*A)
  theta = eta/delta
  se = sqrt( (1/delta^2)*(var(A*XY)+var(A*XX)*theta^2-2*theta*cov(A*XY,A*XX)))
  return(c(eta,delta,theta,se))
}

se1 <- newMRSE(XX=mydata$XX, XY=mydata$XY, A=mydata$A)
se1 =sqrt((sum((mydata$A^2))/(t(mydata$A)%*%mydata$XX)^2)*(sum((mydata$XY-theta1*mydata$XX)^2)/number))


## other method
sigSBPkG <- read.table("fastGWASBPLD01.FASTGWA",header = T)
sigBMIkG <- read.table("fastGWABMILD01.FASTGWA",header = T)
sigSNPBMI005 <- sigBMIkG[sigBMIkG$P<0.005,]
sigSNPBMI0005 <- sigBMIkG[sigBMIkG$P<0.0005,]
sigSNPBMI00005 <- sigBMIkG[sigBMIkG$P<0.00005,]
write.table(sigSNPBMI005$SNP,"sigBMI005.txt", row.names = F, col.names = F, quote = F)

write.table(sigSNPBMI0005$SNP,"sigBMI0005.txt", row.names = F, col.names = F, quote = F)

write.table(sigSNPBMI00005$SNP,"sigBMI00005.txt", row.names = F, col.names = F, quote = F)

plot(sigSNPBMI$BETA)

betaYG = sigSBPkG$BETA
betaXG = sigBMIkG$BETA
sebetaYG = sigSBPkG$SE
sebetaXG = sigBMIkG$SE
pYG = sigSBPkG$P
pXG = sigBMIkG$P

Top <- order(pXG)[1:56]
MRTop<-mr_allmethods(mr_input(bx = betaXG[sigBMIkG$P<0.00005],
                              bxse = sebetaXG[sigBMIkG$P<0.00005],
                              by = betaYG[sigBMIkG$P<0.00005],
                              byse = sebetaYG[sigBMIkG$P<0.00005]), method="main", iterations = 5000)



myGRM = ReadGRMBin(prefix="allsig005")
BMISBP <- read.table("BMISBP.txt", header = F)
data <- data.frame(ID = BMISBP[,1], BMI = BMISBP[,3], SBP = BMISBP[,4])

newdata <- data.frame(id=myGRM$id[,1], BMI=rep(0,nrow(myGRM$id)),SBP =rep(0,nrow(myGRM$id))) 
A=snpproduct=myGRM$off
for (i in 1:nrow(myGRM$id)) {
  newdata$BMI[i]= data$BMI[data$ID==myGRM$id[i,1]]
  newdata$SBP[i]=data$SBP[data$ID==myGRM$id[i,1]]
}
newdata$BMI <- scale(newdata$BMI)
newdata$SBP <- scale(newdata$SBP)

n=length(newdata$id)
k=1
number=n*(n-1)/2

colno=(k+1)*k
X=newdata[,2:3]

pheno=matrix(NA,nrow=number,3)
l=1

XX= c()
T1=X[,1]%*%t(X[,1])
for (i in 2:n) {
  a = T1[i,1:(i-1)]
  XX=c(XX,a)
}
pheno[,1]=XX

T2=X[,2]%*%t(X[,1])
XY= c()
for (i in 2:n) {
  a = T2[i,1:(i-1)]
  XY=c(XY,a)
}
pheno[,2]=XY

T3=X[,2]%*%t(X[,2])
YY= c()
for (i in 2:n) {
  a = T3[i,1:(i-1)]
  YY=c(YY,a)
}
pheno[,3]=YY
A=snpproduct
mydata = data.frame(XX=pheno[,1], XY=pheno[,2] ,A=snpproduct)
mydata = mydata[!is.na(mydata$XY),]
mydata = mydata[!is.na(mydata$XX),]
summary(lm(XX~A))
summary(lm(XY~A))
summary(lm(YY~A))

theta2=
  cor(mydata$XY,mydata$A)/cor(mydata$XX,mydata$A)

newMRSE <- function(XX,XY,A){
  eta = mean(XY*A)
  delta = mean(XX*A)
  theta = eta/delta
  se = sqrt( (1/delta^2)*(var(A*XY)+var(A*XX)*theta^2-2*theta*cov(A*XY,A*XX)))
  return(c(eta,delta,theta,se))
}

se1 <- newMRSE(XX=mydata$XX, XY=mydata$XY, A=mydata$A)
se1 =sqrt((sum((mydata$A^2))/(t(mydata$A)%*%mydata$XX)^2)*(sum((mydata$XY-theta2*mydata$XX)^2)/number))


