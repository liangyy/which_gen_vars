
library(CENTIPEDE)

annoFolder="pwms/seed/pwmScan.sites/"

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
pwmName<-cargs[1];
if(length(cargs)>=2)
expName<-cargs[2];
if(length(cargs)>=3)
outDir<-cargs[3];

cat(pwmName," - ",expName,"\n");

fileAnno <- paste(annoFolder,pwmName,'.bed.gz',sep="");
anno <- read.table(fileAnno,sep='\t',as.is=T)
colnames(anno) <- c("chr","start","end","id","PwmScore","strand","nothing")

Mlen <- anno$end[1]-anno$start[1]+1

# Load the DNase data
fileX8b <- paste('/dev/shm/',expName,'.x8b',sep="")
con <- pipe(paste("bedXbFileXtract ",fileX8b,fileAnno," stdout -window=100"))
cutsite <- read.table(con,sep='\t',as.is=T,fill=T)
cutsite <- as.matrix(cutsite[,-1])
r <- rowSums(cutsite)
good.rows <- !is.na(r);
cutsite <- cutsite[good.rows,]
anno <- anno[good.rows,]
cutsite[cutsite > 32] <- 32

## Fit the model with pwm annotation and an intercept 
c.fit <- fitCentipede(Xlist=list(DNase1=cutsite),Y=as.matrix(data.frame(Ict=1,Pwm=anno$PwmScore)), DampLambda = 0.001, DampNegBin = 0.001)

res <- anno[,1:6]
res$PostPr <- c.fit$PostPr
res$LogRatio <- c.fit$LogRatios

lfit <- glm(Post ~ Pwm, data.frame(Post=c.fit$PostPr,Pwm=anno$PwmScore),family=binomial())
aux=summary(lfit)$coefficients
Zscore <- aux[,3]
BetaLogit <- aux[,1]

lambda <- c.fit$LambdaParList$DNase1;

system(paste("mkdir -p ",outDir,"/",expName,"/",sep=""))

fileName <- paste(outDir,"/",expName,"/",pwmName,".Beta.txt",sep="")
cat(pwmName,BetaLogit,Zscore,"\n",file=fileName,sep="\t")

fileName <- paste(outDir,"/",expName,"/",pwmName,".lambda.txt",sep="")
cat(lambda,"\n",file=fileName,sep="\t")

aux <- good.rows
aux[!good.rows] <- NA
aux[good.rows] <- res$LogRatio

fileName <- paste(outDir,"/",expName,"/",pwmName,".PostPr",sep="")
write.table(aux,file=fileName,col.names=F,row.names=F,quote=F,sep="\t")
system(paste("gzip -f ",fileName))

fileName <- paste(outDir,"/",expName,"/",pwmName,".Model.Rd",sep="")
save(c.fit,anno,good.rows,file=fileName)
