library(CENTIPEDE)

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
pwmName<-cargs[1];
if(length(cargs)>=2)
expName<-cargs[2];

cat(pwmName," - ",expName,"\n");

## Open the model from the initial fit
myFileName <- paste('initial/',expName,'/',pwmName,'.Model.Rd',sep="")
load(myFileName)
annoTop5K <- anno

# Get the positions of the PWM rescan
fileAnno=paste('pwms/recal/pwmRescan.sites/',pwmName,'.bed.gz',sep="");
if(!file.exists(fileAnno)){
	quit(save="no")
}

## Filter all the locations that overlap a SNP, to fit the Centipede mode on sites not overlapping SNPs
system(paste("intersectBed -a pwms/recal/pwmRescan.sites/",pwmName,".bed.gz -b pwms/recal/pwmSNPs.sites/",pwmName,".bed.gz -wa -c | awk '$7==0' | gzip >/tmp/",pwmName,".bed.gz",sep=""))

fileAnno <- paste("/tmp/",pwmName,".bed.gz",sep="")
anno <- read.table(fileAnno,sep='\t',as.is=T)
colnames(anno) <- c("chr","start","end","id","PwmScore","strand","snp")

Mlen <- anno$end[1]-anno$start[1]+1
fileX8b <- paste('/dev/shm/',expName,'.x8b',sep="")
con <- pipe(paste("bedXbFileXtract ",fileX8b,fileAnno," stdout -window=100"))
cutsite <- read.table(con,sep='\t',as.is=T,fill=T)
cutsite <- as.matrix(cutsite[,-1])

r <- rowSums(cutsite)
good.rows <- !is.na(r);
cutsite <- cutsite[good.rows,]
anno <- anno[good.rows,]

## Filter sites with and extreme number of cutsites (i.e., PCR-artifacts)
Tpcr <- quantile(cutsite,0.9995)
cat("Tpcr: ",Tpcr," #: ",sum(cutsite>Tpcr),"\n");
cutsite[cutsite > Tpcr] <- Tpcr

## Running CENTIPEDE with the model that we learned from Top5K
# source("~/piquelab/centipede/centipede/R/fitCentipedeV2.R")
source("src/fitCentipedeV2.R")

new.fit <- fitCentipede3(Xlist=list(DNase1=cutsite),Y=as.matrix(data.frame(Ict=1,Pwm=anno$PwmScore)), 
												 DampLambda = 0.1, DampNegBin = 0.001,sweeps=200);
rm(cutsite)

res <- anno[,1:6]
res$exp <- expName
res$PostPr <- new.fit$PostPr
res$PostLogOdds <- new.fit$LogRatios
res$DataLogRatio <- new.fit$NegBinLogRatio + new.fit$MultiNomLogRatio # LogRatioData?
res$PriorLogRatio <- new.fit$PriorLogRatio

## How well PWM score matches DNase profile
ct <- cor.test(jitter(res$PwmScore),jitter(res$DataLogRatio),method='spearman')
cat("#Spearman_p.val=",ct$p.value,"_rho=",ct$estimate,"\n")

lfit <- glm(Post ~ Pwm, data.frame(Post=plogis(res$DataLogRatio),Pwm=anno$PwmScore) ,family=binomial())
aux=summary(lfit)$coefficients
Zscore <- aux[,3]
BetaLogit2 <- aux[,1]
cat("#New fit Z-score=",Zscore,"\n")

model <- list(LambdaParList=new.fit$LambdaParList, 
              BetaLogit=new.fit$BetaLogit, 
              NegBinParList=new.fit$NegBinParList)

lambda <- model$LambdaParList$DNase1;

## Quit if the PWM doesn't predict the footprints well
if(ct$p.value>6E-7){
	quit(save="no")
}

###############
## save res  ##
###############
resFolder <- paste("full/nonVar/",expName,"/",sep="")
system(paste("mkdir -p ",resFolder,sep=""))
fileName <- paste(resFolder,pwmName,".PostPr.gz",sep="")
write.table(res[res$PostPr>0.90,],file=gzfile(fileName),col.names=F,row.names=F,quote=F,sep="\t")

resFolder <- paste("full/model/",expName,"/",sep="")
system(paste("mkdir -p ",resFolder,sep=""))
fileName <- paste(resFolder,pwmName,".Rdata",sep="")
save(model,file=fileName)


##############################################################
##  Now analyze SNP positions 
##############################################################
fileAnnoSnp <- paste("pwms/recal/pwmSNPs.sites/",pwmName,".bed.gz",sep="")
annoSnp <- read.table(fileAnnoSnp,sep='\t',as.is=T)
colnames(annoSnp) <- c("chr","start","end","id","PwmScoreRef","strandRef","offRef","PwmScoreAlt","strandAlt","offAlt","SnpPos","RefAllele","VarAllele","SeqContext")
con <- pipe(paste("bedXbFileXtract ",fileX8b,fileAnnoSnp," stdout -window=100"))
cutsite2 <- read.table(con,sep='\t',as.is=T,fill=T)
cutsite2 <- as.matrix(cutsite2[,-1])
r2 <- rowSums(cutsite2)
good.rows2 <- !is.na(r2);
cutsite2 <- cutsite2[good.rows2,]
annoSnp <- annoSnp[good.rows2,]
cutsite2[cutsite2 > Tpcr] <- Tpcr

var.fit <- fitCentipede3(Xlist=list(DNase1=cutsite2),Y=as.matrix(data.frame(Ict=1,Pwm=annoSnp$PwmScoreRef)), 
	DampLambda = 0.1, DampNegBin = 0.001,
	LambdaParList=new.fit$LambdaParList, 
	BetaLogit=new.fit$BetaLogit, 
	NegBinParList=new.fit$NegBinParList,
	sweeps=0)
varLogOdds <- var.fit$NegBinLogRatio + var.fit$MultiNomLogRatio


###############################

annoSnp$id=pwmName;
BetaLogit <- new.fit$BetaLogit;

Y <- as.matrix(data.frame(Ict=1,Pwm=annoSnp$PwmScoreRef))
annoSnp$PriorLogRatioRef <- Y %*% BetaLogit

annoSnp$PwmScoreAlt <- as.numeric(annoSnp$PwmScoreAlt)
Y <- as.matrix(data.frame(Ict=1,Pwm=annoSnp$PwmScoreAlt))
annoSnp$PriorLogRatioAlt <- Y %*% BetaLogit

annoSnp$LogRatiosData <- varLogOdds;
annoSnp$PostPrRef <- plogis(annoSnp$LogRatiosData + annoSnp$PriorLogRatioRef)
annoSnp$PostPrAlt <- plogis(annoSnp$LogRatiosData + annoSnp$PriorLogRatioAlt)

aux <-  annoSnp$PwmScoreRef - annoSnp$PwmScoreAlt
summary(aux)

#############################################################
## Save results 

maxPostPr <- pmax(annoSnp$PostPrRef,annoSnp$PostPrAlt)

resFolder <- paste("full/SNPs/",expName,"/",sep="")
system(paste("mkdir -p ",resFolder,sep=""))
fileName <- paste(resFolder,pwmName,".PostPr.gz",sep="")
write.table(annoSnp[maxPostPr>0.99,],file=gzfile(fileName),col.names=F,row.names=F,quote=F,sep="\t")
