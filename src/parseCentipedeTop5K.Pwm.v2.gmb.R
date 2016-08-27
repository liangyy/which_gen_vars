
baseFolder='initial/combo/'

library(centipede)
library(MASS)

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
pwmName<-cargs[1];
if(length(cargs)>=2)
seedFolder<-cargs[2]; # pwms/seed/
if(length(cargs)>=3)
outFolder<-cargs[3];

system(paste("mkdir -p",outFolder));

aux <- read.table("factorNames.txt",sep="\t",as.is=T)
TfNames <- aux$V2
names(TfNames) <- aux$V1

aux <- read.table("dnase/numReads.txt",sep="\t",as.is=T,quote="")
Nreads <- aux$V2
names(Nreads) <- sub(".x8b.bz2.Qsub.e","",aux$V1)

ii <- pwmName
aux <- paste(baseFolder,"/",ii,'.beta.txt.gz',sep="")
aux <- (read.table(aux,as.is=T,sep="\t"))
sNames <- aux$V1

beta <- aux[,2:5]
rownames(beta) <- sNames

jj <- which.min(beta$V2)
tmax <- sNames[jj]
tmax

Tcut <- (qlogis(0.01)-beta[jj,1])/beta[jj,2]


## Load Sequences
seqFile <- paste(seedFolder,'/pwmScan.seq/',pwmName,'.seq.gz',sep="")
seqs <- scan(gzfile(seqFile),what=character(0), sep="")
seqs <- strsplit(seqs,character(0))
seqs <- t(sapply(seqs,cbind))
Sseqs<-dim(seqs)[2]

## Load model max
myFileName <- paste('initial/',tmax,'/',pwmName,'.Model.Rd',sep="")
load(myFileName)
LogRatio <- c.fit$NegBinLogRatio+c.fit$MultiNomLogRatio

PostPr <- plogis(LogRatio);
PostPr[PostPr<0.90] <- 0.0;
pwmmat2 <- ComputePWM(seqs,PostPr)
pwmmat2 <- apply(pwmmat2+1E-6,2,function(col){col/sum(col)})

bkgModel <- table(seqs[PostPr>0.90,])+1E-6
bkgModel <- bkgModel/sum(bkgModel)
ic <- apply(pwmmat2,2,function(col){sum(col*log2(col/bkgModel))})
pwm.ic.range <- range(which(ic>0.1))
pwm.ic.min <- min(pwm.ic.range)
pwm.ic.max <- max(pwm.ic.range)

new.pwm <- rbind(bkgModel,t(pwmmat2[,pwm.ic.min:pwm.ic.max]),bkgModel)
row.names(new.pwm) <- 1:nrow(new.pwm)

## Recalculate PWMscore with new PWM. 
matseq <- cbind(seqs=="A",seqs=="C",seqs=="G",seqs=="T")+0.0;
W2 <- ncol(matseq)

pwmvec2 <- as.numeric(t(pwmmat2))
pwmvec2 <- log2(pwmvec2+1E-6)+2
PwmScore2 <- matseq %*% pwmvec2

bkgModel.vec <- log2(as.numeric(sapply(1:4,function(ii){rep(bkgModel[ii],W2/4)})))+2

PwmScore2b <- matseq %*% (pwmvec2-bkgModel.vec)

rho.0 <- cor.test(anno$PwmScore,LogRatio,method="spearman")
str(rho.0)

rho.1 <- cor.test(PwmScore2,LogRatio,method="spearman")
str(rho.1)

rho.1b <- cor.test(PwmScore2b,LogRatio,method="spearman")
str(rho.1b)

cat("#rho:",pwmName,tmax,rho.0$estimate,rho.1$estimate,rho.1b$estimate,
    rho.0$p.value,rho.1$p.value,rho.1b$p.value,"\n",sep="\t")

## Write new PWM model, and logistic 
## New logistic seq method
fLogit <- function(BetaLogit,Y,Ez){
    -sum(Ez*(Y %*% BetaLogit) - log(1+exp(Y %*% BetaLogit)))
  }
gLogit <- function(BetaLogit,Y,Ez){
    myPi <- plogis(Y %*% BetaLogit)
    -t(Ez-myPi) %*% Y
  }

simpleLogit <-
function(X,Ez){
   BetaFit <- optim(c(0,0),fLogit,gLogit,Y=as.matrix(data.frame(IntCept=1,X=X)),Ez=Ez,method="BFGS",control=list(maxit=500),hessian=TRUE);
   logitSE <- sqrt(diag(ginv(as.matrix(BetaFit$hessian))))
   Zlogit <- BetaFit$par/logitSE
   c(BetaFit$par,Zlogit)
 }

newLogit <- function(X,Ez){
  BetaFit0 <- optim(0,fLogit,gLogit,Y=matrix(1,length(Ez),1),Ez=Ez,method="BFGS",control=list(maxit=500),hessian=TRUE);
  BetaFit <- optim(c(0,1),fLogit,gLogit,Y=as.matrix(data.frame(IntCept=1,X=X)),Ez=Ez,method="L-BFGS-B",lower=c(-100,0.999),upper=c(0,1.001),control=list(maxit=500));
  D <- 2*(BetaFit$value-BetaFit0$value)
  c(BetaFit$par,D) ## Not nested model
}


###############################

beta.new2b <- newLogit(PwmScore2b,plogis(LogRatio))
beta.new2 <- newLogit(PwmScore2,plogis(LogRatio))
beta.new0 <- newLogit(anno$PwmScore,plogis(LogRatio))

cat("#Dbeta:",pwmName,tmax,beta.new0[1],beta.new2[1],beta.new2b[1],
    beta.new0[3],beta.new2[3],beta.new2b[3],"\n",sep="\t")

tpwm <- function(p,beta){
  qlogis(p)-beta[1]/beta[2]
}

p <- c(0.01,0.05,0.1,0.25,0.5,0.75,0.9,0.95,0.99)
Tcut2 <- sapply(p,function(p){tpwm(p,beta.new2)})
Tcut2b <- sapply(p,function(p){tpwm(p,beta.new2b)})

cat("#p.cut",pwmName,tmax,p,"\n",sep="\t")
cat("#p.Tcut2",pwmName,tmax,round(Tcut2,digits=2),"\n",sep="\t")
cat("#p.Tcut2b",pwmName,tmax,round(Tcut2b,digits=2),"\n",sep="\t")

## Output the updated PWM motif model
resFile <- paste(outFolder,pwmName,'.pwm',sep="");
con <- file(resFile,"w");
cat("#",pwmName,"\t",TfNames[pwmName],"\n",sep="",file=con)
cat("#p.cut",pwmName,tmax,p,"\n",sep="\t",file=con)
cat("#p.Tcut2",pwmName,tmax,round(Tcut2,digits=2),"\n",sep="\t",file=con)
cat("#p.Tcut2b",pwmName,tmax,round(Tcut2b,digits=2),"\n",sep="\t",file=con)
write.table(round(new.pwm,digits=8),quote=F,sep="\t",row.names=F,file=con)
close(con)

cat("#COMPLETED.1: ",pwmName,"\n");

#####################
##  PLOT NEW LOGO  ##
#####################
resFile <- paste(outFolder,"/",pwmName,".logo.png",sep="")
png(resFile);
pwm.v0 <- read.table(paste(seedFolder,"/pwmFiles/",pwmName,".pwm",sep=""),skip=1,as.is=T,sep="\t")
par(mfrow=c(2,1))

pwmLogo(t(new.pwm))
title(main=paste(pwmName," New --",TfNames[pwmName]))
ic2 <- apply(new.pwm,1,function(col){sum(col*log2(col/bkgModel))})
#lines(1:length(ic2),ic2,lwd=2)
pwmLogo(t(pwm.v0))
title(main=paste(pwmName," Old --",TfNames[pwmName]))

dev.off()

########################
## PLOT TOP FOOTPRINT ##
########################
resFile <- paste(outFolder,"/",pwmName,".lambda.png",sep="")
png(resFile);
Mlen <- anno$end[1]-anno$start[1]
par(mfrow=c(1,1))
plotProfile(c.fit$LambdaParList$DNase1,Mlen=Mlen,legTitle=paste(tmax,"\n",pwmName,"--",TfNames[pwmName]))
dev.off()

###########################
##  PWM model prediction ##
###########################
resFile <- paste(outFolder,"/",pwmName,".logit.png",sep="")
png(resFile);
x <- seq(0,50,0.1)
plot(NA,xlim=range(x),ylim=c(0,1),xaxs="i",yaxs="i",axes=F,xlab="",ylab="")
ii <- which.min(beta$V2)
mycol <- c("lightblue","blue","red","darkgreen","orange","purple","magenta")
pal <- colorRampPalette(c("purple","blue","light green", "yellow", "orange", "red","darkred"))
o <-  order(-beta$V2)
mycol <- pal(length(o))
sapply(1:length(o),function(jj){
  ii <- o[jj]
  b0 <- beta$V2[ii]; b1 <- beta$V3[ii];
  y <- plogis(b0+b1*x);
  lines(x,y,col=mycol[jj],lwd=3)
  ii
})
b0 <- beta$V2[ii]; b1 <- beta$V3[ii];
y <- plogis(b0+b1*x);
lines(x,y,col='black',lwd=4,lty=3)
abline(v=10,lty=3)
axis(1)
axis(2,las=1) 
abline(h=1,lty=3)
title(xlab="log2 PWM score")
title(ylab="Predicted proportion bound")
title(main=paste(pwmName,"--",TfNames[pwmName]))

b0 <- beta.new2b[1]; b1 <- beta.new2b[2];
y <- plogis(b0+b1*x);
lines(x,y,col='black',lwd=4,lty=2)

dev.off()