
baseFolder='initial/combo/'

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
thresh<-cargs[1];

thresh <- as.numeric(thresh)
cat(thresh)

# two-column file with PWM ID and factor names
aux <- read.table("factorNames.txt",sep="\t",as.is=T)
TfNames <- aux$V2
names(TfNames) <- aux$V1

aux <- read.table("dnase/numReads.txt",sep="\t",as.is=T,quote="")
Nreads <- aux$V2
names(Nreads) <- sub(".x8b.bz2.Qsub.e","",aux$V1)

aux <- scan(pipe(paste("find ",baseFolder," -name '*.lmat.txt.gz'",sep="")),what=character(0))
aux <- sub(baseFolder,"",aux)
pwmNames <- sub(".lmat.txt.gz","",aux)

ii <- "M00256"
aux <- paste(baseFolder,"/",ii,'.beta.txt.gz',sep="")
aux <- (read.table(aux,as.is=T,sep="\t"))
sNames <- aux$V1

system.time(
betaZ1 <- t(sapply(pwmNames,function(ii){
  aux <- paste(baseFolder,"/",ii,'.beta.txt.gz',sep="")
  aux <- (read.table(aux,as.is=T,sep="\t"))
  ret <- aux$V5
  names(ret) <- aux$V1
  ret[sNames]
}))
)
str(betaZ1)
sum(is.na(betaZ1))
which(is.na(betaZ1),arr.ind=T)

## RPR: I had to do this to fix manually the samples that did not work. 
## Does for which the which(is.na returned errors)
## row 213: M00289,  
## col 396: EpiUwRmapDNaseFetal.Kidney.Right.SRX040397
pwmNames <- pwmNames[-213]
sNames <- sNames[-396]
Nreads <- Nreads[sNames]

betaZ1 <- betaZ1[pwmNames,sNames]
which(is.na(betaZ1),arr.ind=T)

zrange <- range(betaZ1,na.rm=T)
Nact <- rowSums(betaZ1>6,na.rm=T)

## This is where we would apply whatever criteria we'd like to use
## to call a TF "active" in a tissue. Could be blanket rule (zscore > 5),
## or more complex (determine cutoff where x% of TFs are "active", etc)
actTfs <- betaZ1 > thresh
tmpFolder <- "/wsu/tmp/"
junk <- sapply(sNames, function(name,betas) {
	aux <- betas[,name]
	aux <- names(which(aux))
	fileRes=paste(tmpFolder,"/",name,".pwms.txt",sep="")
	cat(file=fileRes,aux[order(aux)],sep="\n")
	},
	actTfs)
rm(junk)
