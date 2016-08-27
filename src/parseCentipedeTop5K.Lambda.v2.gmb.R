
baseFolder = 'initial/'
outFolder  = "initial/combo"

cargs<-commandArgs(trail=TRUE);
if(length(cargs)>=1)
pwmName<-cargs[1];

cat(pwmName,"\n");

aux <- scan(pipe(paste("ls ",baseFolder,"*/",pwmName,".lambda.txt",sep="")),what=character(0))
aux <- sub(baseFolder,"",aux)
snames <- sub("/.*","",aux)

system.time(
beta <- t(sapply(snames,function(ii){
  aux <- paste(baseFolder,ii,"/",pwmName,'.Beta.txt',sep="")
  as.matrix(read.table(aux,as.is=T,sep="\t")[-1])
}))
)
beta <- beta[,-ncol(beta)]

system.time(
  lmat <-  t(sapply(snames,function(ii){
    aux <- paste(baseFolder,ii,"/",pwmName,'.lambda.txt',sep="")
    as.matrix(read.table(aux,as.is=T,sep="\t"))
  }))
  )
lmat <- as.matrix(lmat[,-ncol(lmat)]);
row.names(lmat) <- snames;
S <- ncol(lmat)

system.time(
  pmat <-  sapply(snames,function(ii){
    aux <- paste(baseFolder,ii,"/",pwmName,'.PostPr.gz',sep="")
    cat(aux,"\n")
    as.matrix(read.table(aux,as.is=T,sep="\t"))
  }))
L <- nrow(pmat)


## Output aggregate data from all the CENTIPEDE tissues
resFile=paste(outFolder,"/",pwmName,".lmat.txt.gz",sep="")
write.table(lmat,gzfile(resFile),quote=FALSE,col.names=F,row.names=T,sep="\t")

resFile=paste(outFolder,"/",pwmName,".beta.txt.gz",sep="")
write.table(beta,gzfile(resFile),quote=FALSE,col.names=F,row.names=T,sep="\t")

resFile=paste(outFolder,"/",pwmName,".pmat.txt.gz",sep="")
write.table(pmat,gzfile(resFile),quote=FALSE,col.names=T,row.names=F,sep="\t")


## Plot the lambda profiles for each of the tissues
probs <- seq(0,1,0.2)
lab <- paste("(",probs[-length(probs)]*100,",",probs[-1]*100,"]%",sep="")
steps <- c("blue4", "blue","green", "red4")
pal <- colorRampPalette(c("purple","blue","darkgreen","gold","red","darkred"))
mycols <- pal(length(lab))
W <- 100;
Mlen <- S/2-W*2;
tmpFolder <- Sys.getenv("TMPDIR","/tmp")
Sys.setenv("DISPLAY"="localhost:10.0")
sapply(c("OpenChrom","UwDnase","UwDgf","EpiUwRmap"),function(myh){
  sel <- grep(myh,snames)
  aux <- -beta[sel,1]
  breaks <- quantile(aux,probs=probs)
  cc <- cut(aux,breaks,labels=lab)
  fname <- paste(tmpFolder,"/",myh,".",pwmName,".png",sep="")
  png(fname)
  par(cex=1.0)
  par(mar=c(2,2,2,2))
  plot(NA,xlim=c(0,S/2),ylim=c(0,max(lmat)/2),xlab="",ylab="",axes="F")
  abline(h=1/S,lty=2)
  sapply(1:length(lab),function(jj){
    sel2 <- sel[cc==lab[jj]];##intersect(grep(myh,snames),which(cc==lab[jj]))
    sapply(sel2,function(ii){
      points(lmat[ii,1:(S/2)]*0.5+lmat[ii,(1:(S/2))+S/2]*0.5,col=mycols[jj],pch='.',cex=2.0)
      0
    })
    lines(colMeans(lmat[sel2,1:(S/2)]*0.5+lmat[sel2,(1:(S/2))+S/2]*0.5),col=mycols[jj],lwd=2,pch='.',cex=2.0)
    length(sel2)
  })
  legend("topright",lab,lwd=2,lty=1,col=mycols,title=expression(paste(beta[1]," Quantiles")))
  axis(1, at = seq(1, W, len = 3),
       labels = -(W + 1 - seq(1, W + 1, len = 3)), padj = -0.5, tck = -0.01)
  axis(1, at = W + Mlen + seq(1, W, len = 3),
       labels = seq(0, W, len = 3), padj = -0.5, tck = -0.01)
  axis(1, at = W + seq(1, Mlen), labels = NA, padj = -0.5, 
       tck = +0.01, col = "purple4")
  axis(2, padj = 0.5)
  abline(v = c(W, W + Mlen + 1), lty = 2)
  title(paste(myh,pwmName))
  dev.off()
  system(paste("scp -r ",fname,outFolder))
  fname
})

