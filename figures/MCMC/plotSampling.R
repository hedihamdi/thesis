#!/usr/bin/Rscript


cexmain <- 1.5
cexlab <- 1.5
cexaxis <- 1.5


wb <- c(expression(w[A]), expression(w[C]), expression(w[G]), expression(w[T]))
wblog <- c(expression(log[10](C[A](t))), 
           expression(log[10](C[C](t))), 
           expression(log[10](C[G](t))), 
           expression(log[10](C[T](t))))


# INFO
resfile<-"sampling-info.dat"
t<-read.table(resfile)
pseud<-NULL # alpha,alpha,beta,beta
pseud[c(1,4)] <- t[1,2]
pseud[2:3] <- t[2,2]
consInst <- t[3,2]
Ndir <- t[4,2]
Nb <- t[5:8,2]
winit<-t[9:12,2]
wopti<-t[13:16,2]
acceptRate <- t[17,2]

# Independent species model
wmeaninde <- NULL
wvarinde <- NULL
for (i in 1:4){
   wmeaninde[i] <- ( Nb[i] + pseud[i] )/( sum(Nb) + sum(pseud) )
   wvarinde[i] <- wmeaninde[i]*( (Nb[i]+pseud[i]+1)/(sum(Nb)+sum(pseud)+1) - wmeaninde[i] )
}


#PLOTS
pch<-19
cex<-0.6

# SAMPLE POINTS
resfile<-"sampling-w.dat"
t<-read.table(resfile)
Ntot<-length(t$V1)
W<-as.matrix(t[,1:4])
Wmean<-as.matrix(t[,5:8])
Wvar<-as.matrix(t[,9:12])

pdf("plotSampling.pdf",family = "Palatino")
#par(mfrow=c(2,2), oma = c(0,2,0,0))
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))
for (i in 1:4){
   name<-wb[i]
   plot(W[1:500,i],
        xlab="Nombre d'itérations",
        ylab=name,
        main="",#,paste("Échantillonage de",name),
        pch=pch,
        cex=cex,
        ylim=c(0,1),
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
}
dev.off()

# Step correlation
Corr<-NULL
nmax<-200
for (n in seq(1,nmax)){
   W1<-W[1:(Ntot-n),]
   W2<-W[(n+1):Ntot,]
   corr<-colMeans(W1*W2)-colMeans(W1)*colMeans(W2)
   Corr<-rbind(Corr,corr)
}
   
pdf("plotAutocorr.pdf",family = "Palatino")
#par(mfrow=c(2,4))
#par(mfrow=c(2,2))
#par(mfrow=c(2,2), oma = c(0,2,0,0))
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))
#for (i in 1:4){
   #name<-colnames(W)[i]
   #plot(seq(1,nmax),Corr[,i],ylab=paste(name),xlab="Step interval",main=paste(name,"autocorrelation"),pch=pch,cex=cex,
        #cex.main = cexmain,
        #cex.lab = cexlab,
        #cex.axis = cexaxis
        #)
#}
taus<-NULL
for (i in 1:4){
   name<-wb[i]
   Corr[Corr==0]<-min(Corr[Corr!=0])
   x<-(1:20)
   y<-abs(Corr[x,i])
   g<-glm(log10(y)~x)
   tau<-5*(-1/(g$coeff[2]*log(10)))
   taus<-c(taus,tau)
   
   plot(
        seq(1,nmax),
        log10(abs(Corr[,i])),
        ylab=wblog[i],
        xlab="t",
        pch = 16,
        main="",#paste(name,"autocorrelation"),pch=pch,cex=cex,
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
   abline(g$coeff)

   usr <- par( "usr" )
   x <- 2/3 * usr[1] + 1/3 * usr[2]
   y <- 1/8 * usr[3] + 7/8 * usr[4]
   
   text(x,y,paste("T ~",1/100*floor(100*tau/5)), cex = 1.5)
  # legend("topright",paste("T ~",1/100*floor(100*tau/5)), cex = 2)
}
mtext("Corrélation entre échantillons", 3, outer = TRUE,font = 2, cex = cexmain)
dev.off()

nsample<-floor(max(taus))
nsample<-30
iter<-nsample*seq(Ntot/nsample)

if (nsample>100){
   print("Problems with nsample, using value of 100")
   nsample <- 100
}

Niter<-length(iter)

pdf("plotCVmean.pdf",family = "Palatino")
# MEAN
#par(mfrow=c(2,2), oma = c(0,0,3,0))
#par(mfrow=c(2,2), oma = c(0,2,3,0))
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))
Ncv<-0
for (i in 1:4){
   name<-colnames(Wmean)[i]
   samplemean<-NULL
   for (j in 1:Niter){
      samplemean<-c(samplemean,mean(W[iter[1:j],i]))
   }
   #plot(iter,Wmean[iter,i],xlab="Number of iteration",ylab=name,main=paste(name),pch=pch,cex=cex,ylim=c(0,1))
   #lines(iter,Wmean[Ntot,i]+1.96*sqrt(wvarinde[i]/(iter/nsample)),lty='dashed',lwd=3,col="red")
   #lines(iter,Wmean[Ntot,i]-1.96*sqrt(wvarinde[i]/(iter/nsample)),lty='dashed',lwd=3,col="red")
   plot(1:Niter,
        samplemean,
        log="x",
        xlab="n",
        ylab=wb[i],
        main="",
        pch=pch,
        cex=cex,
        ylim=c(0,1),
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
   lines(iter,samplemean[Niter]+1.96*sqrt(wvarinde[i]/(iter/nsample)),lty='dashed',lwd=3,col="red")
   lines(iter,samplemean[Niter]-1.96*sqrt(wvarinde[i]/(iter/nsample)),lty='dashed',lwd=3,col="red")

}
#par(mfrow=c(1,1), oma = c(0,2,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)
#mtext("Convergence de la moyenne", 3, outer = TRUE,font = 2,cex=1.5)
dev.off()

pdf("plotCVvar.pdf",family = "Palatino")
#par(mfrow=c(2,2), oma = c(0,0,3,0))
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))
# VAR
sb <- c(expression(sigma[A]), expression(sigma[C]), expression(sigma[G]), expression(sigma[T]))
for (i in 1:4){
   # i use max value = uniform distribution
   ymax <- max(1/3-1/2^2,wvarinde[i],Wvar[,i])
   name<-colnames(Wvar)[i]
   plot(1:Niter,
        log="x",
        Wvar[iter,i],
        xlab="n",
        ylab=sb[i],
        main="",
        pch=pch,
        cex=cex,
        ylim=c(0,ymax),
         xaxt="n",
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
    axis(1, at =c(1,10,100,1000), cex.axis=cexaxis)
   #abline(v=cutoff,col="blue",lwd=3)
   #abline(v=Ncv,col="yellow",lwd=3)
   #legend(pos,c(paste("cutoff DKL=",dklcut,"N=",cutoff),"cutoff 95% CI with var inde"),col=c("blue","yellow"),lty=c(1,1),lwd=c(3,3))
}
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)
#mtext("Convergence de la variance", 3, outer = TRUE,font = 2,cex=1.5)

