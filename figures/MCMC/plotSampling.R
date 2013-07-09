#!/usr/bin/Rscript


cexmain <- 1.5
cexlab <- 1.5
cexaxis <- 1.5


wb <- c(expression(w[A]), expression(w[C]), expression(w[G]), expression(w[T]))
wblog <- c(expression(log[10](C[A](t))), 
           expression(log[10](C[C](t))), 
           expression(log[10](C[G](t))), 
           expression(log[10](C[T](t))))

bases <- c("A","C","T","G")

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

pdf("plotDirichlet.pdf",family = "Palatino")
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))


# DIRICHLET
resfile<-"sampling-dirichlet.dat"
t<-read.table(resfile)
W<-as.matrix(t[,1:4])
for (i in 1:4){
   hist(W[,i],
        xlab=wb[i],
        ylab="Nombre d'éléments",
        main="",#paste("Distribution de Dirichlet pour",parse(text=wb[i]),"\n winit=",winit[i]),
        breaks=seq(0,1,0.02),
        xlim=c(0,1),
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
   abline(v=winit[i],col="red",lwd=3)
}
#par(mfrow=c(1,1), oma = c(0,2,2,0))
#mtext(paste("Ndir=",Ndir,",",consInst,"conserved instances"), 3, outer = TRUE,font = 2)
#mtext("Distribution de Dirichlet", 3, outer = TRUE,font = 2, cex= cexmain)

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
#par(mfrow=c(1,1), oma = c(0,2,2,0))
#mtext(paste("Acceptation rate=",acceptRate), 3, outer = TRUE,font = 2)
#mtext("Echantillonnage", 3, outer = TRUE,font = 2, cex = cexmain)

## LogLIKELIHOOD
#resfile<-"sampling-loglikelihood-w-dirichlet.dat"
#t<-na.omit(read.table(resfile))
#fmean<-t$V1
#fmax<-t$V2
#par(mfrow=c(2,2), oma = c(0,0,3,0))
#for (i in 1:4){
   #name<-colnames(W)[i]
   #norm<-log(sum(exp(min(fmean[1:1000])-fmean[1:1000])))-min(fmean[1:1000])
   #plot(W[1:1000,i],fmean[1:1000]-norm,xlab=name,ylab=paste("Loglikelihood"),main=paste("Loglikelihood values sampled"),xlim=c(0,1),pch=pch,cex=cex)
   #norm<-log(sum(exp(min(fmax[1:1000])-fmax[1:1000])))-min(fmax[1:1000])
   #points(W[1:1000,i],fmax[1:1000]-norm,pch=pch,cex=cex,col="green")
   #abline(v=winit[i],col="red",lwd=3)
   #abline(v=Wmean[Ntot,i],col="black",lwd=3)
   #abline(v=wopti[i],col="green",lwd=3)
   #abline(v=wmeaninde[i],col="yellow",lwd=3)
   #if (Wmean[Ntot,i]>0.5) pos <- "topleft"
   #else pos <- "topright"
   #legend(pos,c("Likelihood for mean","Likelihood for max","w init","w mean","w max (SD)","w mean inde"),col=c("black","green","red","black","green","yellow"),pch=c(pch,pch,-1,-1,-1,-1),lwd=c(-1,-1,3,3,3,3),lty=c(-1,-1,1,1,1,1))
#}
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Acceptation rate=",acceptRate), 3, outer = TRUE,font = 2)

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

## WINDOW
#icW<-NULL
#for (d in seq(1,200,5)){
   #ws<-NULL
   #xs<-1:((Ntot-d)/d)
   #for (i in xs){
      #beg<-i*d
      #sto<-(i+1)*d
      #ws<-rbind(ws,colMeans(W[beg:sto,]))
   #}
   #icW<-rbind(icW,1.96*sqrt(apply(ws,2,var)/length(xs)))
#}
#par(mfrow=c(2,2), oma = c(0,0,3,0))
#for(i in 1:4){
   #name<-colnames(W)[i]
   #plot(seq(1,200,5),icW[,i],main=paste("Variance of",name,"averaged over windows of size d"),xlab="d",ylab="IC(d) ~ sqrt(var(w)/d)",pch=pch,cex=cex)
#}
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext("Another way to do the same thing: check independance of drawings using Central Limit Theorem", 3, outer = TRUE,font = 2)

nsample<-floor(max(taus))
nsample<-30
iter<-nsample*seq(Ntot/nsample)

if (nsample>100){
   print("Problems with nsample, using value of 100")
   nsample <- 100
}

## DKL
#par(mfrow=c(1,1), oma = c(0,0,3,0))
Niter<-length(iter)
#dkl<-rowSums(Wmean[iter[2:Niter],]*log2(Wmean[iter[2:Niter],]/Wmean[iter[1:(Niter-1)],]))
#dkl[ dkl<=0 ] <- min(dkl[dkl>0])
#plot(log10(nsample*(1:(Niter-1))),log10(dkl),xlab="Number of iteration (log10 scale)",ylab="dkl (bits, log10 scale)",pch=pch,cex=cex,main="Kullback-Leibler Divergence to previous mean")
#g<-glm(log10(dkl)~log10(nsample*(1:(Niter-1))))
#abline(g$coeff)

## CUTOFF
#xs<-nsample*(1:(Niter-1))
#dklcut <- 1e-5
#cutoff <- min(xs[which(dkl<dklcut)])
#abline(v=log10(cutoff),col="blue",lwd=3)

#legend("topright",c(paste("Dkl ~ 1/N^",-1/10*floor(10*g$coeff[2])),paste("cutoff DKL=",dklcut,"N=",cutoff)),col=c("black","blue"),lwd=c(1,3))
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)

## Euclidean distance
#par(mfrow=c(1,1), oma = c(0,0,3,0))
#Niter<-length(iter)
#dist<-rowSums((Wmean[iter[2:Niter],]-Wmean[iter[1:(Niter-1)],])^2)
#dist[ dist==0 ] <- min(dist[dist>0])
#plot(log10(nsample*(1:(Niter-1))),log10(dist),xlab="Number of iteration (log10 scale)",ylab="Euclidean distance (log10 scale)",pch=pch,cex=cex,main="Euclidean distance to previous mean")
#estim <- rowSums(Wvar[nsample*(1:(Niter-1)),])*(nsample/(nsample+nsample*(1:(Niter-1))))^2
#lines(log10(nsample*(1:(Niter-1))),log10(estim),col="red")
#g<-glm(log10(dist)~log10(nsample*(1:(Niter-1))))
#abline(g$coeff)
#legend("topright",c(paste("Dist ~ 1/N^",-1/10*floor(10*g$coeff[2])),"estimation sum(Var(w))*(n/(n+N))^2"),col=c("black","red"),lty=1)
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)

pdf("plotCVmeanvar.pdf",family = "Palatino")
# MEAN
#par(mfrow=c(2,2), oma = c(0,0,3,0))
par(mfrow=c(2,2), oma = c(0,2,3,0))
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
   plot(iter,
        samplemean,
        log="x",
        xlab="Itération",
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

   #ncv <- 10
   #itercv <- 0
   #for (j in (ncv+1):Niter){
      ##if ( sum(abs(Wmean[iter[(j-ncv):(j-1)],i]-Wmean[iter[j],i])<1.96*sqrt(wvarinde[i]/((j-ncv):(j-1)))) == ncv ){
      #if ( sum(abs(samplemean[(j-ncv):(j-1)]-samplemean[j])<1.96*sqrt(wvarinde[i]/((j-ncv):(j-1)))) == ncv ){
         #itercv <- j
         #break
      #}
   #}
  #abline(v=iter[itercv],col="yellow",lwd=3)
   #Ncv<-max(Ncv,iter[itercv])

   #if (Wmean[Ntot,i]>0.5) pos <- "bottomright"
   #else pos <- "topright"
   #abline(v=cutoff,col="blue",lwd=3)
   #legend(pos,c("95% confindence with inde variance",paste("cutoff DKL=",dklcut,"N=",cutoff),"cutoff with var inde"),col=c("red","blue","yellow"),lty=c(2,1,1),lwd=c(3,3,3))
}
par(mfrow=c(1,1), oma = c(0,2,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)
mtext("Convergence de la moyenne", 3, outer = TRUE,font = 2,cex=1.5)

## LogLIKELIHOOD
#resfile<-"sampling-loglikelihood-w-dirichlet.dat"
#t<-na.omit(read.table(resfile))
#f<-t$V1
#par(mfrow=c(2,2), oma = c(0,0,3,0))
#for (i in 1:4){
   #name<-colnames(W)[i]
   #plot(W[1:1000,i],f[1:1000],xlab=name,ylab=paste("Loglikelihood"),main=paste("Loglikelihood values sampled"),xlim=c(0,1),pch=pch,cex=cex)
   #abline(v=winit[i],col="red",lwd=3)
   #abline(v=Wmean[Ntot,i],col="black",lwd=3)
   #abline(v=Wmean[cutoff,i],col="blue",lwd=3)
   #abline(v=Wmean[Ncv,i],col="yellow",lwd=3)
   #if (Wmean[Ntot,i]>0.5) pos <- "topleft"
   #else pos <- "topright"
   #legend(pos,c("w init","w mean",paste("w mean at dkl cutoff=",dklcut),"w mean at 95% CI"),col=c("red","black","blue","yellow"),lwd=3)
#}
#par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Acceptation rate=",acceptRate), 3, outer = TRUE,font = 2)




# VAR
par(mfrow=c(2,2), oma = c(0,0,3,0))
for (i in 1:4){
   # i use max value = uniform distribution
   ymax <- max(1/3-1/2^2,wvarinde[i],Wvar[,i])
   name<-colnames(Wvar)[i]
   plot(iter,
        log="x",
        Wvar[iter,i],
        xlab="Itération",
        ylab=wb[i],
        main="",
        pch=pch,
        cex=cex,
        ylim=c(0,ymax),
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
   #abline(v=cutoff,col="blue",lwd=3)
   #abline(v=Ncv,col="yellow",lwd=3)
   #legend(pos,c(paste("cutoff DKL=",dklcut,"N=",cutoff),"cutoff 95% CI with var inde"),col=c("blue","yellow"),lty=c(1,1),lwd=c(3,3))
}
par(mfrow=c(1,1), oma = c(0,0,2,0))
#mtext(paste("Sample step = 5*tau ~",nsample), 3, outer = TRUE,font = 2)
mtext("Convergence de la variance", 3, outer = TRUE,font = 2,cex=1.5)

