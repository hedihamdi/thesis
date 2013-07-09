#!/usr/bin/Rscript
cexmain <- 1.5
cexlab <- 1.5
cexaxis <- 1.5

wb <- c(expression(w[A]), expression(w[C]), expression(w[G]), expression(w[T]))

# INFO
resfile<-"sampling-info.dat"
t<-read.table(resfile)
pseud<-NULL # alpha,beta,beta,alpha
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
for (i in 1:4){
   wmeaninde[i] <- ( Nb[i] + pseud[i] )/( sum(Nb) + sum(pseud) )
}


#PLOTS
pch<-19
cex<-0.7
lwd <- 2
pdf("plotConvergence.pdf",family="Palatino")
par(mfrow=c(2,2), mar = c(5.1,5.1,2.1,2.1))


# SAMPLE POINTS
resfile<-"sampling-w.dat"
t<-read.table(resfile)
Ntot<-length(t$V1)
#Wmean<-as.matrix(rbind(winit,t[,5:8]))
Wmean<-as.matrix(rbind(wmeaninde,t[,5:8]))
colnames(Wmean)<-c("wa mean","wc mean","wg mean","wt mean")

resfile<-"sampling-w-opti.dat"
t<-read.table(resfile)
Wopti<-as.matrix(t[,1:4])
colnames(Wopti)<-c("wa opti","wc opti","wg opti","wt opti")

# MEAN
Ncv<-0
#iter <- seq(1, 1500, 10)
iter <- floor(10^(seq(0, log10(nrow(Wmean)), length.out=100)))
names <- c('wa','wc','wg','wt')
for (i in 1:4){
   name<-wb[i]
   plot(
        iter,
        Wmean[iter,i],
        xlab="Itérations",
        ylab=name,
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis,
        pch=pch,
        cex=cex,
        ylim=c(0,1),
        log="x"
        )
   points(Wopti[,i], pch = pch, cex = cex, col = "red")
   abline(h=Wmean[nrow(Wmean),i], col = "black", lty = 2, lwd = lwd)
   abline(h=Wopti[nrow(Wopti),i], col = "red", lty = 2, lwd = lwd)
   abline(h=wmeaninde[i], col = "purple", lty = 2, lwd = lwd)
   legend("topright",leg = c('MCMC','Gradient','Inde'), pch=c(pch,pch,-1),lty=c(-1,-1,2),col = c('black','red', 'purple'), cex = 1.2,pt.cex=cex,lwd=c(-1,-1,lwd))
}
for (i in 1:4){
   name<-wb[i]
   plot(
        iter,
        Wmean[iter,i],
        xlab="Itérations",
        ylab=name,
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis,
        pch=pch,
        cex=cex,
        ylim=c(min(Wmean[,i],Wopti[,i],wmeaninde[i]),max(Wmean[,i],Wopti[,i],wmeaninde[i])),
        log="x"
        )
   points(Wopti[,i], pch = pch, cex = cex, col = "red")
   abline(h=Wmean[nrow(Wmean),i], col = "black", lty = 2, lwd = lwd)
   abline(h=Wopti[nrow(Wopti),i], col = "red", lty = 2, lwd = lwd)
   abline(h=wmeaninde[i], col = "purple", lty = 2, lwd = lwd)
   legend("topright",leg = c('Mean','Max','Inde'), lty = 2,col = c('black','red', 'purple'), lwd = lwd, cex = 1.5)
}

