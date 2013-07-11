#!/usr/bin/Rscript
cexmain <- 1.5
cexlab <- 1.5
cexaxis <- 1.5
alpha <- 0.05 # precision for convergence

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
wvarinde <- NULL
for (i in 1:4){
   wmeaninde[i] <- ( Nb[i] + pseud[i] )/( sum(Nb) + sum(pseud) )
   wvarinde[i] <- wmeaninde[i]*( (Nb[i]+pseud[i]+1)/(sum(Nb)+sum(pseud)+1) - wmeaninde[i] )
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
W<-as.matrix(t[,1:4])
Wmean<-as.matrix(rbind(wmeaninde,t[,5:8]))
Wvar<-as.matrix(t[,9:12])

nsample<-30
iter<-nsample*seq(Ntot/nsample)
Niter<-length(iter)
sdW <- list()
sdTLC <- list()
hascv <- NULL
offset <- 11
for (i in 1:4){
    sdw<-NULL
    sdtlc <- NULL
    comp <- NULL
    isok <- NULL
    for (j in offset:Niter){
        #mm <- max(abs(Wmean[iter[(j-(offset-1)):(j-1)],i]-Wmean[iter[j],i]))
        #mm <- sqrt(1/offset*sum((Wmean[iter[(j-(offset-1)):(j-1)],i]-Wmean[iter[j],i])^2))
        mm <- sqrt(Wvar[iter[j],i]/j)
        sdw <- c(sdw,mm)
        sdtlc <- c(sdtlc, sqrt(wvarinde[i]/j))
        isok <- c(isok, mm <= alpha*Wmean[iter[j],i])
    }
    sdW[[i]] <- sdw
    sdTLC[[i]] <- sdtlc
    #hascv <- cbind(hascv, sdw <= alpha * sdtlc)
    hascv <- cbind(hascv,isok)
}
inds <- which(rowSums(hascv) == 4)
cv <- inds[1]+offset-1
   
delta <- c(expression(sigma[A]/n),
           expression(sigma[C]/n),
           expression(sigma[G]/n),
           expression(sigma[T]/n)
           )
for (i in 1:4){
    #plot(iter[offset:Niter],
    plot(offset:Niter,
         sdW[[i]],
         #sdTLC[[i]]/Wmean[iter[offset:Niter],i],
         pch=16, 
         log="x",
         xlab = "n",
         ylab = delta[i],
         main="",
         ylim = c(min(sdTLC[[i]],sdW[[i]]),max(sdTLC[[i]],sdW[[i]])),
         xaxt="n",
         cex.main=cexmain,
         cex.lab=cexlab,
         cex.axis=cexaxis
         )
    axis(1, at =c(1,10,100,1000), cex.axis=cexaxis)
    #lines(iter[offset:Niter], sdTLC[[i]], col = 'red')
    #lines(iter[offset:Niter], alpha*Wmean[iter[offset:Niter],i], col = 'grey',lwd=lwd)
    lines(offset:Niter, alpha*Wmean[iter[offset:Niter],i], col = 'grey',lwd=lwd)
    #abline(v=iter[cv], col = 'blue',lwd=lwd)
    #legend("topright", c('Inde','CV'),col=c('red','blue'),lty=1,cex=1.5)
    #legend("topright", c('5% of mean', 'CV'),col=c('grey','blue'),lty=1,cex=1.2,lwd=lwd)
    legend("topright", c('5% of mean'),col=c('grey'),lty=1,cex=1.2,lwd=lwd)
}

dev.off()


# DIRICHLET
pdf("plotDirichlet.pdf",family = "Palatino")
par(mfrow=c(2,2), mar = c(5.1,5.1,4.1,2.1))

resfile<-"sampling-dirichlet.dat"
t<-read.table(resfile)
Wdir<-as.matrix(t[,1:4])
for (i in 1:4){
   hist(Wdir[,i],
        xlab=wb[i],
        ylab="Nombre d'éléments",
        main="",#paste("Distribution de Dirichlet pour",parse(text=wb[i]),"\n winit=",winit[i]),
        breaks=seq(0,1,0.02),
        xlim=c(0,1),
        cex.main = cexmain,
        cex.lab = cexlab,
        cex.axis = cexaxis
        )
   #abline(v=winit[i],col="blue",lwd=3)
   #abline(v=wmeaninde[i],col="purple",lwd=3)
   abline(v=Wmean[iter[cv],i],col="black",lwd=3, lty =2)
}
dev.off()

pdf("plotConvergenceCompare.pdf",family="Palatino")
par(mfrow=c(2,2), mar = c(5.1,5.1,2.1,2.1))

resfile<-"sampling-w-opti.dat"
t<-read.table(resfile)
Wopti<-as.matrix(t[,1:4])
colnames(Wopti)<-c("wa opti","wc opti","wg opti","wt opti")

# MEAN
Ncv<-0
#iter <- seq(1, 1500, 10)
#xs <- floor(10^(seq(0, log10(nrow(Wmean)), length.out=100)))
xs <- floor(10^(seq(0, log10(iter[cv]), length.out=100)))
names <- c('wa','wc','wg','wt')
for (i in 1:4){
   name<-wb[i]
   plot(
        xs,
        Wmean[xs,i],
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
   #abline(h=Wmean[nrow(Wmean),i], col = "black", lty = 2, lwd = lwd)
   abline(h=Wmean[iter[cv],i], col = "black", lty = 2, lwd = lwd)
   abline(h=Wopti[nrow(Wopti),i], col = "red", lty = 2, lwd = lwd)
   #abline(h=wmeaninde[i], col = "purple", lty = 2, lwd = lwd)
   #abline(v=iter[cv], col = "blue", lty = 1, lwd = lwd)
   #legend("topright",leg = c('MCMC','Gradient','Inde','CV'), pch=c(pch,pch,-1,-1),lty=c(-1,-1,2,1),col = c('black','red', 'purple','blue'), cex = 1.2,pt.cex=cex,lwd=c(-1,-1,lwd,lwd))
   #legend("topright",leg = c('Mean','Max','Inde'), lty = 2,col = c('black','red', 'purple'), lwd = lwd, cex = 1.2)
   legend("topright",leg = c('MCMC','Gradient'), pch=pch,col = c('black','red'), cex = 1.2, pt.cex=cex)
}
#for (i in 1:4){
   #name<-wb[i]
   #plot(
        #xs,
        #Wmean[xs,i],
        #xlab="Itérations",
        #ylab=name,
        #cex.main = cexmain,
        #cex.lab = cexlab,
        #cex.axis = cexaxis,
        #pch=pch,
        #cex=cex,
        #ylim=c(min(Wmean[,i],Wopti[,i],wmeaninde[i]),max(Wmean[,i],Wopti[,i],wmeaninde[i])),
        #log="x"
        #)
   #points(Wopti[,i], pch = pch, cex = cex, col = "red")
   ##abline(h=Wmean[nrow(Wmean),i], col = "black", lty = 2, lwd = lwd)
   #abline(h=Wmean[iter[cv],i], col = "black", lty = 2, lwd = lwd)
   #abline(h=Wopti[nrow(Wopti),i], col = "red", lty = 2, lwd = lwd)
   #abline(h=wmeaninde[i], col = "purple", lty = 2, lwd = lwd)
   #legend("topright",leg = c('Mean','Max','Inde'), lty = 2,col = c('black','red', 'purple'), lwd = lwd, cex = 1.2)
#}

