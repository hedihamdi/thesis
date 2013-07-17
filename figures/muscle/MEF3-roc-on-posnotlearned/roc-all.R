#!/usr/bin/Rscript
pdf("roc-all.pdf", family="Palatino")
par(mar=c(5,5,5,5))

filesevol <- list.files(path='evol',pattern='*.dat',full.names=TRUE)
#filesevol <- 'evol/roc-h6.dat'
filesref <- list.files(path='.',pattern='*-[0-9].dat',full.names=FALSE)
#legend<-c("MEF3 reference","MEF3 avec evolution (h6)")

files <- filesref
cols <- rainbow(length(files))
for (i in 1:length(files)){
    inp<-scan(files[[i]],list(thres=0,FP=0,TP=0))
    attach(inp)
    if (i==1)
        plot(FP,
             TP,
             lwd = 3,
             type = "o",
             col = cols[i],
             xlab = "Faux Positifs",
             ylab = "Vrais Positifs", 
             main = "Sans evolution",
             cex.main = 1.5,
             cex.axis = 1.5,
             cex.lab = 1.5,
             xlim = c(0,1),
            ylim = c(0,1)
            )
    else
        points(FP,TP,lwd=2,type="o",col=cols[i])
    detach(inp)
}
abline(0,1,lty="dashed")
legend <- files
legend("bottomright",legend,fill=cols, cex=1)


files <- filesevol
cols <- rainbow(length(files))
for (i in 1:length(files)){
    inp<-scan(files[[i]],list(thres=0,FP=0,TP=0))
    attach(inp)
    if (i==1)
        plot(FP,
             TP,
             lwd = 3,
             type = "o",
             col = cols[i],
             xlab = "Faux Positifs",
             ylab = "Vrais Positifs", 
             main = "Avec evolution",
             cex.main = 1.5,
             cex.axis = 1.5,
             cex.lab = 1.5,
             xlim = c(0,1),
            ylim = c(0,1)
            )
    else
        points(FP,TP,lwd=2,type="o",col=cols[i])
    detach(inp)
}
abline(0,1,lty="dashed")
legend <- files
legend("bottomright",legend,fill=cols, cex=1)
q()
