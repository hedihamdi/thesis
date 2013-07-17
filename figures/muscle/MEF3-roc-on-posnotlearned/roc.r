#!/usr/bin/Rscript

files <- c('roc-known-woevol-coordsposavailable-8.dat','evol/roc-h6.dat')
legend<-c("MEF3 reference","MEF3 avec evolution")
cols <- c("black",'red')

pdf("roc.pdf", family="Palatino")
par(mar=c(5,5,5,5))

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
             main = "Courbe ROC de MEF3",
             cex.main = 1.5,
             cex.axis = 1.5,
             cex.lab = 1.5,
             xlim = c(0,1),
            ylim = c(0,1)
            )
    else
        points(FP,TP,lwd=3,type="o",col=cols[i])
    detach(inp)
}
abline(0,1,lty="dashed")
legend("bottomright",legend,fill=cols, cex=1.3)


#for (t in 5:9){

    #file <- paste('roc-known-woevol-coordsposavailable-',t,'.dat',sep='')
    #inp<-scan(file,list(thres=0,FP=0,TP=0))
    #attach(inp)
    #lines(FP,TP,lwd=1,col="grey")
    #detach(inp)
#}
    #file <- paste('roc-known-woevol-coordsposavailable-8.dat',sep='')
    #inp<-scan(file,list(thres=0,FP=0,TP=0))
    #attach(inp)
     #points(FP,TP,lwd=3,type="o",col="black")
    #detach(inp)

#for (t in c('f5','f6','f7','h5','h6','h7')){

    #file <- paste('evol/roc-',t,'.dat',sep='')
    #inp<-scan(file,list(thres=0,FP=0,TP=0))
    #attach(inp)
    #lines(FP,TP,lwd=1,col="pink")
    #detach(inp)
#}
#file <- 'evol/roc-h6.dat'
#inp<-scan(file,list(thres=0,FP=0,TP=0))
#attach(inp)
#points(FP,TP,lwd=3,type="o",col="red")
#detach(inp)

#q()
