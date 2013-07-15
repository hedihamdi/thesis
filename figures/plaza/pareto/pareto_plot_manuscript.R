#!/usr/bin/Rscript
pdf("plot_pareto.pdf", family="Palatino")
par(mar=c(5,5,5,5))

jitfact=0.5


sizepts=0.25
offset=0.15

pchc=16
pchnc=20
cexc=1
cexnc=1

y0=25
x0=14

#skip=1 to begin reading after line 1
ig=read.table("ovoQ6/roc-enhancers-info.dat")
#legend(2/3*x0,2/3*y0,c("Felsen","Halpern"),fill=c("blue","red"))

name<-ig[,1]
tmot<-ig[,2][name=="ovo_Q6"]
tdet<-ig[,3][name=="ovo_Q6"]
moti<-ig[,4][name=="ovo_Q6"]
nmot<-ig[,5][name=="ovo_Q6"]
auc<-ig[,6][name=="ovo_Q6"]
FP<-jitter(ig[,7][name=="ovo_Q6"],factor=jitfact)
FN<-jitter(ig[,8][name=="ovo_Q6"],factor=jitfact)
xf<-FN
yf<-FP
plot(xf,yf,
     xlab="Faux Negatifs",
     ylab="Faux Positifs",
     main="Pareto Plot",
     pch=pchc,
     xlim=c(0,x0),
     ylim=c(0,y0),
     col="green3",
     cex.main=2,
     cex.axis=2,
     cex.lab=2
     
     )
name<-"c"
auc=round(auc*100)/100
#text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
#text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")


#skip=1 to begin reading after line 1
ig=read.table("TTATGGAA/roc-enhancers-info.dat")
#legend(2/3*x0,2/3*y0,c("Felsen","Halpern"),fill=c("blue","red"))

name<-ig[,1]
tmot<-ig[,2][name=="TTATGGAA"]
tdet<-ig[,3][name=="TTATGGAA"]
moti<-ig[,4][name=="TTATGGAA"]
nmot<-ig[,5][name=="TTATGGAA"]
auc<-ig[,6][name=="TTATGGAA"]
FP<-jitter(ig[,7][name=="TTATGGAA"],factor=jitfact)
FN<-jitter(ig[,8][name=="TTATGGAA"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=pchc,col="yellow")
name<-"c"
auc=round(auc*100)/100
#text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
#text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")



#skip=1 to begin reading after line 1
ig=read.table("rouge-bleu/roc-enhancers-info.dat")
#legend(2/3*x0,2/3*y0,c("Felsen","Halpern"),fill=c("blue","red"))

name<-ig[,1]
tmot<-ig[,2][name=="rouge"]
tdet<-ig[,3][name=="rouge"]
moti<-ig[,4][name=="rouge"]
nmot<-ig[,5][name=="rouge"]
auc<-ig[,6][name=="rouge"]
FP<-jitter(ig[,7][name=="rouge"],factor=jitfact)
FN<-jitter(ig[,8][name=="rouge"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=pchc,xlim=c(0,x0),ylim=c(0,y0),col="red")
name<-"c"
auc=round(auc*100)/100
#text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
#text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")


name<-ig[,1]
tmot<-ig[,2][name=="bleu"]
tdet<-ig[,3][name=="bleu"]
moti<-ig[,4][name=="bleu"]
nmot<-ig[,5][name=="bleu"]
auc<-ig[,6][name=="bleu"]
FP<-jitter(ig[,7][name=="bleu"],factor=jitfact)
FN<-jitter(ig[,8][name=="bleu"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=pchc,xlim=c(0,x0),ylim=c(0,y0),col="blue")
name<-"c"
auc=round(auc*100)/100
#text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
#text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")

cexc2=2
cexnc2=2

name<-ig[,1]
yf<-jitter(ig[,7][name=="combi"],factor=jitfact)
xf<-jitter(ig[,8][name=="combi"],factor=jitfact)
points(xf,yf,col="purple",cex=cexc2,pch=pchnc)

legend("topright",
   leg = c("OvoQ6","svbF7","blue motif","combi svbF7+blue","yellow motif"),
   cex = 1.5,
   pch =pchc,
   #pt.cex = c(cexc,cexnc,cexc,cexnc,cexc,cexnc,cexc2,cexnc2,cexc2,cexnc2),
   col = c("green3","red","blue","purple","yellow"))
