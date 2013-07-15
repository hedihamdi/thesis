#!/usr/bin/Rscript
pdf("plot_pareto.pdf")

jitfact=0.5

sizepts=0.25
offset=0.15

y0=18
x0=14

#skip=1 to begin reading after line 1
ig=read.table("roc-enhancers-info.dat",skip=1)
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
plot(xf,yf,xlab="False Negatives",ylab="False Positives",main="Pareto Plot for RED motif (Enhancers)",pch=46,xlim=c(0,x0),ylim=c(0,y0),col="red")
name<-""
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")

name<-ig[,1]
tmot<-ig[,2][name=="rouge_nocons"]
tdet<-ig[,3][name=="rouge_nocons"]
moti<-ig[,4][name=="rouge_nocons"]
nmot<-ig[,5][name=="rouge_nocons"]
auc<-ig[,6][name=="rouge_nocons"]
FP<-jitter(ig[,7][name=="rouge_nocons"],factor=jitfact)
FN<-jitter(ig[,8][name=="rouge_nocons"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=20,col="pink")
name<-"nocons"
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=offset)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=offset)


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
plot(xf,yf,xlab="False Negatives",ylab="False Positives",main="Pareto Plot for BLUE motif (Enhancers)",pch=46,xlim=c(0,x0),ylim=c(0,y0),col="blue")
name<-""
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")

name<-ig[,1]
tmot<-ig[,2][name=="bleu_nocons"]
tdet<-ig[,3][name=="bleu_nocons"]
moti<-ig[,4][name=="bleu_nocons"]
nmot<-ig[,5][name=="bleu_nocons"]
auc<-ig[,6][name=="bleu_nocons"]
FP<-jitter(ig[,7][name=="bleu_nocons"],factor=jitfact)
FN<-jitter(ig[,8][name=="bleu_nocons"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=20,col="cyan")
name<-"nocons"
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=offset)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=offset)

