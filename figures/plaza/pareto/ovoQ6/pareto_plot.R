#!/usr/bin/Rscript
pdf("plot_pareto.pdf")

jitfact=0.5

y0=72
x0=30

sizepts=0.25
offset=0.15

y0=18
x0=14

#skip=1 to begin reading after line 1
ig=read.table("roc-enhancers-info.dat",skip=1)
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
plot(xf,yf,xlab="False Negatives",ylab="False Positives",main="Pareto Plot for OVOQ6 (Enhancers)",pch=46,xlim=c(0,x0),ylim=c(0,y0),col="darkgreen")
name<-"c"
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=0.15)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=0.15)
abline(y0,-y0/x0,lty="dashed")

name<-ig[,1]
tmot<-ig[,2][name=="ovo_Q6_no_cons"]
tdet<-ig[,3][name=="ovo_Q6_no_cons"]
moti<-ig[,4][name=="ovo_Q6_no_cons"]
nmot<-ig[,5][name=="ovo_Q6_no_cons"]
auc<-ig[,6][name=="ovo_Q6_no_cons"]
FP<-jitter(ig[,7][name=="ovo_Q6_no_cons"],factor=jitfact)
FN<-jitter(ig[,8][name=="ovo_Q6_no_cons"],factor=jitfact)
xf<-FN
yf<-FP
points(xf,yf,pch=20,col="green3",cex=0.5)
name<-"nc"
auc=round(auc*100)/100
text(FN,FP,pos=1,labels=paste(name,":",tmot,",t",tdet,",",nmot,"x",moti,sep=""),cex=sizepts,offset=offset)
text(FN,FP,pos=3,labels=auc,cex=sizepts,offset=offset)


