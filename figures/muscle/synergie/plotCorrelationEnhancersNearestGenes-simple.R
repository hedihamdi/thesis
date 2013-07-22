#!/usr/bin/Rscript

library(gplots)

cex.main <- 2
cex.lab <- 2
cex.axis <- 2
lwd <- 5

pdf("plotCorrelationEnhancersNearestGenes.pdf", family = "Palatino")
par(mar=c(6,6,6,6))

maxgenes <- 10
maxnum <- 10000

correlatediff<-function(chipsfull,plotname){

   #chips <- chipsfull[sample(length(chipsfull[,1]), min(length(chipsfull[,1]), maxnum)),]
   chips <- chipsfull
   FCs_sm_g <- NULL
   FCs_sm_s <- NULL
   FCs_sm_m <- NULL
   FCs_s_g <- NULL
   FCs_m_g <- NULL
   FCs_dyn <- NULL
   for (i in 1:length(chips[,1])){
      
      chip <- chips[i,]
      name <- chip$V1
      chr <- as.character(chip$V2)
      pos <- (chip$V3 + chip$V4) / 2

      data <- dattss[tsschrom == chr,]
  
      neargenes <- data[order(abs(data[,1] - pos)),]
      neargenes <- as.character(neargenes[1 : maxgenes, 3])
      #neargenes <- as.character(neargenes[1 : maxgenes, 4])

      

      nearfibros <- fibros[rownames(fibros) %in% neargenes,]
      
     if (length(nearfibros) <= 8) next

      FCsm_g <- 2^(nearfibros[,8] - nearfibros[,5])
      FCsm_s <- 2^(nearfibros[,8] - nearfibros[,6])
      FCsm_m <- 2^(nearfibros[,8] - nearfibros[,7])
      FCs_g <- 2^(nearfibros[,6] - nearfibros[,5])
      FCm_g <- 2^(nearfibros[,7] - nearfibros[,5])
      FCdyn <- 2^(nearfibros[,8] - nearfibros[,4])
      FCs <- NULL
      for (neargene in neargenes){
         ind <- match(neargene,rownames(nearfibros))
         if (!is.na(ind))
            FCs <- rbind(FCs, c(FCsm_g[ind], FCsm_s[ind], FCsm_m[ind], FCs_g[ind], FCm_g[ind], FCdyn[ind]))
         else 
            FCs <- rbind(FCs, rep(NA,6))
      }
      
      FCs_sm_g <- cbind(FCs_sm_g, FCs[,1])
      FCs_sm_s <- cbind(FCs_sm_s, FCs[,2])
      FCs_sm_m <- cbind(FCs_sm_m, FCs[,3])
      FCs_s_g <- cbind(FCs_s_g, FCs[,4])
      FCs_m_g <- cbind(FCs_m_g, FCs[,5])
      FCs_dyn <- cbind(FCs_dyn, FCs[,6])
   
   }
   
   meanFCs_sm_g <- rowMeans(FCs_sm_g, na.rm = TRUE )
   meanFCs_sm_s <- rowMeans(FCs_sm_s, na.rm = TRUE )
   meanFCs_sm_m <- rowMeans(FCs_sm_m, na.rm = TRUE )
   meanFCs_s_g <- rowMeans(FCs_s_g, na.rm = TRUE )
   meanFCs_m_g <- rowMeans(FCs_m_g, na.rm = TRUE )
   meanFCs_dyn <- rowMeans(FCs_dyn, na.rm = TRUE )
   sds <- NULL
   sds <- cbind(sds, apply(FCs_sm_g,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   sds <- cbind(sds, apply(FCs_sm_s,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   sds <- cbind(sds, apply(FCs_sm_m,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   sds <- cbind(sds, apply(FCs_s_g,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   sds <- cbind(sds, apply(FCs_m_g,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   sds <- cbind(sds, apply(FCs_dyn,1,std <- function(x){ sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))}))
   ymin <- min(meanFCs_sm_g, meanFCs_sm_s, meanFCs_dyn, meanFCs_s_g, meanFCs_m_g, meanFCs_dyn) - max(sds)
   ymax <- max(meanFCs_sm_g, meanFCs_sm_s, meanFCs_dyn, meanFCs_s_g, meanFCs_m_g, meanFCs_dyn) + max(sds)
   
   FCsm_g <- 2^(fibros[,4] - fibros[,1])
   FCsm_s <- 2^(fibros[,4] - fibros[,2])
   FCsm_m <- 2^(fibros[,4] - fibros[,3])
   FCs_g <- 2^(fibros[,2] - fibros[,1])
   FCm_g <- 2^(fibros[,3] - fibros[,1])
   FCdyn <- 2^(fibros[,8] - fibros[,4])

      
   x <- c(1:maxgenes)
   plot(x,
        meanFCs_sm_g,
        lwd=5,
        col="purple",
        type="l",
        ylim=c(ymin, ymax),
        #main = plotname,
        xlab="Rang du gene",
        ylab="Fold Change",
        xlim=c(1,10), 
        cex.axis=cex.axis, 
        cex.lab=cex.lab, 
        cex.main=cex.main
        )
   plotCI(x,
          meanFCs_sm_g,
          lwd=1,
          col="purple",
          uiw=sds[,1], 
          add=TRUE,
          gap=0, 
          sfrac=0.003, 
          cex=0.1
          )
   abline(h = mean(FCsm_g), 
          lty = 2, 
          col = "purple")
   
   lines(x,
         meanFCs_s_g,
         lwd=5,
         col="red"
         )
   plotCI(x,
          meanFCs_s_g,
          lwd=1,
          col="red",
          uiw=sds[,4], 
          add=TRUE,
          gap=0, 
          sfrac=0.003, 
          cex=0.1
          )
   abline(h = mean(FCs_g), 
          lty = 2, 
          col = "red"
          )

   lines(x,
         meanFCs_m_g,
         lwd=5,
         col="blue"
         )
   plotCI(x,
          meanFCs_m_g,
          lwd=1,
          col="blue",
          uiw=sds[,5], 
          add=TRUE,
          gap=0, 
          sfrac=0.003, 
          cex=0.1
          )
   abline(h = mean(FCm_g), 
          lty = 2, 
          col = "blue"
          )

   leg <- c("Six1,4+MyoD VS GFP","Six1,4 VS GFP","MyoD VS GFP")
   legend("topright",
          leg=leg,
          col=c("purple","red","blue"),
          lwd=rep(5,3), 
          cex=1.5
          )

}

# TSS data
#dattss<-read.table('/Users/Marc/Desktop/work/these/projects/cochin/pascal/fibros-six+myod-synergy/affy+qPCR/affy/files/genes-n-strand-all-TSS.dat',header=TRUE)
dattss <- read.table('genes-n-strand-only-protein-coding-TSS.dat')
dattss <- read.table('genes-n-strand-all-TSS.dat')
tssgene <- dattss[,3]
tsscoord <- dattss[,1]
tsschrom <- as.character(dattss[,4])
#dattss <- read.table('/Users/Marc/Desktop/work/these/programs/imogene/build/install/share/imogene/eutherian/annot/TSS-coord.dat')
#tssgene <- as.character(dattss[,4])
#tsscoord <- as.numeric(dattss[,1])
#tsschrom <- as.character(dattss[,3])

# fibroblasts
fibroswt<-read.table('/Users/Marc/Desktop/work/these/projects/cochin/pascal/fibros-six+myod-synergy/affy+qPCR/matrice_filtered.txt',header=TRUE)
matfibrowtdiff<-as.matrix(fibroswt[,c(5,10,6,11)])
matfibrowtprolif<-as.matrix(fibroswt[,c(2,7,3,8)])
fibros<-cbind(matfibrowtprolif,matfibrowtdiff)
rownames(fibros)<-fibroswt[,1]
FCsdyn<-2^(fibros[,8]-fibros[,4])

# C2C12 data for comparison
c2c12_wt<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/WT-MB.txt", row.names=1,header=TRUE)
matc2c12_wt<-as.matrix(c2c12_wt)
c2c12_wt<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/WT-T0.txt", row.names=1,header=TRUE)
matc2c12_wt<-cbind(matc2c12_wt,c2c12_wt[match(rownames(c2c12_wt),rownames(matc2c12_wt)),])
c2c12_wt<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/WT-T24.txt",row.names=1,header=TRUE)
matc2c12_wt<-cbind(matc2c12_wt,c2c12_wt[match(rownames(c2c12_wt),rownames(matc2c12_wt)),])
c2c12_wt<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/WT-MT.txt", row.names=1,header=TRUE)
matc2c12_wt<-cbind(matc2c12_wt,c2c12_wt[match(rownames(c2c12_wt),rownames(matc2c12_wt)),])
colnames(matc2c12_wt)<-c("WT MB","WT T0","WT T24","WT MT")

c2c12_ko<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/KO-SIX1.txt", row.names=1,header=TRUE)
matc2c12_ko<-as.matrix(c2c12_ko)
c2c12_ko<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/KO-SIX4.txt", row.names=1,header=TRUE)
matc2c12_ko<-cbind(matc2c12_ko,c2c12_ko[match(rownames(c2c12_ko),rownames(matc2c12_ko)),])
c2c12_ko<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/KO-SIX14.txt", row.names=1,header=TRUE)
matc2c12_ko<-cbind(matc2c12_ko,c2c12_ko[match(rownames(c2c12_ko),rownames(matc2c12_ko)),])
c2c12_ko<-read.table("/Users/Marc/Desktop/work/these/data/affy/six14-blais/KO-MYOG.txt", row.names=1,header=TRUE)
matc2c12_ko<-cbind(matc2c12_ko,c2c12_ko[match(rownames(c2c12_ko),rownames(matc2c12_ko)),])
colnames(matc2c12_ko)<-c("KO SIX1","KO SIX4","KO SIX1,4","KO MYOG")

# MYOD + MEF3
datchip_myod_mef3<-read.table("/Users/Marc/Desktop/work/these/projects/cochin/pascal/fibros-six+myod-synergy/affy+qPCR/affy/files/myod+mef3-coords.dat")

# the best six+myod genes
fcpos<-1.5 # synergy
fcneg<-1.3 # no synergy

print('myod+mef3')
chips <- datchip_myod_mef3
name <- "MyoD + MEF3 enhancers (differentiation)"
correlatediff(chips,name)

## MEF3
#print('mef3')
#datchip_mef3<-read.table("files/mef3-coords.dat")

#chips <- datchip_mef3
#name <- "MEF3 enhancers (differentiation)"
#correlatediff(chips,name)

## MYOD
#print('myod')
#datchip_myod<-read.table("files/myod_all-coords.dat")
#chips <- datchip_myod
#name <- "MyoD enhancers (differentiation)"
#correlatediff(chips,name)

dev.off()
