#!/usr/bin/Rscript
pdf("plotVaryFC.pdf", family = "Palatino")
cex.main <- 2
cex.lab <- 2
cex.axis <- 2
lwd <- 5
fcmax <- 4
step <- 0.1

# fibroblasts
fibroswt<-read.table('/Users/Marc/Desktop/work/these/projects/cochin/pascal/fibros-six+myod-synergy/affy+qPCR/matrice_filtered.txt',header=TRUE)
matfibrowtMT<-as.matrix(fibroswt[,c(5,10,6,11)])
matfibrowtMB<-as.matrix(fibroswt[,c(2,7,3,8)])
fibros<-cbind(matfibrowtMB,matfibrowtMT)
rownames(fibros)<-fibroswt[,1]


# MYOD + MEF3
allratiosdiff <- list()
i <- 1

for (file in c('files/myod+mef3-coords.dat','files/mef3-coords.dat','files/myod_all-coords.dat')){
    print(file)
    datchip_myod_mef3<-read.table(paste('/Users/Marc/Desktop/work/these/projects/cochin/pascal/fibros-six+myod-synergy/affy+qPCR/affy/',file,sep=""))
    chip_myod_mef3<-unique(sort(as.character(datchip_myod_mef3$V1)))

    fibroswchip <- fibros[which(rownames(fibros) %in% chip_myod_mef3),]
    FCsprolifwchip<-2^(fibroswchip[,4]-apply(fibroswchip[,1:3],1,max))
    FCsdiffwchip<-2^(fibroswchip[,8]-apply(fibroswchip[,5:7],1,max))
    FCsprolif<-2^(fibros[,4]-apply(fibros[,1:3],1,max))
    FCsdiff<-2^(fibros[,8]-apply(fibros[,5:7],1,max))

    #FCsprolifwchip[FCsprolifwchip < 1] <- - 1 / FCsprolifwchip[FCsprolifwchip < 1]
    #FCsprolif[FCsprolif < 1] <- - 1 / FCsprolif[FCsprolif < 1]
    #FCsdiffwchip[FCsdiffwchip < 1] <- - 1 / FCsdiffwchip[FCsdiffwchip < 1]
    #FCsdiff[FCsdiff < 1] <- - 1 / FCsdiff[FCsdiff < 1]
    
    #FCsprolifwchip<-(fibroswchip[,4]-apply(fibroswchip[,1:3],1,max))
    #FCsdiffwchip<-(fibroswchip[,8]-apply(fibroswchip[,5:7],1,max))
    #FCsprolif<-(fibros[,4]-apply(fibros[,1:3],1,max))
    #FCsdiff<-(fibros[,8]-apply(fibros[,5:7],1,max))


    #ratiosprolif <- NULL
    #ratiosdiff <- NULL
    #fcneg <- seq(0,fcmax,step)
    #for (fc in fcneg){
        #indchip <- which(FCsprolifwchip <= fc)
        #ind <- which(FCsprolif <= fc)
        #ratiosprolif <- c(ratiosprolif, length(indchip)/length(ind))

        #indchip <- which(FCsdiffwchip <= fc)
        #ind <- which(FCsdiff <= fc)
        #ratiosdiff <- c(ratiosdiff, length(indchip)/length(ind))
    #}
    #ymax <- max(na.omit(c(ratiosprolif,ratiosdiff)))
    #plot(fcneg, ratiosprolif, type = "l", col = "blue", main = file, xlab = "FC", ylab = "Proportion with enhancer", lwd = 4, ylim = c(0, ymax))
    #lines(fcneg, ratiosdiff, type = "l", col = "red", lwd = 4)
    #legend("topleft", leg = c("Prolif","Diff"), lty = c(1,1), col = c("blue","red"), lwd = 4)
    
    ratiosprolif <- NULL
    ratiosdiff <- NULL
    fcpos <- seq(0,fcmax,step)
    for (fc in fcpos){
        indchip <- which(FCsprolifwchip >= fc)
        ind <- which(FCsprolif >= fc)
        ratiosprolif <- c(ratiosprolif, length(indchip)/length(ind))

        indchip <- which(FCsdiffwchip >= fc)
        ind <- which(FCsdiff >= fc)
        ratiosdiff <- c(ratiosdiff, length(indchip)/length(ind))
    }

    #ymax <- max(na.omit(c(ratiosprolif,ratiosdiff)))
    #plot(fcpos, ratiosprolif, type = "l", col = "blue", main = file, xlab = "FC", ylab = "Proportion with enhancer", lwd = 4, ylim = c(0, ymax))
    #lines(fcpos, ratiosdiff, type = "l", col = "red", lwd = 4)
    #legend("topleft", leg = c("Prolif","Diff"), lty = c(1,1), col = c("blue","red"), lwd = 4)
    
    allratiosdiff[[i]] <- ratiosdiff
    i <- i + 1

}

for (i in 1:3)
    print(allratiosdiff[[i]][1])


par(mar=c(6,6,6,6))
plot(fcpos, 
     allratiosdiff[[1]] / allratiosdiff[[1]][1], 
     col = "purple", 
     lwd = lwd, 
     type = "l", 
     main = "", 
     xlab = "Fold-Change en affymetrix", 
     ylab = "Facteur d'enrichissement", 
     cex.main = cex.main,
     cex.lab = cex.lab,
     cex.axis = cex.axis,
     ylim = c(1, 6)
     )
lines(fcpos, 
     allratiosdiff[[2]] / allratiosdiff[[2]][1], 
     lwd = lwd, 
     col = "red"
     )
lines(fcpos, 
     allratiosdiff[[3]] / allratiosdiff[[3]][1], 
     lwd = lwd, 
     col = "blue"
     )
legend("topleft",
       leg = c("MyoD + MEF3", "MEF3 cons","MyoD ChIP-seq" ),
       cex = 1.5,
       lty = 1,
       lwd = lwd,
       col = c("purple","red","blue")
       )
