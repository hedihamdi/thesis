#!/usr/bin/Rscript

library(gplots)
by <- 0.1

dir.create('FIGS-STATISTICS', showWarnings = FALSE)

xlim <- c(-5,1)
ylim <- c(0,10)
breaks <- seq(from = -5, to = 2, by = 0.05)

plotfig <- function(xlab, col, filepdf, filetxt){ 

    print(filetxt)
    pdf(filepdf, family = "Palatino")
    par(mar=c(6,6,4,4))
    t <- read.table(filetxt)
    h <- hist(log10(abs(t$V1)), 
              breaks = breaks, 
              plot = FALSE,
              )

    plot(
         h$mids,
         h$count / 22, # the 22 TFs out of 28 that have better pairwise model than PWM 
         #h$density * by,
         #h$count / sum(h$count),
         #log="x", 
         type='h', 
         lwd=4, 
         lend=2,
         xlab = xlab, 
         ylab = "Nombre par TF",
         main = "",
         cex.axis = 2,
         #cex.main = 3,
         cex.lab = 2.5,
         xlim = xlim,
         ylim = ylim,
         xaxt='n',
         col = col
         )
    ticks <- seq(-4,2,1)
    labs <- c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),'1',expression(10^1),expression(10^2))
    axis(1, at = ticks, labels=labs, cex.axis = 2)
}

plotfig(expression(abs(h[inde])),
        "blue",
        'FIGS-STATISTICS/hinde.pdf',
        "RES-STATISTICS/hinde.txt")
plotfig(expression(abs(h[pw])),
        "red",
        'FIGS-STATISTICS/hpw.pdf',
        "RES-STATISTICS/hpw.txt")
plotfig(expression(abs(J)),
        "red",
        'FIGS-STATISTICS/J.pdf',
        "RES-STATISTICS/J.txt")


# Now for interactions
print('Interactions')

pdf('FIGS-STATISTICS/Jnn.pdf', family = "Palatino")
par(mar=c(6,6,4,4))

ylim <- c(0,10)
t <- read.table('RES-STATISTICS/J.txt')
data <- c(-log10(-t$V1[t$V1<0]),log10(t$V1[t$V1>0]))
data <- t$V1[t$V1!=0]
data <- log10(abs(t$V1))
#breaks <- seq(from = min(data)-1, to = max(data)+1, by = by)
#h <- hist(log10(abs(t$V1)), 
h <- hist(data, 
          breaks = breaks, 
          plot = FALSE,
          )

allcounts <- h$count

plot(
     h$mids,
     h$count / 22, # the 22 TFs out of 28 that have better pairwise model than PWM 
     #h$density * by,
     #h$count / sum(h$count),
     #log="x", 
     type='h', 
     lwd=4, 
     lend=2,
     xlab = expression(J), 
     ylab = "Nombre par TF",
     main = "",
     cex.axis = 2,
     #cex.main = 3,
     cex.lab = 2.5,
     xlim = xlim,
     ylim = ylim,
     xaxt='n',
     col = "red"
     )
ticks <- seq(-4,2,1)
labs <- c(expression(10^-4),expression(10^-3),expression(10^-2),expression(10^-1),'1',expression(10^1),expression(10^2))
axis(1, at = ticks, labels=labs, cex.axis = 2)

t <- read.table('RES-STATISTICS/J_nearest_neighbour.txt')
data <- c(-log10(-t$V1[t$V1<0]),log10(t$V1[t$V1>0]))
data <- t$V1[t$V1!=0]
data <- log10(abs(t$V1))
#h <- hist(log10(abs(t$V1)), 
h <- hist(data, 
          breaks = breaks, 
          plot = FALSE,
          )
hprop <- h$count / sum(allcounts)

points(
       h$mids,
       h$count / 22, # the 22 TFs out of 28 that have better pairwise model than PWM 
       #hprop,
       type = 'h',
       lwd=4,
       col = 'black'
       )

legend('topleft', c('Toutes les interactions', 'Plus proches voisins'), fill=c('red','black'),cex=1.5)
