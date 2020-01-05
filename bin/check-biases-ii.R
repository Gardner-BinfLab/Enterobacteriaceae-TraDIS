library(stringr)
library("MASS")
library(gdata)
library(mgcv)
library("dbscan")
library("ROCR")
library("ggplot2")
cols <- c("#8c510a", "#80cdc1", "#01665e")
# library("plotrix")
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261", "BW25113", "EC958")
dict = c("C. rodentium ICC168", "Escherichia coli ETEC CS17", "E. cloacae NCTC 9394", "Escherichia coli ETEC H10407", "E. coli UPEC ST131",
         "K. pneumoniae RH201207", "K. pneumoniae Ecl8", "S. Enteritidis", "S. Typhimurium A130",
         "S. Typhimurium SL1344", "S. Typhimurium D23580", "S. Typhi Ty2", "S. Typhimurium SL3261",
         "E. coli BW25113", "E. coli UPEC ST131 EC958")
names(dict) <- names
genome_length = c(5346659, 4994793, 4908759, 5153435, 5174631, 5869288, 5324709, 4685848, 4895639, 4878012, 4879400, 4791961, 4878012, 4631469, 5109767)
names(genome_length) <- names
sp = 0.2
biasespath <- "../results/biases/check-biases/"
outdir <- "../results/biases/dbscan"
dir.create(outdir)
list_of_files <- list.files(path=biasespath, full.names=T, recursive=FALSE)
iitotal = c()
gctotal = c()
dtotal = c()
postotal = c()
ii_dnormalisedtotal = c()
ii_dgcnormalisedtotal = c()
pdf("../figures/biases.pdf")
for (filename in list_of_files)
{
  biasestable = read.table(filename, header = FALSE, stringsAsFactors = FALSE)
  colnames(biasestable) <- c("name", "ii", "essentiality", "dist", "gc", "length", "pos")
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  ii <- c(biasestable$ii)
  iitotal <- c(iitotal, ii)
  gc <- c(biasestable$gc)
  gctotal <- c(gctotal , gc)
  d <- c(biasestable$dist)
  dtotal <- c(dtotal , d)
  pos <- c(biasestable$pos)
  
  pos <- sort(pos)
  sortedd <- d[order(pos)]
  mind = which.min(sortedd)
  maxd = which.max(sortedd)
  lfdist = genome_length[strsplit(basename(filename), "\\.")[[1]][1]][[1]] - pos[mind]
  pos[mind:length(pos)] = pos[mind:length(pos)] - pos[mind] + 1
  pos[1:mind-1] = lfdist + pos[1:mind-1]
  postotal <- c(postotal, pos)
  # plot(pos, ii, pch = '.', ylim=c(0,2), xlab = "Gene position", ylab = "insertion index",
  #      main = paste("Distance bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2, xaxt='n')
  # axis(1,at=c(1, pos[maxd], pos[mind-1]),labels = c('Origin', 'Terminus', 'Origin'), cex.axis=2)
  # # abline(mean(ii),0,col='green', lwd=5)
  # lines(loess.smooth(pos,ii, span=sp), col=2, lwd=5)
  dat <- data.frame(pos,ii)
  print(ggplot(dat,aes(pos,ii))+geom_point(color=cols[1], shape='.') +geom_smooth(method = "gam", formula = y ~ s(x), fill=cols[2], color=cols[3]) + ylim(0,2)+
          xlab("Gene position") +ylab("Insertion index") + ggtitle(paste("Distance bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]))+
          theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=20, colour = 'black')) + theme(axis.title.x = element_text(size=24)) +
          theme(axis.title.y = element_text(size=24))+ guides(fill=FALSE) +
          scale_x_continuous(breaks=c(1, pos[maxd], pos[mind-1]), labels=c('Origin', 'Terminus', 'Origin')) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          panel.border = element_rect(colour = "black", fill=NA, size=0.5)))
  
  plot(d, ii, pch = '.', ylim=c(0,2), xlab = "Gene position", ylab = "Insertion index",
       main = paste("Distance bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2)
  # abline(mean(ii),0,col='green', lwd=5)
  lines(loess.smooth(d,ii, span=sp), col=2, lwd=5)
  
  # plot(gc, ii, pch = '.', ylim=c(0,5), xlab = "GC content", ylab = "insertion index",
  #      main = paste("GC bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2)
  # abline(mean(ii),0,col='green', lwd=5)
  # lines(loess.smooth(gc,ii, span=sp), col=2, lwd=5)
  dat <- data.frame(gc,ii)
  print(ggplot(dat,aes(gc,ii))+geom_point(color=cols[1], shape='.') +geom_smooth(method = "gam", formula = y ~ s(x), fill=cols[2], color=cols[3]) + ylim(0,2)+
    xlab("GC content") +ylab("Insertion index") + ggtitle(paste("GC bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]))+ theme_bw()+
      theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=20, colour = 'black')) + theme(axis.title.x = element_text(size=24)) +
      theme(axis.title.y = element_text(size=24))+ guides(fill=FALSE) + 
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")))
  
  newgc = gc[biasestable$essentiality != 'essential' & ii>0.1]
  newii = ii[biasestable$essentiality != 'essential'&ii>0.1]
  # plot(newgc, newii, pch = '.', ylim=c(0,5), xlab = "GC content", ylab = "insertion index",
  #      main = paste("GC bias without essential genes -", dict[strsplit(basename(filename), "\\.")[[1]][1]]), cex.lab = 2, cex.axis = 2, cex.main =2)
  # lines(loess.smooth(newgc,newii, span=sp), col=2, lwd=5)
  dat <- data.frame(newgc,newii)
  print(ggplot(dat,aes(newgc,newii))+geom_point(color=cols[1], shape='.') +geom_smooth(method = "gam", formula = y ~ s(x), fill=cols[2], color=cols[3]) + ylim(0,2)+
          xlab("GC content") +ylab("Insertion index") + ggtitle(paste("GC bias -", dict[strsplit(basename(filename), "\\.")[[1]][1]]))+ theme_bw()+
          theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=20, colour = 'black')) + theme(axis.title.x = element_text(size=22)) +
          theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) + 
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black")))
  
  fit <- loess(ii~d, span = sp)
  loessprediction <- predict(fit, d)
  ii_dnormalised = ii/(loessprediction/mean(ii))
  ii_dnormalisedtotal <- c(ii_dnormalisedtotal, ii_dnormalised)
  plot(d, ii_dnormalised, pch='.', xlab = "Gene position", ylab = "Insertion index",
       main = paste(dict[strsplit(basename(filename), "\\.")[[1]][1]], '- normalised distance'), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,5))
  lines(loess.smooth(d,ii_dnormalised, span=sp), col=2, lwd=5)
  
  fit <- loess(ii_dnormalised~gc, span =sp)
  loessprediction <- predict(fit, gc)
  ii_dgcnormalised = ii_dnormalised/(loessprediction/mean(ii_dnormalised))
  ii_dgcnormalisedtotal <- c(ii_dgcnormalisedtotal, ii_dnormalised)
  plot(gc, ii_dgcnormalised, pch='.', xlab = "GC content", ylab = "Insertion index",
       main = paste(dict[strsplit(basename(filename), "\\.")[[1]][1]], '- normalised GC'), cex.lab = 2, cex.axis = 2, cex.main =2,
       ylim=c(0,5))
  lines(loess.smooth(gc,ii_dgcnormalised, span=sp), col=2, lwd=5)
}

plot(gctotal, iitotal, pch = '.', ylim=c(0,5), xlab = "GC content", ylab = "Insertion index", main = "GC bias", 
     cex.lab = 2, cex.axis = 2, cex.main =2)
lines(loess.smooth(gctotal,iitotal, span=sp), col=2, lwd=5)

dev.off()

pdf("~/EnTrI/figures/insertion-indices-normalised.pdf")
s = 1
for (filename in list_of_files)
{
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  biasestable = read.table(filename, header = FALSE, stringsAsFactors = FALSE)
  colnames(biasestable) <- c("name", "ii", "essentiality", "dist", "gc")
  biasestable$ii <- ii_dgcnormalisedtotal[s:(s+length(biasestable$ii)-1)]
  # for (i in seq(1,length(biasestable$ii)))
  # {
  #   if (biasestable$ii[i] < 0)
  #   {
  #     print(biasestable$ii[i])
  #     biasestable$ii[i] = 0
  #   }
  # }
  ii=biasestable$ii
  # nG = length(ii)
  # 
  # #identify second maxima
  # h <- hist(ii, breaks=0:(max(ii)*50+1)/50,plot=FALSE)
  # maxindex <- which.max(h$density[3:length(h$density)])
  # maxval <- h$mids[maxindex+2]
  # 
  # #find inter-mode minima with loess
  # r <- floor(maxval *1000)
  # I = ii < r / 1000
  # h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
  # lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  # m1 = h1$mids[which.min(predict(lo))]
  # m2 = h$mids[max(which(h$counts>5))]
  # I1 = ((ii < m1)&(ii > 0))
  # I2 = ((ii >= m1)&(ii < m2))
  # 
  # f1 = (sum(I1) + sum(ii == 0))/nG
  # f2 = (sum(I2))/nG
  # 
  # # if (strsplit(basename(filename),"\\.")[[1]][1] == 'SEN')
  # # {
  # #   m3 = h$mids[min(which(diff(sign(diff(h$counts)))==-2)+1)] # the index of the first local maximum
  # #   I1 = ((ii < m1)&(ii >= m3))
  # # }
  # # 
  # d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
  # d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves
  # 
  # #plots
  # if (strsplit(basename(filename),"\\.")[[1]][1] != 'SEN')
  # {
  #   hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  #   lines(0:2000/500, f1*dgamma(0:2000/500, 1, d1$estimate[2])) # was [2]
  #   lines(0:2000/500, f2*dgamma(0:2000/500, d2$estimate[1], d2$estimate[2]))
  #   lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
  #   upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
  #   essen <- lower/10000
  #   ambig <- upper/10000
  #   noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  # }
  # else
  # {
  #   hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  #   lines(0:2000/500, f1*dgamma(0:2000/500, d1$estimate[1], d1$estimate[2])) # was [2]
  #   lines(0:2000/500, f2*dgamma(0:2000/500, d2$estimate[1], d2$estimate[2]))
  #   lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
  #   upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
  #   essen <- lower/10000
  #   ambig <- upper/10000
  #   noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  # }
  # 
  # lines(c(essen, essen), c(0,20), col="red")
  # lines(c(ambig, ambig), c(0,20), col="red")
  # lines(c(noness, noness), c(0,20), col="red")
  # 
  # mtext(paste(essen, ":", "Essential changepoint"), side=3, adj=1, padj=2)
  # mtext(paste(ambig, ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
  # mtext(paste(noness, ":", "Non-essential changepoint"), side=3, adj=1, padj=5.5)
  # 
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  res <- dbscan(as.matrix(ii), minPts = 200, eps = 0.05)
  # res <- optics(as.matrix(ii), minPts = 200, eps = 0.05)
  # res <- extractDBSCAN(res, 0.05)
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(ii,col=res$cluster+1, xlab='Genes', ylab='Insertion index', cex.axis = 1.5, cex.main = 2, cex.lab = 2, main=dict[locusid])
  ess <- res$cluster[which.min(ii)]
  print(sum(res$cluster == ess))
  nes <- getmode(res$cluster)
  belthr <- max(ii[res$cluster==nes])
  essthr <- max(ii[res$cluster==ess])
  nesthr <- min(ii[res$cluster==nes])
  
  bel <- ii[ii>belthr]
  res2 <- dbscan(as.matrix(bel), minPts = 100, eps = 0.1)
  ambig <- res2$cluster[which.min(bel)]
  belam <- max(bel[res2$cluster==ambig])
  
  for (i in (1:length(ii)))
  {
    if (res$cluster[i] == ess)
    {
      biasestable$essentiality[i] = "essential"
    } else if(res$cluster[i] == nes)
    {
      biasestable$essentiality[i] = "non-essential"
    } else if(ii[i] > belthr)
    {
      if (ii[i] < belam)
      {
        biasestable$essentiality[i] = "benloss-ambiguous"
      }
      else
      {
        biasestable$essentiality[i] = "beneficial-loss"
      }
    } else
    {
      biasestable$essentiality[i] = "essential-ambiguous"
    }
  }
  
  h <- hist(ii, breaks =seq(min(ii),max(ii)+1,0.02), plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,essthr, nesthr, belthr, belam, Inf))
  par(mar = mar.default + c(0, 1, 0, 0))
  max1 = max(h$counts)
  max2 = sort(h$counts,partial=length(h$counts)-1)[length(h$counts)-1]
  plot(h, col=c(cols[1], cols[2], cols[3], cols[3], cols[3])[cuts], xlab = "Insertion index", main =dict[locusid], cex.lab = 2,
       cex.axis = 2, cex.main = 2, xlim=c(0,4), ylim=c(0,max1), lty= "blank")
  box()
  text(2.5,max1-20, paste("n =", length(ii)), lty=1, lwd=4, cex=1.5, bty="n")
  # legend(2,max1-50, c("Essential", "Ambiguous", "Non-essential", "Beneficial loss"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.5,
  #        col=c("sienna4", "gray", "midnightblue", "red"), bty="n")
  legend(2,max1-50, c("Essential", "Ambiguous", "Non-essential"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.5,
        col=c(cols[1], cols[2], cols[3]), bty="n")
  
  if (dict[locusid]=='E. coli BW25113')
  {
    max1 = 100
    plot(h, col=c(cols[1], cols[2], cols[3], cols[3], cols[3])[cuts], xlab = "Insertion index", main =dict[locusid], cex.lab = 2,
         cex.axis = 2, cex.main = 2, xlim=c(0,4), ylim=c(0,max1), lty= "blank")
    box()
    text(2.5,max1-20, paste("n =", length(ii)), lty=1, lwd=4, cex=1.5, bty="n")
    legend(2,max1-30, c("Essential", "Ambiguous", "Non-essential"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.5,
           col=c(cols[1], cols[2], cols[3]), bty="n")
  }
  
  # hist(ii,breaks=0:(max(ii)*50+1)/50, xlim=c(0,4), freq=FALSE,xlab="Insertion index", main=dict[locusid])
  # lines(c(essthr, essthr), c(0,20), col="red")
  # lines(c(nesthr, nesthr), c(0,20), col="red")
  # lines(c(belthr, belthr), c(0,20), col="red")
  # mtext(paste(round(essthr,digits = 3), ":", "Essential changepoint"), side=3, adj=1, padj=2)
  # mtext(paste(round(nesthr,digits = 3), ":", "Ambiguous changepoint"), side=3, adj=1, padj=3.75)
  # mtext(paste(round(belthr,digits = 3), ":", "Non-essential changepoint"), side=3, adj=1, padj=5.5)
  
  esslevel <- ifelse(ii > 0, log(ii/essthr), -4.5)
  
  write.table(cbind(biasestable[1:3], esslevel), file = paste(outdir, '/', basename(filename), sep=''), quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
  s = s + length(biasestable$ii)
}

dev.off()