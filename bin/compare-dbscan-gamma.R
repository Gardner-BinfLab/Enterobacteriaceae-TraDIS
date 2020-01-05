library("dbscan")
library("MASS")
library("ROCR")
datapath <- "../results/biases/dbscan"
list_of_files <- list.files(path=datapath, full.names=T, recursive=FALSE)
names = c("ROD", "CS17", "ENC", "ETEC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261", "BW25113", "EC958")
dict = c("C. rodentium ICC168", "Escherichia coli ETEC CS17", "E. cloacae NCTC 9394", "Escherichia coli ETEC H10407", "E. coli UPEC ST131",
         "K. pneumoniae RH201207", "K. pneumoniae Ecl8", "S. Enteritidis", "S. Typhimurium A130",
         "S. Typhimurium SL1344", "S. Typhimurium D23580", "S. Typhi Ty2", "S. Typhimurium SL3261",
         "E. coli BW25113", "E. coli UPEC ST131 EC958")
names(dict) <- names
sumdb = 0
sumgamma = 0
counter = 0
pdf('../figures/gamma-vs-dbscan.pdf')
for (filename in list_of_files)
{
  ### read ecogene counterparts
  locus = sub(pattern = "(.*?)\\..*$", replacement = "\\1", basename(filename))
  real = read.table(paste('../results/ecogenecounterparts/',locus,'.txt',sep = ''), as.is=TRUE, header=FALSE, sep="\t")
  names(real) <- c('gene', 'essentiality')
  
  ### dbscan
  data <- read.table(filename)[,1:2]
  ii <- data[,2]
  
  dbres <- dbscan(as.matrix(ii), minPts = 200, eps = 0.05)
  ess <- dbres$cluster[which.min(ii)]
  dbthr <- max(ii[dbres$cluster==ess])
  
  ### gamma fit
  nG = length(ii)
  # 
  #identify second maxima
  h <- hist(ii, breaks=0:(max(ii)*50+1)/50,plot=FALSE)
  maxindex <- which.max(h$density[3:length(h$density)])
  maxval <- h$mids[maxindex+2]

  #find inter-mode minima with loess
  r <- floor(maxval *1000)
  I = ii < r / 1000
  h1 = hist(ii[I],breaks=(0:r/1000),plot=FALSE)
  lo <- loess(h1$density ~ c(1:length(h1$density))) #loess smothing over density
  m1 = h1$mids[which.min(predict(lo))]
  m2 = h$mids[max(which(h$counts>5))]
  I1 = ((ii < m1)&(ii > 0))
  I2 = ((ii >= m1)&(ii < m2))

  f1 = (sum(I1) + sum(ii == 0))/nG
  f2 = (sum(I2))/nG

  d1 = fitdistr(ii[I1], "gamma", lower=min(ii[I1]))
  d2 = fitdistr(ii[I2], "gamma", lower=min(ii[I2])) #fit curves

  #plots
  if (strsplit(basename(filename),"\\.")[[1]][1] != 'SEN')
  {
    lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
    upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, 1,d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
    essen <- lower/10000
    ambig <- upper/10000
    noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  }else
  {
    lower <- max(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) < -2))
    upper <- min(which(log((pgamma(1:20000/10000, d2$e[1],d2$e[2])*(1-pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)))/(pgamma(1:20000/10000, d1$e[1],d1$e[2], lower.tail=FALSE)*(1-pgamma(1:20000/10000, d2$e[1],d2$e[2]))) , base=2) > 2))
    essen <- lower/10000
    ambig <- upper/10000
    noness <- min(ii[pgamma(ii, d2$e[1],d2$e[2])>=0.99])
  }
  gammathr <- essen
  ### comparison
  pred <- prediction(ii, 1-real$essentiality)
  perf <- performance(pred,"tpr","fpr")
  auc <- performance(pred,measure = "auc")@y.values[[1]]
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2)
  plot(perf,col="#5ab4ac",lty=1,lwd=4,cex.lab = 2,cex.axis = 2, cex.main =2, print.cutoffs.at = c(dbthr, gammathr),
       cutoff.label.function=function(x) {c('                  DBSCAN','                Gamma')}, ylim=c(0.92,1), main=dict[locus], cex=1.5)
  perf <- performance(pred,'mat')
  par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2)
  plot(perf,col="#5ab4ac",lty=1,lwd=4,cex.lab = 2,cex.axis = 2, cex.main =2, main=dict[locus], ylim=c(0.5,0.9), xlim=c(0,0.5), xlab='Insertion index threshold')
  points(dbthr,perf@y.values[[1]][sum(perf@x.values[[1]]>dbthr)], cex=1.5)
  text(dbthr,perf@y.values[[1]][sum(perf@x.values[[1]]>dbthr)], '                  DBSCAN',cex = 1.5)
  points(gammathr,perf@y.values[[1]][sum(perf@x.values[[1]]>gammathr)], cex=1.5)
  text(gammathr,perf@y.values[[1]][sum(perf@x.values[[1]]>gammathr)], '                Gamma', cex=1.5)
  if (dict[locus]=='E. coli BW25113')
  {
    par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2)
    plot(perf,col="#5ab4ac",lty=1,lwd=4,cex.lab = 2,cex.axis = 2, cex.main =2, main=dict[locus], ylim=c(0.5,0.75), xlim=c(0,0.4), xlab='Insertion index threshold')
    points(dbthr,perf@y.values[[1]][sum(perf@x.values[[1]]>dbthr)], cex=1.5)
    text(dbthr,perf@y.values[[1]][sum(perf@x.values[[1]]>dbthr)], '                  DBSCAN',cex = 1.5)
    points(gammathr,perf@y.values[[1]][sum(perf@x.values[[1]]>gammathr)], cex=1.5)
    text(gammathr,perf@y.values[[1]][sum(perf@x.values[[1]]>gammathr)], '                Gamma', cex=1.5)
  }
  sumdb = sumdb + perf@y.values[[1]][sum(perf@x.values[[1]]>dbthr)]
  sumgamma = sumgamma + perf@y.values[[1]][sum(perf@x.values[[1]]>gammathr)]
  counter = counter + 1
}
dev.off()
# print average MCC values for DBSCAN and GAMMA
print(sumdb/counter)
print(sumgamma/counter)