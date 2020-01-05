library(stringr)
library("dbscan")
library("ggplot2")
fasta_dir <- "~/EnTrI/data/fasta-protein/chromosome/Salmonella_enterica_subsp_enterica_serovar_Typhimurium_SL1344_FQ312003_v4.fasta"
plots_dir <- "~/EnTrI/data/plot-files/density/"
output_dir1 <- "~/EnTrI/results/density/"
ecogene <- "~/EnTrI/results/ecogenecounterparts/SL1344.txt"

getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
mar.default <- c(5,4,4,2) + 0.1

ecoessentiality <- read.table(ecogene, row.names = 1, col.names = c('name','essentiality'))
list_of_files <- list.files(path=plots_dir, full.names=T, recursive=FALSE)
plots = list()
num_inserts = list()
den = list()
for (filename in list_of_files)
{
  plotfile = as.matrix(read.table(filename, as.is=TRUE))
  locusid = strsplit(basename(filename),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  num_inserts[[locusid]] = sum(plots[[locusid]]>0)
  len = length(plots[[locusid]])
  den[[locusid]] = len/num_inserts[[locusid]]
}

iitable = data.frame(locus.tag=NA, gene.length=NA ,SL1344_1=NA, SL1344_2=NA, SL1344_3=NA, SL1344_4=NA, SL1344_5=NA, SL1344_6=NA, SL1344_7=NA,
                     SL1344_8=NA, SL1344_9=NA, SL1344_10=NA)
noninsfree = list(SL1344_1=0, SL1344_2=0, SL1344_3=0, SL1344_4=0, SL1344_5=0, SL1344_6=0, SL1344_7=0,
               SL1344_8=0, SL1344_9=0, SL1344_10=0)
geneinfo <- matrix(ncol=3)
colnames(geneinfo) <- c('locustag', 'start', 'end')
fastafile = readLines(fasta_dir)
for (line in fastafile)
{
  if (startsWith(line, ">"))
  {
    matchresult = str_match(line, ">([[:graph:]]+)[[:blank:]]\\[[[:graph:]]+\\/([[:digit:]]+)\\-([[:digit:]]+)[[:print:]]+\\(([[:alpha:]]+)\\)")
    if (!(is.na(matchresult[5])))
    {
      locustag = matchresult[2]
      start = as.numeric(matchresult[3])
      end = as.numeric(matchresult[4])
      direction = matchresult[5]
      leng = end - start + 1
      starttrim = floor(leng*5/100)
      endtrim = floor(leng*20/100)
      genelength = "short"
      if (direction == "Forward")
      {
        newstart = start + starttrim
        newend = end - endtrim
        if (100 < leng)
        {
          genelength = "long"
        }
      }
      else
      {
        newstart = start + endtrim
        newend = end - starttrim
        if (100 < leng)
        {
          genelength = "long"
        }
      }
      newlen = newend - newstart + 1
      iirow <- c(locustag, genelength)
      for (locusid in names(iitable)[c(-1,-2)])
      {
        ii = (sum(plots[[locusid]][newstart:newend]>0) / newlen) / (num_inserts[[locusid]] / len)
        iirow <- c(iirow, ii)
        if (sum(plots[[locusid]][newstart:newend]>0))
        {
          noninsfree[[locusid]] = noninsfree[[locusid]] + 1
        }
      }
      iitable <- rbind(iitable, iirow)
      geneinfo <- rbind(geneinfo, c(locustag, newstart, newend))
    }
  }
}
iitable <- iitable[-1,]
geneinfo <- geneinfo[-1,]
# insfree <- noninsfree
# for (locustag in names(noninsfree))
# {
#   insfree[[locustag]] = length(iitable$locus.tag) - noninsfree[[locustag]]
# }

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_1)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_1)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr1 = fp/(fp+tn)
tpr1 = tp/(tp+fn)
mcc1 = (tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
fps1 = fp
essens1 = sum(res$cluster==2)
unif <- rep(0, len)
ind <- sample(len, num_inserts$SL1344_1, replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/num_inserts$SL1344_1/len
}
essthr <- max(as.numeric(iitable$SL1344_1)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif1 = fpunif/(fpunif+tnunif)
essensunif1 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_3)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_3)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr3 = fp/(fp+tn)
tpr3 = tp/(tp+fn)
mcc3 = (tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
fps3 = fp
essens3 = sum(res$cluster==2)
unif <- rep(0, len)
ind <- sample(len, num_inserts$SL1344_3, replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/num_inserts$SL1344_3/len
}
essthr <- max(as.numeric(iitable$SL1344_3)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif3 = fpunif/(fpunif+tnunif)
essensunif3 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_7)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_7)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr7 = fp/(fp+tn)
tpr7 = tp/(tp+fn)
mcc7 = (tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
fps7 = fp
essens7 = sum(res$cluster==2)
unif <- rep(0, len)
ind <- sample(len, num_inserts$SL1344_7, replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/num_inserts$SL1344_7/len
}
essthr <- max(as.numeric(iitable$SL1344_7)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif7 = fpunif/(fpunif+tnunif)
essensunif7 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_5)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_5)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr5 = fp/(fp+tn)
tpr5 = tp/(tp+fn)
mcc5 = (tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
fps5 = fp
essens5 = sum(res$cluster==2)
unif <- rep(0, len)
ind <- sample(len, num_inserts$SL1344_5, replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/num_inserts$SL1344_5/len
}
essthr <- max(as.numeric(iitable$SL1344_5)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif5 = fpunif/(fpunif+tnunif)
essensunif5 = sum(iisunif <= essthr)

res <- dbscan(as.matrix(as.numeric(iitable$SL1344_9)), minPts = 200, eps = 0.05)
# plot(as.matrix(as.numeric(iitable$SL1344_9)),col=res$cluster+1, pch=20)
tp = sum(res$cluster==2 & ecoessentiality=="1")
fp = sum(res$cluster==2 & ecoessentiality!="1")
fn = sum(res$cluster!=2 & ecoessentiality=="1")
tn = sum(res$cluster!=2 & ecoessentiality!="1")
fpr9 = fp/(fp+tn)
tpr9 = tp/(tp+fn)
mcc9 = (tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))
fps9 = fp
essens9 = sum(res$cluster==2)
unif <- rep(0, len)
ind <- sample(len, num_inserts$SL1344_9, replace = FALSE)
unif[ind] = 1
iisunif <- rep(0, nrow(geneinfo))
for (j in seq(1, length(iisunif)))
{
  iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
  iisunif[j] <- iisunif[j]/num_inserts$SL1344_9/len
}
essthr <- max(as.numeric(iitable$SL1344_9)[res$cluster==2])
fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
tnunif = sum(iisunif > essthr & ecoessentiality!="1")
fprunif9 = fpunif/(fpunif+tnunif)
essensunif9 = sum(iisunif <= essthr)

one_rep_names <- c("SL1344_1", "SL1344_3", "SL1344_5", "SL1344_7", "SL1344_9")
one_rep_den <- c(den$SL1344_1,den$SL1344_3[1],den$SL1344_5[1],den$SL1344_7[1],den$SL1344_9[1])
sorted_one_rep_den <- sort(one_rep_den)
sorted_one_rep_names <- one_rep_names[order(one_rep_den)]

fprs = c()
fprsunif = c()
essens = c()
essensunif = c()
dens = c()
# ins = c()
ins = seq(2e4,55e4,5e3)
tprs = c()
mccs = c()
fps = c()
noninsfrees = c()
# for (i in seq(9, 200, 1))
pdf('~/EnTrI/figures/false-positive-rate_density-histograms.pdf')
for (num in ins)
{
  # num = round(len/i)
  i = len/num
  dens <- c(dens, len/num)
  # ins <- c(ins, num)
  samplefrom <- min(sorted_one_rep_names[sorted_one_rep_den<i])
  sampledplot <- plots[[samplefrom]]
  iis <- rep(0, nrow(geneinfo))
  ind <- sample(which(sampledplot > 0), sum(sampledplot>0)-num, replace = FALSE)
  sampledplot[ind] = 0
  unif <- rep(0, length(sampledplot))
  ind <- sample(length(sampledplot), num, replace = FALSE)
  unif[ind] = 1
  iisunif <- rep(0, nrow(geneinfo))
  counter = 0 # counts the number of genes with at least one insertion
  for (j in seq(1, length(iis)))
  {
    iis[j] <- sum(sampledplot[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
    iis[j] <- iis[j]/(num/length(sampledplot))
    iisunif[j] <- sum(unif[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0)/(as.numeric(geneinfo[j,3])-as.numeric(geneinfo[j,2])+1)
    iisunif[j] <- iisunif[j]/(num/length(sampledplot))
    if (sum(sampledplot[as.numeric(geneinfo[j,2]):as.numeric(geneinfo[j,3])]>0))
    {
      counter = counter + 1
    }
  }
  noninsfrees <- c(noninsfrees, counter)
  res <- dbscan(as.matrix(iis), minPts = 200, eps = 0.05)
  ess <- res$cluster[which.min(iis)]
  essthr <- max(iis[res$cluster==ess])
  nes <- getmode(res$cluster)
  belthr <- max(iis[res$cluster==nes])
  nesthr <- min(iis[res$cluster==nes])
  fp = sum(res$cluster==ess & ecoessentiality!="1")
  fpunif = sum(iisunif <= essthr & ecoessentiality!="1")
  tn = sum(res$cluster!=ess & ecoessentiality!="1")
  tnunif = sum(iisunif > essthr & ecoessentiality!="1")
  tp = sum(res$cluster==ess & ecoessentiality=="1")
  fn = sum(res$cluster!=ess & ecoessentiality=="1")
  fpr = fp/(fp+tn)
  fprunif = fpunif/(fpunif+tnunif)
  fprs <- c(fprs, fpr)
  fprsunif <- c(fprsunif, fprunif)
  essens <- c(essens, sum(res$cluster==ess))
  essensunif <- c(essensunif, sum(iisunif <= essthr))
  tprs = c(tprs,(tp/(tp+fn)))
  mccs = c(mccs,((tp*tn-fp*fn)/(sqrt(tp+fp)*sqrt(tp+fn)*sqrt(tn+fp)*sqrt(tn+fn))))
  fps = c(fps,fp)
  if (num %in% seq(2e4,52e4,1e5))
  {
    h <- hist(iis, breaks =seq(min(iis),max(iis)+1,0.02), plot=FALSE)
    cuts <- cut(h$breaks, c(-Inf,essthr, nesthr, belthr, Inf))
    par(mar = mar.default + c(0, 1, 0, 0))
    max1 = max(h$counts)
    max2 = sort(h$counts,partial=length(h$counts)-1)[length(h$counts)-1]
    plot(h, col=c("darkgoldenrod4", "black", "midnightblue", "red")[cuts], xlab = "Insertion index",
         main =paste("density:", as.character(format(num/len, digits = 3))),
         cex.lab = 2, cex.axis = 2, cex.main = 2, xlim=c(0,4), ylim=c(0,max1), lty= "blank")
    # text(2.5,max1-20, paste("n =", length(iis)), lty=1, lwd=4, cex=1.5, bty="n")
    legend(2,max1-50, c("Essential", "Ambiguous", "Non-essential", "Beneficial loss"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.5,
           col=c("darkgoldenrod4", "black", "midnightblue", "red"), bty="n")
  }
}
dev.off()

dens <- append(dens, den$SL1344_9, after=107)
dens <- append(dens, den$SL1344_5, after=65)
dens <- append(dens, den$SL1344_7, after=56)
dens <- append(dens, den$SL1344_3, after=15)
dens <- append(dens, den$SL1344_1, after=14)

ins <- append(ins, num_inserts$SL1344_9, after=107)
ins <- append(ins, num_inserts$SL1344_5, after=65)
ins <- append(ins, num_inserts$SL1344_7, after=56)
ins <- append(ins, num_inserts$SL1344_3, after=15)
ins <- append(ins, num_inserts$SL1344_1, after=14)

fprs <- append(fprs, fpr9, after=107)
fprs <- append(fprs, fpr5, after=65)
fprs <- append(fprs, fpr7, after=56)
fprs <- append(fprs, fpr3, after=15)
fprs <- append(fprs, fpr1, after=14)

essens <- append(essens, essens9, after=107)
essens <- append(essens, essens5, after=65)
essens <- append(essens, essens7, after=56)
essens <- append(essens, essens3, after=15)
essens <- append(essens, essens1, after=14)

tprs <- append(tprs, tpr9, after=107)
tprs <- append(tprs, tpr5, after=65)
tprs <- append(tprs, tpr7, after=56)
tprs <- append(tprs, tpr3, after=15)
tprs <- append(tprs, tpr1, after=14)

mccs <- append(mccs, mcc9, after=107)
mccs <- append(mccs, mcc5, after=65)
mccs <- append(mccs, mcc7, after=56)
mccs <- append(mccs, mcc3, after=15)
mccs <- append(mccs, mcc1, after=14)

fps <- append(fps, fps9, after=107)
fps <- append(fps, fps5, after=65)
fps <- append(fps, fps7, after=56)
fps <- append(fps, fps3, after=15)
fps <- append(fps, fps1, after=14)

noninsfrees <- append(noninsfrees, noninsfree$SL1344_9, after=107)
noninsfrees <- append(noninsfrees, noninsfree$SL1344_5, after=65)
noninsfrees <- append(noninsfrees, noninsfree$SL1344_7, after=56)
noninsfrees <- append(noninsfrees, noninsfree$SL1344_3, after=15)
noninsfrees <- append(noninsfrees, noninsfree$SL1344_1, after=14)
insfrees <- length(iitable$locus.tag) - noninsfrees

pdf('~/EnTrI/figures/false-positive-rate_density.pdf')

par(mar=c(6.1,5.1,4.1,2.1))
plx<-predict(loess(fprs ~ 1/dens, span = 0.2), se=T)
plot(1/dens,fprs, type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = 'False positive rate', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,0.050)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)
# axis(side = 3, cex.axis = 2, labels = seq(5e5,0,-1e5), at=seq(-5e5,0,1e5))
# mtext(side = 3, line = 3, '# insertion sites', cex = 2)

# dat <- data.frame(1/dens, fprs)
# colnames(dat) <- c('dns', 'fprs')
# print(ggplot(dat,aes(dns,fprs))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("False positive rate")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black"))
#       )


par(mar=c(6.1,5.1,4.1,2.1))
plx<-predict(loess(tprs ~ 1/dens, span = 0.2), se=T)
plot(1/dens,tprs, type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = 'True positive rate', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,0.050)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)

# dat <- data.frame(1/dens, tprs)
# colnames(dat) <- c('dns', 'tprs')
# print(ggplot(dat,aes(dns,tprs))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("True positive rate")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) + 
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black")))

# plx<-predict(loess(fprsunif ~ ins, span = 0.2), se=T)
# points(ins, fprsunif)
# lines(ins,plx$fit, col=3)

plx<-predict(loess(essens ~ 1/dens, span = 0.2), se=T)
plot(1/dens,essens, type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = '# predicted essential genes', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,430)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)

# dat <- data.frame(1/dens, essens)
# colnames(dat) <- c('dns', 'essens')
# print(ggplot(dat,aes(dns,essens))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("# predicted essential genes")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) + 
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black")))

plx<-predict(loess(mccs ~ 1/dens, span = 0.2), se=T)
plot(1/dens,mccs, type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = 'Matthews correlation coefficient', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,430)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)

# plx<-predict(loess(essensunif ~ ins, span = 0.2), se=T)
# points(ins,essensunif)
# lines(ins,plx$fit, col=3)

# dat <- data.frame(1/dens, mccs)
# colnames(dat) <- c('dns', 'mccs')
# print(ggplot(dat,aes(dns,fprs))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("Matthews correlation coefficient")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black")))

plx<-predict(loess(insfrees ~ 1/dens, span = 0.2), se=T)
plot(1/dens,insfrees, type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = '# insertion free genes', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,430)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)

# dat <- data.frame(1/dens, insfrees)
# colnames(dat) <- c('dns', 'insfrees')
# print(ggplot(dat,aes(dns,insfrees))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("# insertion free genes")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black")))

plx<-predict(loess(noninsfrees*100/length(iitable$locus.tag) ~ 1/dens, span = 0.2), se=T)
plot(1/dens,noninsfrees*100/length(iitable$locus.tag), type = 'p', #ylim=c(0,0.1),
     col=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "#8c510a", "#d8b365"),
     pch=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), 17, 20),
     xlab = 'Insertion density', ylab = '% genes with insertion(s)', cex.lab = 2, cex.axis = 2, cex=1.5#, ylim=c(0,430)
)
lines(1/dens,plx$fit, lwd=2, col="#01665e")
# lines(1/dens,plx$fit + qt(0.025,plx$df)*plx$se, lty=2, lwd=2)
# lines(1/dens,plx$fit + qt(0.975,plx$df)*plx$se, lty=2, lwd=2)
x <- c(1/dens, rev(1/dens))
y <- c(plx$fit + qt(0.025,plx$df)*plx$se, rev(plx$fit + qt(0.975,plx$df)*plx$se))
color <- adjustcolor("#c7eae5",alpha.f=0.5)
polygon(x,y,col=color, border = NA)
axis(side =3, at=1/seq(10,300,10), labels=seq(10,300,10), cex.axis=2)
mtext(side = 3, line = 2.5, 'Insertion resolution', cex = 2)

# dat <- data.frame(1/dens, noninsfrees*100/length(iitable$locus.tag))
# colnames(dat) <- c('dns', 'noninsfree')
# print(ggplot(dat,aes(dns,noninsfree))+
#         geom_point(color=ifelse(dens %in% c(den$SL1344_1, den$SL1344_3, den$SL1344_7, den$SL1344_5, den$SL1344_9), "red", "midnightblue"), size=3) +
#         geom_smooth(fill='gray47', color='black', method = "loess", span=0.2)+
#         xlab("Insertion density") +ylab("% genes with insertion(s)")+ theme_bw()+
#         theme(plot.title = element_text(size=24, hjust=0.5, face = 'bold'), axis.text=element_text(size=18)) + theme(axis.title.x = element_text(size=22)) +
#         theme(axis.title.y = element_text(size=22))+ guides(fill=FALSE) +
#         theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
#               panel.background = element_blank(), axis.line = element_line(colour = "black")))

dev.off()