library("dbscan")
library("ROCR")
library("ggplot2")
getmode <- function(v) {
  uniqv <- unique(v)
  uniqv[which.max(tabulate(match(v, uniqv)))]
}
cols <- c("#01665e", "#5ab4ac", "#8c510a", "#d8b365", "#c51b7d")
notrimpath = '../results/biases/check-biases_not-trimmed/BW25113.txt'
trimpath = '../results/biases/check-biases/BW25113.txt'
infotab <- read.table(notrimpath, header = FALSE, stringsAsFactors = FALSE)
trimtab <- read.table(trimpath, header = FALSE, stringsAsFactors = FALSE)
colnames(infotab) <- c("name", "ii", "dist", "gc", "length", "pos")
colnames(trimtab) <- c("name", "ii", "essentiality", "dist", "gc", "length", "pos")
labpath <- '../results/ecogenecounterparts/BW25113.txt'
essentiality <- read.table(labpath, header = FALSE, stringsAsFactors = FALSE)
sp = 0.2

biasestable <- as.data.frame(matrix(0, ncol = 6, nrow = length(infotab$ii)))
colnames(biasestable) <- c("expected", "none", "dist", "gc", "pos", "all")
biasestable$expected <- essentiality$V2
labels <- essentiality$V2
row.names(biasestable) <- essentiality$V1

ii <- infotab$ii
res <- dbscan(as.matrix(infotab$ii), minPts = 200, eps = 0.05)
ess <- res$cluster[which.min(infotab$ii)]
for (i in (1:length(infotab$ii)))
{
  if (res$cluster[i] == ess)
  {
    biasestable$none[i] = 1
  } else
  {
    biasestable$none[i] = 0
  }
}

fit <- loess(infotab$ii~infotab$dist, span = sp)
loessprediction <- predict(fit, infotab$dist)
ii_dnormalised = infotab$ii/(loessprediction/mean(infotab$ii))
res <- dbscan(as.matrix(ii_dnormalised), minPts = 200, eps = 0.05)
ess <- res$cluster[which.min(ii_dnormalised)]
for (i in (1:length(ii_dnormalised)))
{
  if (res$cluster[i] == ess)
  {
    biasestable$dist[i] = 1
  } else
  {
    biasestable$dist[i] = 0
  }
}

fit <- loess(infotab$ii~infotab$gc, span = sp)
loessprediction <- predict(fit, infotab$gc)
ii_gcnormalised = infotab$ii/(loessprediction/mean(infotab$ii))
res <- dbscan(as.matrix(ii_gcnormalised), minPts = 200, eps = 0.05)
ess <- res$cluster[which.min(ii_gcnormalised)]
for (i in (1:length(ii_gcnormalised)))
{
  if (res$cluster[i] == ess)
  {
    biasestable$gc[i] = 1
  } else
  {
    biasestable$gc[i] = 0
  }
}

ii_trimmed <- trimtab$ii
res <- dbscan(as.matrix(trimtab$ii), minPts = 200, eps = 0.05)
ess <- res$cluster[which.min(trimtab$ii)]
for (i in (1:length(trimtab$ii)))
{
  if (res$cluster[i] == ess)
  {
    biasestable$pos[i] = 1
  } else
  {
    biasestable$pos[i] = 0
  }
}

allpath <- '../results/biases/dbscan/BW25113.txt'
alltab <- read.table(allpath, header = FALSE, stringsAsFactors = FALSE)
ii_total <- alltab$V2
res <- dbscan(as.matrix(ii_total), minPts = 200, eps = 0.05)
ess <- res$cluster[which.min(ii_total)]
for (i in (1:length(ii_total)))
{
  if (res$cluster[i] == ess)
  {
    biasestable$all[i] = 1
  } else
  {
    biasestable$all[i] = 0
  }
}

pdf('../figures/compare-bias-predictions.pdf')
par(mar=c(5.1,5.1,4,2), cex.axis=2)
predii <- prediction(-ii, labels)
perfii <- performance(predii,"tpr","fpr")
aucii <- performance(predii,measure = "auc")@y.values[[1]]
plot(perfii,col=cols[1],lty=1,lwd=4,cex.lab=2,cex.axis=2, ylim=c(0.7,1), cex.main=2,
     main = 'Essentiality predictors in E. coli BW25113')

prediid <- prediction(-ii_dnormalised, labels)
perfiid <- performance(prediid,"tpr","fpr")
auciid <- performance(prediid,measure = "auc")@y.values[[1]]
plot(perfiid,col=cols[2],lty=1,lwd=4,add=TRUE)

prediigc <- prediction(-ii_gcnormalised, labels)
perfiigc <- performance(prediigc,"tpr","fpr")
auciigc <- performance(prediigc,measure = "auc")@y.values[[1]]
plot(perfiigc,col=cols[3],lty=1,lwd=4,add=TRUE)

prediitr <- prediction(-ii_trimmed, labels)
perfiitr <- performance(prediitr,"tpr","fpr")
auciitr <- performance(prediitr,measure = "auc")@y.values[[1]]
plot(perfiitr,col=cols[4],lty=1,lwd=4,add=TRUE)

prediito <- prediction(-ii_total, labels)
perfiito <- performance(prediito,"tpr","fpr")
auciito <- performance(prediito,measure = "auc")@y.values[[1]]
plot(perfiito,col=cols[5],lty=1,lwd=4,add=TRUE)

lab = c(paste('Unnorm. AUC = ', format(round(aucii, 4), nsmall = 4)),
        paste('Distance norm. AUC = ', format(round(auciid, 4), nsmall = 4)),
        paste('GC norm. AUC = ', format(round(auciigc, 4), nsmall = 4)),
        paste('Trimmed genes AUC = ', format(round(auciitr, 4), nsmall = 4)),
        paste('Norm. for all AUC = ', format(round(auciito, 4), nsmall = 4)))
legend("bottomright", inset=.05, legend=lab, lwd=4, col=cols, cex = 1.5)
dev.off()