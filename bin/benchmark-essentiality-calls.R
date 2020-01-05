library("ROCR")
library(stringr)
library("mclust")
library(mixtools)
library("plotrix")

locus = c('BN373','ENC','ROD','SL1344','STMMW','t','ERS227112','NCTC13441','SEN','SL3261','STM', "EC958", "BW25113")
address = c('Klebsiella_pneumoniae_subsp_pneumoniae_Ecl8_HF536482_v1.fasta','Enterobacter_cloacae_subsp_cloacae_NCTC_9394_v1.fasta',
            'Citrobacter_rodentium_ICC168_FN543502_v1.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhimurium_SL1344_FQ312003_v4.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhimurium_str_D23580_v1.fasta',
            'Salmonella_enterica_subsp_enterica_serovar_Typhi_Ty2_v1.fasta', 'Klebsiella_pneumoniae_RH201207_v0.fasta',
            'Escherichia_coli_UPEC_ST131_chromosome_v0.fasta','P125109.fasta',
            'SL3261.fasta','Salmonella_enterica_subsp_enterica_serovar_Typhimurium_A130_v0.fasta', 'HG941718.fasta',
            'CP009273.fasta')
names = c("ROD", "ENC", "NCTC13441", "ERS227112", "BN373", "SEN", "STM", "SL1344", "STMMW", "t", "SL3261", "BW25113", "EC958")
dict = c("C. rodentium ICC168", "E. cloacae NCTC 9394", "E. coli UPEC ST131",
         "K. pneumoniae RH201207", "K. pneumoniae Ecl8", "S. Enteritidis", "S. Typhimurium A130",
         "S. Typhimurium SL1344", "S. Typhimurium D23580", "S. Typhi Ty2", "S. Typhimurium SL3261",
         "E. coli BW25113", "E. coli UPEC ST131 EC958")
names(dict) <- names

contingency <- array(0, dim=c(2,2,length(locus)))
dimnames(contingency)[[3]]=locus
dimnames(contingency)[[2]]=c('Essential', 'Non-essential')
dimnames(contingency)[[1]]=c('Essential', 'Non-essential')
names(dimnames(contingency))=c('Real', 'Predicted', 'Bacterium')

outdir = '../results/maximise_MCC/'
dir.create(outdir)
outdir_montecarlo = paste(outdir, 'monte-carlo/', sep='')
dir.create(outdir_montecarlo)
outdir_pca = paste(outdir, 'pca/', sep='')
dir.create(outdir_pca)
outdir_pcaeq = paste(outdir, 'pca-eq.txt', sep='')
file.create(outdir_pcaeq)
# colors=c('blue', 'darkslategrey', 'limegreen', 'red', 'cyan', 'black', 'orange', 'purple', 'gray', 'brown', 'goldenrod4')
colors=c("#c51b7d", "#e9a3c9", "#8c510a", "#d8b365", "#01665e", "#5ab4ac")
avgaucii=c(0,0)
avgaucmc=c(0,0)
avgaucconz=0
avgaucmeandist=0
avgaucpca=0
avgaucpca2=0
ranks=c(0,0,0,0,0,0)
for (i in seq(length(locus)))
{
  real_old = read.table(paste('../results/ecogenecounterparts/',locus[i],'.txt',sep = ''), as.is=TRUE, header=FALSE, sep="\t")
  names(real_old) <- c('gene', 'essentiality')
  
  # montecarlo = read.table(paste('../results/monte-carlo/',locus[i],'.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  # montecarlo = montecarlo[,c(2,3,4,5,6)]
  # montecarlo[,1] = -montecarlo[,1]
  # montecarlo[,3] = -montecarlo[,3]
  # montecarlo[,5] = -montecarlo[,5]
  # names(montecarlo) <- c('-DESeqPval','DESeqLFC','-LFC','DESeqLFCDist', '-LFCDist')
  montecarlo = read.table(paste('../results/Tn-seq/',locus[i],'.out.DESeq.tsv',sep=''), as.is=TRUE, header=TRUE, sep="\t")
  montecarlo = montecarlo[,c(7,6,2)]
  montecarlo[,2] = -montecarlo[,2]
  montecarlo[,3] = -montecarlo[,3]
  montecarlo = montecarlo[montecarlo$id %in% real_old$gene,]
  real_new = real_old[real_old$gene %in% montecarlo$id,]
  montecarlo = montecarlo[,c(2,3)]
  ## montecarlo$logfoldchange = -montecarlo$logfoldchange
  
  biotradis = read.table(paste('../results/insertion-indices/gamma/', locus[i], '.txt',sep=''), as.is=TRUE, header=FALSE, sep="\t")
  biotradis = biotradis[,c(2,4)]
  biotradis = -biotradis
  biotradis = biotradis[real_old$gene %in% real_new$gene,]
  
  pdf(paste('../figures/essential-call-comparison-', locus[i], '.pdf', sep=''))
  
  aucbiotradis = c()
  cutoffbiotradis = c()
  for (j in seq(2))
  {
    predbiotradis <- prediction(biotradis[,j], real_new$essentiality)
    perfbiotradis <- performance(predbiotradis,"tpr","fpr")
    aucbiotradis <- c(aucbiotradis, performance(predbiotradis,measure = "auc")@y.values[[1]])
    if (j==2)
    {
      # plot(perfbiotradis,col=colors[j],lty=1,lwd=4,cex.lab=2,cex.axis=2, cex.main=2, add=TRUE)
    }
    else
    {
      mar.default <- c(5,4,4,2) + 0.1
      par(mar = mar.default + c(0, 1, 0, 0), cex.axis=2, family='sans')
      plot(perfbiotradis,col=colors[1],lty=1,lwd=4,cex.lab=2, main=dict[locus[i]], cex.main=2, ylim=c(0.6,1))
    }
    
    perfbiotradismcc <- performance(predbiotradis,'mat')
    maxmcc = max(perfbiotradismcc@y.values[[1]][!is.na(perfbiotradismcc@y.values[[1]])])
    cutoff = perfbiotradismcc@x.values[[1]][perfbiotradismcc@y.values[[1]]==maxmcc & !is.na(perfbiotradismcc@y.values[[1]])]
    cutoffbiotradis = c(cutoffbiotradis, cutoff)
    avgaucii[j] = avgaucii[j] + aucbiotradis[j]
  }
  
  aucmontecarlo = c()
  cutoffmontecarlo = c()
  # for (j in seq(5))
  for (j in seq(2))
  {
    predmontecarlo <- prediction(montecarlo[,j], real_new$essentiality)
    perfmontecarlo <- performance(predmontecarlo,"tpr","fpr")
    aucmontecarlo <- c(aucmontecarlo, performance(predmontecarlo,measure = "auc")@y.values[[1]])
    if (j==2)
    {
      plot(perfmontecarlo,col=colors[2],lty=1,lwd=4,cex.lab=2,cex.axis=2, cex.main=2, add=TRUE)
    }
    
    perfmontecarlomcc <- performance(predmontecarlo,'mat')
    maxmcc = max(perfmontecarlomcc@y.values[[1]][!is.na(perfmontecarlomcc@y.values[[1]])])
    cutoff = perfmontecarlomcc@x.values[[1]][perfmontecarlomcc@y.values[[1]]==maxmcc & !is.na(perfmontecarlomcc@y.values[[1]])]
    cutoffmontecarlo = c(cutoffmontecarlo,cutoff)
    avgaucmc[j] = avgaucmc[j] + aucmontecarlo[j]
  }
  
  fastas_dir <- paste("~/EnTrI/data/fasta-protein/chromosome",address[i],sep='/')
  plots_dir <- paste("~/EnTrI/data/plot-files/chromosome/",locus[i],".plot",sep="")
  plots = list()
  sumlength = list()
  plotfile = as.matrix(read.table(plots_dir, as.is=TRUE))
  locusid = strsplit(basename(plots_dir),"\\.")[[1]][1]
  plots[[locusid]] = plotfile[,1] + plotfile[,2]
  
  
  #counttable <- data.frame()
  consecutivezeros = c()
  meandist = c()
  fastafile = readLines(fastas_dir)
  iitable = list()
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
        locusid = str_match(locustag, "([[:graph:]]+)\\_+[[:alnum:]]+")[2]
        if (is.na(locusid))
        {
          locusid = str_match(locustag, "([[:alpha:]]+)[[:digit:]]+")[2]
        }
        if (!(is.null(plots[[locusid]])))
        {
          len = end - start + 1
          consecutives = rle(plots[[locusid]][start:end])
          # meanconsecutives = mean(consecutives$lengths[consecutives$values==0])
          meanconsecutives = len / (sum(plots[[locusid]][start:end]>0)+1)
          maxconsecutivezeros = max(consecutives$lengths[consecutives$values==0]) / len
          meandist = c(meandist, meanconsecutives)
          consecutivezeros = c(consecutivezeros, maxconsecutivezeros)
        }
      }
    }
  }
  consecutivezeros = consecutivezeros[real_old$gene %in% real_new$gene]
  predconz <- prediction(consecutivezeros, real_new$essentiality)
  perfconz <- performance(predconz,"tpr","fpr")
  aucconz <- performance(predconz,measure = "auc")@y.values[[1]]
  # plot(perfconz,col=colors[8],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  plot(perfconz,col=colors[3],lty=1,lwd=4,cex.lab=2,cex.axis=2,cex.main=2, add=TRUE)
  
  perfconzmcc <- performance(predconz,'mat')
  maxmcc = max(perfconzmcc@y.values[[1]][!is.na(perfconzmcc@y.values[[1]])])
  cutoff = perfconzmcc@x.values[[1]][perfconzmcc@y.values[[1]]==maxmcc & !is.na(perfconzmcc@y.values[[1]])]
  cutoffconz = cutoff
  avgaucconz = avgaucconz + aucconz
  
  meandist = meandist[real_old$gene %in% real_new$gene]
  predmeandist <- prediction(meandist, real_new$essentiality)
  perfmeandist <- performance(predmeandist,"tpr","fpr")
  aucmeandist <- performance(predmeandist,measure = "auc")@y.values[[1]]
  # plot(perfmeandist,col=colors[9],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  plot(perfmeandist,col=colors[4],lty=1,lwd=4,cex.lab=2,cex.axis=2, cex.main=2, add=TRUE)
  
  perfmeandistmcc <- performance(predmeandist,'mat')
  maxmcc = max(perfmeandistmcc@y.values[[1]][!is.na(perfmeandistmcc@y.values[[1]])])
  cutoff = perfmeandistmcc@x.values[[1]][perfmeandistmcc@y.values[[1]]==maxmcc & !is.na(perfmeandistmcc@y.values[[1]])]
  cutoffmeandist = cutoff
  avgaucmeandist = avgaucmeandist + aucmeandist
  
  data = cbind(biotradis$V2, montecarlo$log2FoldChange, consecutivezeros, meandist)
  data = apply(data, 2, function(x){(x-mean(x))/sd(x-mean(x))})
  # for (j in seq(ncol(data)))
  # {
  #   data[,j]=(data[,j]-mean(data[,j]))/sd(data[,j]-mean(data[,j]))
  # }
  data.pca <- prcomp(data, center = TRUE, scale. = TRUE)
  if (data.pca$rotation[1,1] < 0)
  {
    data.pca$x = -data.pca$x
  }
  
  write(paste(locus[i], '\n', 'PCA-all', sep=''),file=outdir_pcaeq,append=TRUE)
  write('Eigen vectors:',file=outdir_pcaeq,append=TRUE)
  write(data.pca$rotation[1,],file=outdir_pcaeq,append=TRUE)
  write(data.pca$rotation[2,],file=outdir_pcaeq,append=TRUE)
  write(data.pca$rotation[3,],file=outdir_pcaeq,append=TRUE)
  write(data.pca$rotation[4,],file=outdir_pcaeq,append=TRUE)
  write('Standard deviation:',file=outdir_pcaeq,append=TRUE)
  write(data.pca$sdev,file=outdir_pcaeq,append=TRUE)
  
  predpca <- prediction(data.pca$x[,1], real_new$essentiality)
  perfpca <- performance(predpca,"tpr","fpr")
  aucpca <- performance(predpca,measure = "auc")@y.values[[1]]
  # plot(perfpca,col=colors[10],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  plot(perfpca,col=colors[5],lty=1,lwd=4,cex.lab=2,cex.axis=2, cex.main=2, add=TRUE)
  
  perfpcamcc <- performance(predpca,'mat')
  maxmcc = max(perfpcamcc@y.values[[1]][!is.na(perfpcamcc@y.values[[1]])])
  cutoff = perfpcamcc@x.values[[1]][perfpcamcc@y.values[[1]]==maxmcc & !is.na(perfpcamcc@y.values[[1]])]
  cutoffpca = cutoff
  avgaucpca = avgaucpca + aucpca
  
  data2 = cbind(biotradis$V2, consecutivezeros, meandist)
  data2 = apply(data2, 2, function(x){(x-mean(x))/sd(x-mean(x))})
  # for (j in seq(ncol(data)))
  # {
  #   data[,j]=(data[,j]-mean(data[,j]))/sd(data[,j]-mean(data[,j]))
  # }
  data2.pca <- prcomp(data2, center = TRUE, scale. = TRUE)
  if (data2.pca$rotation[1,1] < 0)
  {
    data2.pca$x = -data2.pca$x
  }
  
  write(paste('\n', 'PCA-ex-montecarlo', sep=''),file=outdir_pcaeq,append=TRUE)
  write('Eigen vectors:',file=outdir_pcaeq,append=TRUE)
  write(data2.pca$rotation[1,],file=outdir_pcaeq,append=TRUE)
  write(data2.pca$rotation[2,],file=outdir_pcaeq,append=TRUE)
  write(data2.pca$rotation[3,],file=outdir_pcaeq,append=TRUE)
  write('Standard deviation:',file=outdir_pcaeq,append=TRUE)
  write(data2.pca$sdev,file=outdir_pcaeq,append=TRUE)
  write('\n\n',file=outdir_pcaeq,append=TRUE)
  
  predpca2 <- prediction(data2.pca$x[,1], real_new$essentiality)
  perfpca2 <- performance(predpca2,"tpr","fpr")
  aucpca2 <- performance(predpca2,measure = "auc")@y.values[[1]]
  # plot(perfpca,col=colors[10],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  plot(perfpca2,col=colors[6],lty=1,lwd=4,cex.lab=2,cex.axis=2, cex.main=2, add=TRUE)
  
  perfpcamcc2 <- performance(predpca2,'mat')
  maxmcc = max(perfpcamcc2@y.values[[1]][!is.na(perfpcamcc2@y.values[[1]])])
  cutoff = perfpcamcc2@x.values[[1]][perfpcamcc2@y.values[[1]]==maxmcc & !is.na(perfpcamcc2@y.values[[1]])]
  cutoffpca2 = cutoff
  avgaucpca2 = avgaucpca2 + aucpca2
  
  # data.sum = rowSums(data)
  # 
  # predsum <- prediction(data.sum, real$essentiality)
  # perfsum <- performance(predsum,"tpr","fpr")
  # aucsum <- performance(predsum,measure = "auc")@y.values[[1]]
  # plot(perfsum,col=colors[10],lty=1,lwd=4,cex.lab=1.5,xaxis.cex.axis=1.7,yaxis.cex.axis=1.7, add=TRUE)
  # 
  # perfsummcc <- performance(predsum,'mat')
  # maxmcc = max(perfsummcc@y.values[[1]][!is.na(perfsummcc@y.values[[1]])])
  # cutoff = perfsummcc@x.values[[1]][perfsummcc@y.values[[1]]==maxmcc & !is.na(perfsummcc@y.values[[1]])]
  # cutoffsum = cutoff
  
  labels <- c(paste("Insertion index, AUC = ", format(round(aucbiotradis[1], 4), nsmall = 4))
              # , paste("BioTraDIS logodds, AUC = ", format(round(aucbiotradis[2], 4), nsmall = 4))
              # , paste("Monte Carlo Pval, AUC = ", format(round(aucmontecarlo[1], 4), nsmall = 4))
              , paste("Monte Carlo DESeq, AUC = ", format(round(aucmontecarlo[2], 4), nsmall = 4))
              # , paste("Monte Carlo LFC, AUC = ", format(round(aucmontecarlo[3], 4), nsmall = 4))
              # , paste("Monte Carlo DESeq LFC Distances, AUC = ", format(round(aucmontecarlo[4], 4), nsmall = 4))
              # , paste("Monte Carlo LFC Distances, AUC = ", format(round(aucmontecarlo[5], 4), nsmall = 4))
              , paste("Largest uninterrupted fraction, AUC = ", format(round(aucconz, 4), nsmall = 4))
              , paste("Mean distance between inserts, AUC = ", format(round(aucmeandist, 4), nsmall = 4))
              , paste("PCA, AUC = ", format(round(aucpca, 4), nsmall = 4))
              # , paste("SUM, AUC = ", format(round(aucsum, 4), nsmall = 4))
              , paste("PCA without Monte Carlo, AUC = ", format(round(aucpca2, 4), nsmall = 4))
              )
  legend("bottomright", inset=.05, labels, lwd=2, col=colors)
  ranks = ranks + c(aucbiotradis[1],aucmontecarlo[2],aucconz,aucmeandist,aucpca,aucpca2)
  
  plot(-biotradis[,2],montecarlo$log2FoldChange, pch=20, col=real_new$essentiality+1, xlab = "BioTraDIS logodds", ylab = "Monte Carlo logfoldchange")
  labels <- c("Essential","Non-essential")
  legend("topright", inset=.05, labels, pch=20, col=c("red","black"))
  
  plot(-biotradis[,2],log2(-montecarlo$padj), pch=20, col=real_new$essentiality+1, xlab = "BioTraDIS logodds", ylab = "log2(Monte Carlo P-value)")
  labels <- c("Essential","Non-essential")
  legend("bottomright", inset=.05, labels, pch=20, col=c("red","black"))
  
  # plot(-montecarlo$DESeqLFC,-montecarlo$DESeqLFCDist, pch=20, col=real$essentiality+1, xlab = "Monte Carlo DESeq LFC", ylab = "Monte Carlo DESeq LFC Distances")
  # labels <- c("Essential","Non-essential")
  # legend("bottomright", inset=.05, labels, pch=20, col=c("red","black"))
  
  # plot(2^(-montecarlo$DESeqLFC),2^(-montecarlo$DESeqLFCDist), pch=20, col=real$essentiality+1, xlab = "Monte Carlo DESeq FC", ylab = "Monte Carlo DESeq FC Distances")
  # labels <- c("Essential","Non-essential")
  # legend("topright", inset=.05, labels, pch=20, col=c("red","black"))
  
  hist(montecarlo$log2FoldChange, breaks=150, xlab = "Monte Carlo logFC", main = "Histogram of Monte Carlo logFC")
  hist(log2(-montecarlo$padj), breaks=150, xlab = "log2(Monte Carlo P-value)", main = "Histogram of Monte Carlo P-value")
  hist(log2(consecutivezeros+1e-5), breaks = 150, xlab = "log2(Largest uninterrupted fraction+1e-5)", main = "Histogram of largest uninterrupted fraction")
  hist(log2(meandist+1e-5), breaks = 150, xlab = "log2(Mean distance between inserts+1e-5)", main = "Histogram of mean distance between inserts")
  
  
  # essentiality = ifelse(montecarlo$DESeqLFC >= cutoffmontecarlo[2], 'essential', ifelse(montecarlo$DESeqLFC <= 
  #                                                                                         mean(montecarlo$DESeqLFC)-0.8*cutoffmontecarlo[2],
  #                                                                                       'beneficial-loss', 'non-essential'))
  essentiality = ifelse(montecarlo$log2FoldChange >= cutoffmontecarlo[2], 'essential', ifelse(montecarlo$log2FoldChange > 0 & -montecarlo$padj< 0.01,
                                                                                        'beneficial-loss', 'non-essential'))
  to_print = cbind(real_new$gene, montecarlo$log2FoldChange, essentiality)
  outpath = paste(outdir_montecarlo, locusid, ".txt", sep="")
  write.table(to_print, file=outpath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  
  tt = length(real_new$gene[montecarlo$log2FoldChange >= cutoffmontecarlo[2] & real_new$essentiality == "1"])
  ft = length(real_new$gene[montecarlo$log2FoldChange >= cutoffmontecarlo[2] & real_new$essentiality == "0"])
  ff = length(real_new$gene[montecarlo$log2FoldChange <= cutoffmontecarlo[2] & real_new$essentiality == "0"])
  tf = length(real_new$gene[montecarlo$log2FoldChange <= cutoffmontecarlo[2] & real_new$essentiality == "1"])
  contingency[,,locus[i]]=rbind(c(tt,tf),c(ft,ff))
  # if (locus[i] == 'CS17')
  # {
  #   cs17fps = real$gene[montecarlo$DESeqLFC >= cutoffmontecarlo[2] & real$essentiality == "0"]
  #   write.table(cs17fps, file='../results/cs17fps.txt
  #               ', quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  # }
  
  # mod = Mclust(montecarlo$DESeqLFC, G=1:2)
  # plot(mod, what = "classification")
  
  # normalmix <- normalmixEM(montecarlo$DESeqLFC, k=3)
  # print(sum(apply(normalmix$posterior, 1, function(x) min(which(x == max(x, na.rm = TRUE))))==1))
  # print(sum(apply(normalmix$posterior, 1, function(x) min(which(x == max(x, na.rm = TRUE))))==2))
  # print(sum(apply(normalmix$posterior, 1, function(x) min(which(x == max(x, na.rm = TRUE))))==3))
  # print(cutoffmontecarlo[4])
  ############################# plot
  #hist(data.pca$x[,1], breaks = 200)
  print(locus[i])
  ################# Stoufer's method
  # zscores=(data[,1]*abs(data.pca$rotation[1,1])+data[,2]*abs(data.pca$rotation[2,1])+data[,3]*abs(data.pca$rotation[3,1]))/
  #   sqrt(sum(data.pca$rotation[1,1]^2+data.pca$rotation[2,1]^2+data.pca$rotation[3,1]^2)) #It is actually equal to data.pca$x[,1]
  normalisedpca = (data.pca$x[,1]- mean(data.pca$x[,1]))/sd(data.pca$x[,1]- mean(data.pca$x[,1]))
  # hist(normalisedpca, breaks = 200)
  pvalue2sided=2*pnorm(-abs(normalisedpca))
  #print(length(real$gene[pvalue2sided<=0.05 & data.pca$x[,1]>0]))
  #print(length(real$gene[pvalue2sided<=0.05 & data.pca$x[,1]<0]))
  print((cutoffpca- mean(data.pca$x[,1]))/sd(data.pca$x[,1]- mean(data.pca$x[,1])))
  pcacutoff = 1.644854 #pnorm(1.644854) = 0.5 #The average of all pcacutoffs defined by maximising MCC is 1.650449
  print(length(real_new$gene[normalisedpca>pcacutoff]))
  print(length(real_new$gene[normalisedpca< -pcacutoff]))
  pca_essentiality = ifelse(normalisedpca >= pcacutoff, 'essential', ifelse(normalisedpca < -pcacutoff,
                                                                                        'beneficial-loss', 'non-essential'))
  pca_print = cbind(real_new$gene, data.pca$x[,1], pca_essentiality)
  # pca_print = cbind(real$gene, normalisedpca, pca_essentiality)
  pcapath = paste(outdir_pca, locusid, ".txt", sep="")
  write.table(pca_print, file=pcapath, quote = FALSE, sep = "\t", col.names = FALSE, row.names = FALSE)
  # #intercept=normalisedpca[normalisedpca>0 & pvalue2sided==max(pvalue2sided[normalisedpca>0 & pvalue2sided<0.05])]
  # intercept=1.5
  # abline(v=intercept,col='red')
  # #intercept=normalisedpca[normalisedpca<0 & pvalue2sided==max(pvalue2sided[normalisedpca<0 & pvalue2sided<0.05])]
  # intercept=-1.5
  # abline(v=intercept,col='red')
  
  h <- hist(normalisedpca, breaks = 200, plot=FALSE)
  cuts <- cut(h$breaks, c(-Inf,-pcacutoff, pcacutoff, Inf))
  plot(h, col=c("darkgoldenrod4", "turquoise4", "darkmagenta")[cuts], xlab = "NPEQ", main =dict[locus[i]], cex.lab = 2,
       cex.axis = 2, cex.main = 2, lty= "blank")
  legend(1,180, c("Essential","Non-essential", "Beneficial loss"), lty=c(1,1,1), lwd=c(4,4,4),cex=1.15,
         col=c("darkmagenta","turquoise4", "darkgoldenrod4"), bty="n")
  
  ############ qq plot
  #qq=qqnorm(data.pca$x[,1])
  #qqline(data.pca$x[,1], col = 2)
  #dist = qq$x-qq$y
  #print(sum(dist>0.5))
  #print(sum(dist< -0.5))
  # nnn=(data.pca$x[,1]-mean(data.pca$x[,1]))/sd(data.pca$x[,1]-mean(data.pca$x[,1]))
  # hist(nnn,breaks=200)
  # print(sum(nnn>1.5))
  # print(sum(nnn < -1.5))
  
  dev.off()
}

pdf('../figures/essentiality-call-accuracy.pdf')
fourfoldplot(contingency, conf.level=0, std='i', space=0.2)
dev.off()

avgaucii = avgaucii / length(locus)
avgaucmc = avgaucmc / length(locus)
avgaucconz = avgaucconz / length(locus)
avgaucmeandist = avgaucmeandist / length(locus)
avgaucpca = avgaucpca / length(locus)
avgaucpca2 = avgaucpca2 / length(locus)
pdf('../figures/average-auc.pdf')
par(las=2, mar= mar.default + c(8,0,0,0))
barplot(c(avgaucii[1],avgaucmc[2],avgaucconz,avgaucmeandist,avgaucpca, avgaucpca2),ylim=c(0.95,0.961),names.arg = c('Insertion Index',
                                                                                                                   'Monte Carlo DESeq',
                                                                                                 'Largest Uninterrupted Fraction',
                                                                                                 'Mean Distance', 'PCA',
                                                                                                 'PCA Excluding Monte Carlo'), xpd = FALSE,
        main = "Average AUC")
dev.off()

############### compare TnSeq:
# tnseqtest = read.csv('~/program-bank/TnSeq/example data&code/ROD-essentiality.csv')
# library(dplyr)
# tnseqtest$Gene=gsub(" ", "", tnseqtest$Gene)
# real2=semi_join(real,tnseqtest, by=c("gene"="Gene"))
# tnseqtest2=semi_join(tnseqtest, real2, by=c("Gene"="gene"))
# predtnseq <- prediction(-tnseqtest2[,1], real2$essentiality)
# perftnseq <- performance(predtnseq,"tpr","fpr")
# auctnseq <- performance(predtnseq,measure = "auc")@y.values[[1]]

print(ranks)