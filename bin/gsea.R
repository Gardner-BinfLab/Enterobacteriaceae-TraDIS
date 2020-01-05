#!/usr/bin/env Rscript

# num_genomes = 3803
# threshold = num_genomes/2
#threshold=0
library("getopt")

opt = getopt(matrix( c('help', 'h', 0, "logical", 
                       'verbose', 'v', 0, "integer",
                       'database', 'd', 1, "character",
                       'genelist', 'l', 1, "character",
                       'sigthresh', 's', 1, "integer",
                       'outdir', 'o', 1, "character"
), ncol=4, byrow=TRUE ) );

if(! is.null(opt$help) || is.null(opt$genelist )  || is.null(opt$database ) )
{
  cat(paste("Usage: gsea.R [-h] [-v] [-e] -d database.txt -l genelist.txt [-s int] [-o outdir]\n"));
  q(status=1);
}

if(is.null(opt$outdir)){	opt$outdir <- getwd()}
if(is.null(opt$sigthresh)){	opt$sigthresh <- 3}


path_db <- read.table(opt$database, header=TRUE, sep="\t")
genelist <- read.table(opt$genelist, header=FALSE, sep="\t")
genelist[,2] = genelist[,2] - median(genelist[,2])
genelist <- genelist[order(genelist[,2], decreasing=TRUE), ]

hash <- function( keys ) {
    result <- new.env( hash = TRUE, parent = emptyenv(), size = length( keys ) )
    for( key in keys ) {
        result[[ key ]] <- NA
    }
    return( result )
}

paths <- unique(path_db$path_id)

length <- length(genelist[,1])

lines <- c()

pval_hash <- hash(paths)
presence <- hash(paths)

for ( i in paths){
  this_path <- path_db[grepl(i, path_db[,2]),]
  labelled_list <- cbind(genelist, sapply(genelist[,1], function (entry) (entry %in% this_path[,1])))
  labelled_list <- labelled_list[order(labelled_list[,2], decreasing=TRUE),]
  n_lab <- sum(labelled_list[,3])
  n_unlab <- length - n_lab
  
  max_p <- 0
  max_m <- 0
  max_p_stat <- 0
  max_p_g <- 0
  max_p_g_up <- 0
  max_m_stat <- 0
  max_m_g <- 0
  max_m_g_down <- 0
  vals <- c()
  for (j in 1:length){
    hit <- sum(labelled_list[1:j, 3])
    #calculate hypergeometric p-value for observing this many or more at this position
    p_val_up <- phyper(hit - 1, n_lab, n_unlab, j, lower.tail=FALSE)
    p_val_down <- phyper(n_lab - hit - 1, n_lab, n_unlab, length - j, lower.tail=FALSE)
    log_p_up <- -log(p_val_up, base=10)
    log_p_down <- -log(p_val_down, base=10)
    if(labelled_list[j,2] >= 0){
      if(log_p_up > max_p){
        max_p <- log_p_up
        max_p_stat <- labelled_list[j,2]
        max_p_g <- hit
        max_p_g_up <- j
      }
    }
    else if(labelled_list[j,2] < 0){
      if(log_p_down > max_m){
        max_m <- log_p_down
        max_m_stat <- labelled_list[j,2]
        max_m_g <- n_lab - hit
        max_m_g_down <- length - j
      }
    }
    vals <- rbind(vals, c(j, log_p_up, log_p_down))
  }
  
  #calculate correlation with phenotype
  up_cor <- cor(labelled_list[,2],vals[,2], method="spearman")
  down_cor <- cor(-labelled_list[,2], vals[,3], method="spearman")
  
  pval_hash[[i]] <- vals
  presence[[i]] <- which(labelled_list[,3])
  lines <- rbind(lines, c(as.character(this_path[1,2]),as.character(this_path[1,3]), n_lab ,max_p,  up_cor, max_p_stat, max_p_g,max_p_g_up, max_m, down_cor, max_m_stat,max_m_g,max_m_g_down))
  
}

lines[,4] <- -log(p.adjust(10^(-as.numeric(lines[,4])), method="BH"), base=10)
lines[,9] <- -log(p.adjust(10^(-as.numeric(lines[,9])), method="BH"), base=10)

lines <- lines[order(apply(lines,1, function(row) max(as.numeric(row[4]), as.numeric(row[9]))), decreasing=TRUE), ]
colnames(lines) <- c("pathway","pathname","path_genes", "up_p","up_cor", "up_max_stat","path_genes_up", "genes_up", "down_p","down_cor","down_max_stat","path_genes_down", "genes_down")

write.table(lines, paste(opt$outdir,"/summary.txt",sep=""), quote=FALSE, row.names=FALSE, sep="\t")

plot_enrich <- function(pathway,direction,stats, rugseq){
  path <- paste(opt$outdir,"/",pathway,sep="")
  maxim = max(stats)
  minim = -maxim/20
  dir.create(path, recursive=TRUE)
  pdf(paste(path,"/plot.",direction,".pdf", sep=""))
  mar.default <- c(5,4,4,2) + 0.1
  par(mar = mar.default + c(0, 1, 0, 0))
  plot(1:length(stats),stats, type='l', xlab="Gene index", ylab="-log10(p-value)", lwd=5, col="#5ab4ac", ylim = c(minim,maxim), cex.lab = 2,
       cex.axis = 2, cex.main =2, main = 'Essentiality of conserved genes')
  lines(c(which(stats == max(stats)),which(stats == max(stats))),c(0,maxim-minim), col="#d8b365", lwd=4)
  # x = max(which(genelist[2] > 0))
  x = nrow(genelist) / 2
  lines(c(x,x),c(0,maxim-minim), col="#8c510a", lwd=4) 
  legend(1300,48, c("Median conservation", "Max P-value"), lty=c(1,1,1,1), lwd=c(4,4,4,4),cex=1.5,
         col=c("#8c510a", "#d8b365"), bty="n")
  rug(rugseq, ticksize = 0.08)
  dev.off()
}

apply(lines,1,function(row) {if(as.numeric(row[4]) > opt$sigthresh){stats <- pval_hash[[row[1]]]; plot_enrich(row[1], "up",stats[,2], presence[[row[1]]])}})
apply(lines,1,function(row) {if(as.numeric(row[9]) > opt$sigthresh){stats <- pval_hash[[row[1]]]; plot_enrich(row[1], "down",stats[,3], presence[[row[1]]])}})

conserved_essential <- labelled_list[1:500,1][labelled_list[1:500,3]]
all_essential <- read.table('~/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials/clusters.info', sep='\t', quote = "")
cons_es_tbl <- all_essential[all_essential[,1]%in%conserved_essential,]
write.table(cons_es_tbl,file = '~/EnTrI/results/walking-hypergeometric-test/essentials-in-top-500-conserved.tsv', quote = FALSE, col.names = FALSE, row.names = FALSE, sep = '\t')
