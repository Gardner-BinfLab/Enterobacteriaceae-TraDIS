insertion_positions <- read.table("../results/insertion-indices/insertion-position-bias.out", row.names=1)
essentialitydir = '../results/biases/dbscan/'
list_of_files = list.files(path=essentialitydir, full.names=T, recursive=FALSE)
nms = c()
essentialitydict = c()
for (filename in list_of_files)
{
  essfile = read.table(filename, as.is=TRUE, header=FALSE, sep="\t")
  nms = c(nms,essfile$V1)
  essentialitydict = c(essentialitydict, essfile$V3)
}
names(essentialitydict) <- nms
pdf("../figures/insertion-position-bias.pdf")
mar.default <- c(5,4,4,2) + 0.1
par(mar = mar.default + c(0, 1, 0, 0))
toplot <- colMeans(insertion_positions[,1:ncol(insertion_positions)-1])
cols = c(rep("#8c510a",5), rep("#80cdc1",75), rep("#01665e",20))
midpoints <- barplot(toplot, main="All", xaxt="n", xlab="Position", ylab="Average insertion index", cex.lab = 2, cex.axis = 2, cex.main =2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=2)
box()
legend(41,0.7, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
essential = c()
nonessential = c()
beneficialloss = c()
# for (i in seq(1,nrow(insertion_positions)))
for (item in row.names(insertion_positions))
{
  # if (insertion_positions[i,ncol(insertion_positions)] < 0.2)
  if (item %in% names(essentialitydict) & essentialitydict[item]=='essential')
  {
    essential = rbind(essential, insertion_positions[item,1:ncol(insertion_positions)-1]) 
  }
  # else if (insertion_positions[i,ncol(insertion_positions)] < 2)
  else if (item %in% names(essentialitydict) & essentialitydict[item]=='non-essential')
  {
    nonessential = rbind(nonessential, insertion_positions[item,1:ncol(insertion_positions)-1])
  }
  # else
  else if (item %in% names(essentialitydict) & essentialitydict[item]=='beneficial-loss')
  {
    beneficialloss = rbind(beneficialloss, insertion_positions[item,1:ncol(insertion_positions)-1]) 
  }
}
toplot <- colMeans(essential)
midpoints <- barplot(toplot, main="Essential", xaxt="n", xlab="Position", ylab="Average insertion index", cex.lab = 2, cex.axis = 2, cex.main =2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=2)
box()
legend(41,0.3, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
toplot <- colMeans(nonessential)
midpoints <- barplot(toplot, main="Non essential", xaxt="n", xlab="Position", ylab="Average insertion index", cex.lab = 2, cex.axis = 2, cex.main =2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=2)
box()
legend(41,0.6, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
toplot <- colMeans(beneficialloss)
midpoints <- barplot(toplot, main="Beneficial loss", xaxt="n", xlab="Position", ylab="Average insertion index", cex.lab = 2, cex.axis = 2, cex.main =2, col=cols, border=NA)
axis(1, at=midpoints[c(1,100)], labels=c('5\'','3\''), cex.axis=2)
box()
legend(41,1.5, c("First 5%","internal", "Last 20%"), lty=c(1,1,1),cex=1.5, lwd=c(4,4, 4),col=c(cols[1],cols[6],cols[81]), bg="white")
dev.off()

# for (item in essential)
# {
#   if (mean(as.numeric(item[2:6])) >= (mean(as.numeric(item[7:81]))+0.02))
#   {
#     print(item)
#   }
# }

# Calculate what percentage of essential genes have insertions in each region of the coding sequence
parts = data.frame(rowSums(essential[,1:5]),rowSums(essential[,6:80]),rowSums(essential[,81:100]))>0
print(colSums(parts)*100/nrow(parts))