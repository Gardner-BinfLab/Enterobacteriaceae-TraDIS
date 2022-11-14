library(rJava)
library(tidyverse)
library(eulerr)
library(grid)

tabpath = '../results/giant-tab/giant-tab.tsv'
table = read.csv(tabpath, sep='\t')[,c(140,142,149)]
table$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential =
  table$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential ==100
table$Endosymbiont..conserved = table$Endosymbiont..conserved==100
table$Enterobacteriaceae..Ancestral.essentiality = table$Enterobacteriaceae..Ancestral.essentiality>0
colnames(table) <- c('Core essential', 'Endosymbiont core', 'Ancestral essential')
table <- table[rowSums(table)>0,]
v <- euler(table)
pdf('../figures/venn-ancestral-core-endosymbiont.pdf')
plot(v, fills=c("#d8b365","#af8dc3","#5ab4ac"), legend = list(cex=2),
     quantities= list(cex=2))
dev.off()

table = read.csv(tabpath, sep='\t')[,c(23,140,149)]
table$EcoGene.Essentiality..Escherichia.coli.BW25113 =
  table$EcoGene.Essentiality..Escherichia.coli.BW25113 == 1
table$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential =
  table$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential ==100
table$Enterobacteriaceae..Ancestral.essentiality = table$Enterobacteriaceae..Ancestral.essentiality>0
colnames(table) <- c('EcoGene essential', 'Core essential', 'Ancestral essential')
table <- table[rowSums(table)>0,]
v <- euler(table)
pdf('../figures/venn-ancestral-core-ecogene.pdf')
plot(v, fills=c("rosybrown4", "#d8b365", "#5ab4ac"), legend = list(cex=2),
     quantities= list(cex=2))
dev.off()
