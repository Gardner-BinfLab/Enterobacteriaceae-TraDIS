library(ComplexHeatmap)
library(stringr)

giant.tab <- read.table('../results/giant-tab/giant-tab.tsv', sep='\t', header = TRUE)

########## ancestral essential
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "Enterobacteriaceae..Ancestral.essentiality",
                   "E..coli...Salmonella...Citrobacter..Ancestral.essentiality",
                   "Salmonella...Citrobacter..Ancestral.essentiality",
                   "Klebsiella..Ancestral.essentiality",
                   "E..coli..Ancestral.essentiality",
                   "Salmonella..Ancestral.essentiality",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Essentiality..Escherichia.coli.ST131.EC958",
                   "TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "EcoGene.Essentiality..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Citrobacter.rodentium.ICC168",
                   "TraDIS.Essentiality..Salmonella.Typhi.Ty2",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.A130",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.D23580",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Essentiality..Salmonella.Enteritidis.P125109")]
tab[, grepl("TraDIS.Essentiality..", names(tab))] <-
  (tab[, grepl("TraDIS.Essentiality..", names(tab))] <= 0 &
     !is.na(tab[, grepl("TraDIS.Essentiality..", names(tab))]))+0
tab <- tab[rowSums(tab[,grepl("TraDIS.Essentiality..", names(tab))]) > 0,]
tab <- tab[tab$Enterobacteriaceae..Ancestral.essentiality==1,
              c("Gene_EGGNOG", "Gene_Keio", "KEGG_acc", "KEGG_pathway",
                "Locus..Klebsiella.pneumoniae.Ecl8",
                "Locus..Klebsiella.pneumoniae.RH201207",
                "Locus..Escherichia.coli.ST131.EC958",
                "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                "Locus..Escherichia.coli.BW25113..Keio.",
                "Locus..Escherichia.coli.BW25113",
                "Locus..Citrobacter.rodentium.ICC168",
                "Locus..Salmonella.Typhi.Ty2",
                "Locus..Salmonella.Typhimurium.A130",
                "Locus..Salmonella.Typhimurium.D23580",
                "Locus..Salmonella.Typhimurium.SL3261",
                "Locus..Salmonella.Typhimurium.SL1344",
                "Locus..Salmonella.Enteritidis.P125109")]
write.table(tab,'../results/giant-tab/ancestral-essential.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## core essential
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential")]
tab <- tab[tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential==100,
           c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
             "Locus..Klebsiella.pneumoniae.Ecl8",
             "Locus..Klebsiella.pneumoniae.RH201207",
             "Locus..Escherichia.coli.ST131.EC958",
             "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
             "Locus..Escherichia.coli.BW25113..Keio.",
             "Locus..Escherichia.coli.BW25113",
             "Locus..Citrobacter.rodentium.ICC168",
             "Locus..Salmonella.Typhi.Ty2",
             "Locus..Salmonella.Typhimurium.A130",
             "Locus..Salmonella.Typhimurium.D23580",
             "Locus..Salmonella.Typhimurium.SL3261",
             "Locus..Salmonella.Typhimurium.SL1344",
             "Locus..Salmonella.Enteritidis.P125109")]
write.table(tab,'../results/giant-tab/core-essential.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## endosymbiont universally conserved
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "Endosymbiont..conserved",
                   "Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential",
                   "Enterobacteriaceae..Ancestral.essentiality")]
tab <- tab[tab$Endosymbiont..conserved==100,
           c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
             "Locus..Klebsiella.pneumoniae.Ecl8",
             "Locus..Klebsiella.pneumoniae.RH201207",
             "Locus..Escherichia.coli.ST131.EC958",
             "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
             "Locus..Escherichia.coli.BW25113..Keio.",
             "Locus..Escherichia.coli.BW25113",
             "Locus..Citrobacter.rodentium.ICC168",
             "Locus..Salmonella.Typhi.Ty2",
             "Locus..Salmonella.Typhimurium.A130",
             "Locus..Salmonella.Typhimurium.D23580",
             "Locus..Salmonella.Typhimurium.SL3261",
             "Locus..Salmonella.Typhimurium.SL1344",
             "Locus..Salmonella.Enteritidis.P125109",
             "Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential",
             "Enterobacteriaceae..Ancestral.essentiality")]
tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential <-
  (tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential == 100)+0
write.table(tab,'../results/giant-tab/endosymbiont-conserved.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## universally essential
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "Bacteria..excluding.symbionts...essential")]
tab <- tab[tab$Bacteria..excluding.symbionts...essential==100,
           c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
             "Locus..Klebsiella.pneumoniae.Ecl8",
             "Locus..Klebsiella.pneumoniae.RH201207",
             "Locus..Escherichia.coli.ST131.EC958",
             "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
             "Locus..Escherichia.coli.BW25113..Keio.",
             "Locus..Escherichia.coli.BW25113",
             "Locus..Citrobacter.rodentium.ICC168",
             "Locus..Salmonella.Typhi.Ty2",
             "Locus..Salmonella.Typhimurium.A130",
             "Locus..Salmonella.Typhimurium.D23580",
             "Locus..Salmonella.Typhimurium.SL3261",
             "Locus..Salmonella.Typhimurium.SL1344",
             "Locus..Salmonella.Enteritidis.P125109")]
write.table(tab,'../results/giant-tab/universally-essential.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## ecogene essentials
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "EcoGene.Essentiality..Escherichia.coli.BW25113")]
tab <- tab[tab$EcoGene.Essentiality..Escherichia.coli.BW25113==1,
           c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
             "Locus..Escherichia.coli.BW25113..Keio.")]
ecogene <- read.table("../results/ecogene-k12.txt", header = FALSE)
missed <- ecogene[!ecogene$V1 %in% tab$Locus..Escherichia.coli.BW25113..Keio.,1]
write.table(tab,'../results/giant-tab/ecogene-essential.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## TraDIS uniquely essential
tab = giant.tab[,c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "Enterobacteriaceae..essential",
                   "Bacteria..excluding.symbionts...essential")]
tab <- tab[tab$Bacteria..excluding.symbionts...essential==35 &
             tab$Enterobacteriaceae..essential==100,
           c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
             "Locus..Klebsiella.pneumoniae.Ecl8",
             "Locus..Klebsiella.pneumoniae.RH201207",
             "Locus..Escherichia.coli.ST131.EC958",
             "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
             "Locus..Escherichia.coli.BW25113..Keio.",
             "Locus..Escherichia.coli.BW25113",
             "Locus..Citrobacter.rodentium.ICC168",
             "Locus..Salmonella.Typhi.Ty2",
             "Locus..Salmonella.Typhimurium.A130",
             "Locus..Salmonella.Typhimurium.D23580",
             "Locus..Salmonella.Typhimurium.SL3261",
             "Locus..Salmonella.Typhimurium.SL1344",
             "Locus..Salmonella.Enteritidis.P125109")]
write.table(tab,'../results/giant-tab/TraDIS-Salmonella_Typhi-uniquely_essenial.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## Endosymbiont heatmap
tab <- giant.tab[giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential ==100
                 | giant.tab$Endosymbiont..conserved==100
                 | giant.tab$Enterobacteriaceae..Ancestral.essentiality>0,
                 grepl("Symbiont..", names(giant.tab))
                 | grepl("Gene_Keio", names(giant.tab))]
tab <- tab[order(tab$Gene_Keio),]
tab <- rbind(tab[tab$Gene_Keio!='-',], tab[tab$Gene_Keio=='-',])
tab$Gene_Keio[nrow(tab)] <- 'EC958_1344'
write.table(tab,'../results/giant-tab/core_ancestral_in_symbionts.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)
rownames(tab) <- tab[,1]
tab <- tab[,-1]
tab[tab == "0"] <- "Absent"
tab[tab == "1"] <- "Present"
pdf("../results/giant-tab/core_ancestral_in_symbionts.pdf", height = 65, width = 25)
colors <- structure(c('white', 'black'), names=c('Absent', 'Present'))
Heatmap(as.matrix(tab), cluster_rows = FALSE, cluster_columns = FALSE,
        column_names_max_height=unit(20, "cm"),
        rect_gp = gpar(col = 'gray'), name = 'Presence',
        col=colors)
dev.off()

########## Core essential: E.coli, Salmonella, Klebsiella
tab <- giant.tab[((giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8" <= 0
                  & giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207" <= 0
                  &giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.Ecl8" == 1
                  & giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.RH201207" == 1)
                 | (giant.tab$"TraDIS.Essentiality..Escherichia.coli.ST131.EC958" <= 0
                    & giant.tab$"TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441" <= 0
                    & giant.tab$"EcoGene.Essentiality..Escherichia.coli.BW25113" == 1
                    & giant.tab$"TraDIS.Essentiality..Escherichia.coli.BW25113" <= 0
                    & giant.tab$"TraDIS.Presence..Escherichia.coli.ST131.EC958" == 1
                    & giant.tab$"TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441" == 1
                    & giant.tab$"EcoGene.Presence..Escherichia.coli.BW25113" == 1
                    & giant.tab$"TraDIS.Presence..Escherichia.coli.BW25113" == 1)
                    | (giant.tab$"TraDIS.Essentiality..Salmonella.Typhi.Ty2" <= 0
                       & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.A130" <= 0
                       & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.D23580" <= 0
                       & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL3261" <= 0
                       & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL1344" <= 0
                       & giant.tab$"TraDIS.Essentiality..Salmonella.Enteritidis.P125109" <= 0
                       & giant.tab$"TraDIS.Presence..Salmonella.Typhi.Ty2" == 1
                       & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.A130" == 1
                       & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.D23580" == 1
                       & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL3261" == 1
                       & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL1344" == 1
                       & giant.tab$"TraDIS.Presence..Salmonella.Enteritidis.P125109" == 1))
                 & giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential < 100,
                 grepl("Gene_EGGNOG", names(giant.tab))
                 | grepl("TraDIS.Essentiality..Klebsiella", names(giant.tab))
                 | grepl("TraDIS.Essentiality..Escherichia", names(giant.tab))
                 | grepl("TraDIS.Essentiality..Salmonella", names(giant.tab))
                 | grepl("EcoGene.Essentiality..Escherichia", names(giant.tab))
                 | grepl("Cluster", names(giant.tab))]

presence <- giant.tab[((giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8" <= 0
                   & giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207" <= 0
                   &giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.Ecl8" == 1
                   & giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.RH201207" == 1)
                  | (giant.tab$"TraDIS.Essentiality..Escherichia.coli.ST131.EC958" <= 0
                     & giant.tab$"TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441" <= 0
                     & giant.tab$"EcoGene.Essentiality..Escherichia.coli.BW25113" == 1
                     & giant.tab$"TraDIS.Essentiality..Escherichia.coli.BW25113" <= 0
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.ST131.EC958" == 1
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441" == 1
                     & giant.tab$"EcoGene.Presence..Escherichia.coli.BW25113" == 1
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.BW25113" == 1)
                  | (giant.tab$"TraDIS.Essentiality..Salmonella.Typhi.Ty2" <= 0
                     & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.A130" <= 0
                     & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.D23580" <= 0
                     & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL3261" <= 0
                     & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL1344" <= 0
                     & giant.tab$"TraDIS.Essentiality..Salmonella.Enteritidis.P125109" <= 0
                     & giant.tab$"TraDIS.Presence..Salmonella.Typhi.Ty2" == 1
                     & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.A130" == 1
                     & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.D23580" == 1
                     & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL3261" == 1
                     & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL1344" == 1
                     & giant.tab$"TraDIS.Presence..Salmonella.Enteritidis.P125109" == 1))
                 & giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential < 100,
                 grepl("TraDIS.Presence..Klebsiella", names(giant.tab))
                 | grepl("TraDIS.Presence..Escherichia", names(giant.tab))
                 | grepl("TraDIS.Presence..Salmonella", names(giant.tab))
                 | grepl("EcoGene.Presence..Escherichia", names(giant.tab))]

presence <- presence[order(tab$Gene_EGGNOG),]
tab <- tab[order(tab$Gene_EGGNOG),]
tab$Gene_EGGNOG[tab$Gene_EGGNOG == '-'] <- str_match(tab$Cluster[tab$Gene_EGGNOG == '-'],'BN373_\\d*')
tab <- tab[,-1]
rownm <- tab[,1]
tab <- tab[,-1]
tab <- as.data.frame(sapply(tab, as.numeric))
rownames(tab) <- rownm
temp <- tab$EcoGene.Essentiality..Escherichia.coli.BW25113
tab[is.na(tab)] <- log(2^6.5)
tab[tab == -4.5] <- NA
tab <- log2(exp(tab))
tab[is.na(tab)] <- -6.5
tab$EcoGene.Essentiality..Escherichia.coli.BW25113 <- temp
tab$EcoGene.Essentiality..Escherichia.coli.BW25113[tab$EcoGene.Essentiality..Escherichia.coli.BW25113 == 0] = 6.5
tab$EcoGene.Essentiality..Escherichia.coli.BW25113[tab$EcoGene.Essentiality..Escherichia.coli.BW25113 == 1] = -6.5
tab[presence == 0] <- NA

nHalf <- 50
Min <- min(tab, na.rm=TRUE)
Max <- max(tab, na.rm=TRUE)
Thresh <- 0
eps <- 1e-10
rc1 <- colorRampPalette(colors = c("#8c510a", "#f6e8c3"), space="Lab")(nHalf)
rc2 <- colorRampPalette(colors = c("#c7eae5", "#01665e"), space="Lab")(nHalf)
ramps <- c(rc1,rc2)
rb1 <- seq(Min, Thresh, length.out=nHalf)
rb2 <- seq(Thresh+eps, Max, length.out=nHalf)
rampbreaks <- c(rb1, rb2)
col_fun = circlize::colorRamp2(rampbreaks, ramps)
pdf("../results/giant-tab/core_ecoli_salmonella_klebsiella.pdf", height = 25)
Heatmap(as.matrix(tab), col = col_fun, name = 'log2(ii/essentiality_threshold)',
        cluster_rows = FALSE, cluster_columns = FALSE, gap = unit(5, "mm"),
        na_col = 'black', column_names_max_height=unit(20, "cm"),
        heatmap_legend_param = list(color_bar = "continuous",
                                    at = c(-6.5,-3,0,3,6.5)))
dev.off()

########## Core essential in E.coli and not all
tab <- giant.tab[(giant.tab$"TraDIS.Essentiality..Escherichia.coli.ST131.EC958" <= 0
                     & giant.tab$"TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441" <= 0
                     & giant.tab$"EcoGene.Essentiality..Escherichia.coli.BW25113" == 1
                     & giant.tab$"TraDIS.Essentiality..Escherichia.coli.BW25113" <= 0
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.ST131.EC958" == 1
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441" == 1
                     & giant.tab$"EcoGene.Presence..Escherichia.coli.BW25113" == 1
                     & giant.tab$"TraDIS.Presence..Escherichia.coli.BW25113" == 1)
                 & giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential < 100,
                 c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "Enterobacteriaceae..essential",
                   "Bacteria..excluding.symbionts...essential",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Presence..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Presence..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Essentiality..Escherichia.coli.ST131.EC958",
                   "TraDIS.Presence..Escherichia.coli.ST131.EC958",
                   "TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "EcoGene.Essentiality..Escherichia.coli.BW25113",
                   "EcoGene.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Escherichia.coli.BW25113",
                   "TraDIS.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Citrobacter.rodentium.ICC168",
                   "TraDIS.Presence..Citrobacter.rodentium.ICC168",
                   "TraDIS.Essentiality..Salmonella.Typhi.Ty2",
                   "TraDIS.Presence..Salmonella.Typhi.Ty2",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.A130",
                   "TraDIS.Presence..Salmonella.Typhimurium.A130",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.D23580",
                   "TraDIS.Presence..Salmonella.Typhimurium.D23580",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Essentiality..Salmonella.Enteritidis.P125109",
                   "TraDIS.Presence..Salmonella.Enteritidis.P125109"
                   )]
tab[,grepl("TraDIS.Essentiality", names(tab))] <-
  lapply(tab[,grepl("TraDIS.Essentiality", names(tab))],
         function(x) ifelse(as.numeric(x) == -4.5, -6.5, log2(exp(as.numeric(x)))))
write.table(tab,'../results/giant-tab/TraDIS_exclusive_core_ecoli.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## Core essential in Klebsiella and not all
tab <- giant.tab[(giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8" <= 0
                  & giant.tab$"TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207" <= 0
                  &giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.Ecl8" == 1
                  & giant.tab$"TraDIS.Presence..Klebsiella.pneumoniae.RH201207" == 1)
                 & giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential < 100,
                 c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "Enterobacteriaceae..essential",
                   "Bacteria..excluding.symbionts...essential",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Presence..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Presence..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Essentiality..Escherichia.coli.ST131.EC958",
                   "TraDIS.Presence..Escherichia.coli.ST131.EC958",
                   "TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "EcoGene.Essentiality..Escherichia.coli.BW25113",
                   "EcoGene.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Escherichia.coli.BW25113",
                   "TraDIS.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Citrobacter.rodentium.ICC168",
                   "TraDIS.Presence..Citrobacter.rodentium.ICC168",
                   "TraDIS.Essentiality..Salmonella.Typhi.Ty2",
                   "TraDIS.Presence..Salmonella.Typhi.Ty2",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.A130",
                   "TraDIS.Presence..Salmonella.Typhimurium.A130",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.D23580",
                   "TraDIS.Presence..Salmonella.Typhimurium.D23580",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Essentiality..Salmonella.Enteritidis.P125109",
                   "TraDIS.Presence..Salmonella.Enteritidis.P125109"
                 )]
tab[,grepl("TraDIS.Essentiality", names(tab))] <-
  lapply(tab[,grepl("TraDIS.Essentiality", names(tab))],
         function(x) ifelse(as.numeric(x) == -4.5, -6.5, log2(exp(as.numeric(x)))))
write.table(tab,'../results/giant-tab/TraDIS_exclusive_core_klebsiella.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)

########## Core essential in Salmonella and not all
tab <- giant.tab[(giant.tab$"TraDIS.Essentiality..Salmonella.Typhi.Ty2" <= 0
                  & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.A130" <= 0
                  & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.D23580" <= 0
                  & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL3261" <= 0
                  & giant.tab$"TraDIS.Essentiality..Salmonella.Typhimurium.SL1344" <= 0
                  & giant.tab$"TraDIS.Essentiality..Salmonella.Enteritidis.P125109" <= 0
                  & giant.tab$"TraDIS.Presence..Salmonella.Typhi.Ty2" == 1
                  & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.A130" == 1
                  & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.D23580" == 1
                  & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL3261" == 1
                  & giant.tab$"TraDIS.Presence..Salmonella.Typhimurium.SL1344" == 1
                  & giant.tab$"TraDIS.Presence..Salmonella.Enteritidis.P125109" == 1)
                 & giant.tab$Enterobacteriaceae..excluding.Salmonella.enterica.serovar.Typhi...essential < 100,
                 c("Gene_EGGNOG","Gene_Keio", "KEGG_acc",	"KEGG_pathway",
                   "Locus..Klebsiella.pneumoniae.Ecl8",
                   "Locus..Klebsiella.pneumoniae.RH201207",
                   "Locus..Escherichia.coli.ST131.EC958",
                   "Locus..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "Locus..Escherichia.coli.BW25113..Keio.",
                   "Locus..Escherichia.coli.BW25113",
                   "Locus..Citrobacter.rodentium.ICC168",
                   "Locus..Salmonella.Typhi.Ty2",
                   "Locus..Salmonella.Typhimurium.A130",
                   "Locus..Salmonella.Typhimurium.D23580",
                   "Locus..Salmonella.Typhimurium.SL3261",
                   "Locus..Salmonella.Typhimurium.SL1344",
                   "Locus..Salmonella.Enteritidis.P125109",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "Enterobacteriaceae..essential",
                   "Bacteria..excluding.symbionts...essential",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Presence..Klebsiella.pneumoniae.Ecl8",
                   "TraDIS.Essentiality..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Presence..Klebsiella.pneumoniae.RH201207",
                   "TraDIS.Essentiality..Escherichia.coli.ST131.EC958",
                   "TraDIS.Presence..Escherichia.coli.ST131.EC958",
                   "TraDIS.Essentiality..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "TraDIS.Presence..Escherichia.coli.UPEC.ST131.NCTC13441",
                   "EcoGene.Essentiality..Escherichia.coli.BW25113",
                   "EcoGene.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Escherichia.coli.BW25113",
                   "TraDIS.Presence..Escherichia.coli.BW25113",
                   "TraDIS.Essentiality..Citrobacter.rodentium.ICC168",
                   "TraDIS.Presence..Citrobacter.rodentium.ICC168",
                   "TraDIS.Essentiality..Salmonella.Typhi.Ty2",
                   "TraDIS.Presence..Salmonella.Typhi.Ty2",
                   "DEG..Salmonella.enterica.serovar.Typhi",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.A130",
                   "TraDIS.Presence..Salmonella.Typhimurium.A130",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.D23580",
                   "TraDIS.Presence..Salmonella.Typhimurium.D23580",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL3261",
                   "TraDIS.Essentiality..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Presence..Salmonella.Typhimurium.SL1344",
                   "TraDIS.Essentiality..Salmonella.Enteritidis.P125109",
                   "TraDIS.Presence..Salmonella.Enteritidis.P125109"
                 )]
tab[,grepl("TraDIS.Essentiality", names(tab))] <-
  lapply(tab[,grepl("TraDIS.Essentiality", names(tab))],
         function(x) ifelse(as.numeric(x) == -4.5, -6.5, log2(exp(as.numeric(x)))))
write.table(tab,'../results/giant-tab/TraDIS_exclusive_core_salmonella.tsv', sep = '\t',
            quote = FALSE, row.names = FALSE)
