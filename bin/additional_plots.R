# script that reads the large table of all essential genes and does a few statistics with
# presence/absence of essentiality, etc.

# load libraries
library(tidyverse)
library(stringr)
library(UpSetR)
library(readxl)
library(cowplot)

setwd("./EnTrI")

# read in the table
giant_table <- read.table("./results/giant-tab/giant-tab.tsv", header = T, sep = "\t")

# read in final table from excel
table_final <- read_xlsx("./data/final_table.xlsx")

# read giant tab all
giant_table_all <- read.table("./results/giant-tab/giant-tab-all.tsv", header = T, sep = "\t")

# select only columns 1,2,3, and the ones that contain "TraDIS" and "Keio"
giant_table <- giant_table %>% select(1,2,3, contains("TraDIS"), contains("Keio"))

# make for all rntries that contain a value <=0 "E", and the ones  >0 "N". remember that now, all the columns are character
# vectors and also contain some characters. So do not change the characters, only the numbers. Use sapply,
# check wether the entry can be made numeric (otherwise itd be "NA"), and then check if it is <=0 or not.
giant_table_adj <- giant_table %>% mutate_all(funs(ifelse(sapply(., function(x) !is.na(as.numeric(x))), ifelse(as.numeric(.) <= 0, "E", "N"), .)))

final_tab <- read.table("./data/giant-tab_final.tsv", header = T, sep = "\t")



# select only used genomes from giant_table_all in colnames :
# # names here: ("BN373"=c(), "ERS227112"=c(), "ROD"=c(), "SL1344"=c(), "SL3261"=c(), "STMMW"=c(), "STM"=c(), "SEN"=c(),
# #                 "t"=c(), "EC958"=c(), "NCTC13441"=c(), "BW25113"=c(), "b"=c())
# new names: (c("Klebsiella pneumoniae Ecl8", "Klebsiella pneumoniae RH201207",
#                      "Citrobacter rodentium ICC168", "Salmonella typhimurium SL1344", "Salmonella typhimurium SL3261",
#                      "Salmonella typhimurium D23580", "Salmonella typhimurium A130", "Salmonella Enteritidis P125109",
#                      "Salmonella typhi Ty2", "Escherichia coli ST131 EC958", "Escherichia coli UPEC ST131"
#                      , "Escherichia coli BW25113", "Escherichia coli K-12 MG1655"))

giant_table_all <- giant_table_all %>% select(1,contains("BN373"), contains("ERS227112"),
                                              contains("ROD"), contains("SL1344"),
                                              contains("SL3261"), contains("STMMW"),
                                              contains("STM"), contains("SEN"), "t",
                                              contains("EC958"), contains("NCTC13441"),
                                              contains("BW25113"), "b")

# rename the the genomes (colnames)
colnames(giant_table_all) <- c("Cluster", "Klebsiella pneumoniae Ecl8", "Klebsiella pneumoniae RH201207",
                               "Citrobacter rodentium ICC168", "Salmonella Typhimurium SL1344", "Salmonella Typhimurium SL3261",
                               "Salmonella Typhimurium D23580", "Salmonella Typhimurium A130", "Salmonella Enteritidis P125109",
                               "Salmonella Typhi Ty2", "Escherichia coli ST131 EC958", "Escherichia coli ST131 NCTC 13441",
                               "Escherichia coli BW25113", "Escherichia coli BW25113 (EcoGene)")

# now create an upset plot with giant_table_all, showing the presence/absence of essentiality in the genomes. essential
# genes contain the number 1 and non-essential genes contain the number 0. Each row is one cluster of genes.
# The plot should show the number ess. genes and their overlap between the genomes.

up1 <- upset(giant_table_all, order.by = "freq", sets = rev(c("Citrobacter rodentium ICC168", "Salmonella Typhi Ty2",
                                                   "Salmonella Enteritidis P125109", "Salmonella Typhimurium SL1344",
                                                   "Salmonella Typhimurium SL3261", "Salmonella Typhimurium D23580",
                                                    "Salmonella Typhimurium A130", "Escherichia coli ST131 NCTC 13441",
                                                    "Escherichia coli ST131 EC958", "Escherichia coli BW25113",
                                                    "Escherichia coli BW25113 (EcoGene)", "Klebsiella pneumoniae RH201207",
                                                    "Klebsiella pneumoniae Ecl8")), nintersects = 20, keep.order = T)

# write as svg
svg("./figures/upsetr.svg")
up1
dev.off()




# OK now try to do a permutaton analysis with e. coli. first, make pivot_longer and add strain column
giant_table_all_long <- giant_table_all %>% pivot_longer(cols = -Cluster, names_to = "Strain", values_to = "Essentiality") %>%
  mutate(genus = gsub(".*Salmonella.*", "Salmonella", Strain)) %>%
    mutate(genus = gsub(".*Escherichia.*", "Escherichia", genus)) %>%
    mutate(genus = gsub(".*Klebsiella.*", "Klebsiella", genus)) %>%
    mutate(genus = gsub(".*Citrobacter.*", "Citrobacter", genus))

# for each CLuster, count the number of essential genes in the Escherichia genus. Then perutation test. E.g. how many
# genes are more significantly more frequently essential within the genus than without when e.g. permuting genus
# labels 1000 times or similar

# first, divide the number of essential genes in the Escherichia genus by the total number of strains that have it as
# essential
giant_table_escherichia <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Escherichia"] == "1")/4) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Escherichia"] == "1")/9) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco > ess_non_eco) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist() %>% length()
print(giant_table_escherichia)

# complete_genus_specific
compl_genus_specific_escerichia <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Escherichia"] == "1")/4) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Escherichia"] == "1")/9) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco == 1 & ess_non_eco == 0) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist()

# same for salmonella
giant_table_salmonella <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Salmonella"] == "1")/5) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Salmonella"] == "1")/8) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco > ess_non_eco) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist() %>% length()

print(giant_table_salmonella)

compl_genus_specific_salmonella <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Salmonella"] == "1")/5) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Salmonella"] == "1")/8) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco == 1 & ess_non_eco == 0) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist()

# same for klebsiella
giant_table_klebsiella <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Klebsiella"] == "1")/4) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Klebsiella"] == "1")/9) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco > ess_non_eco) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist() %>% length()

print(giant_table_klebsiella)

compl_genus_specific_klebsiella <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Klebsiella"] == "1")/2) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Klebsiella"] == "1")/11) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco == 1 & ess_non_eco == 0) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist()

# same for citrobacter
giant_table_citrobacter <- giant_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Citrobacter"] == "1")/4) %>%
  # get the number of essential genes in the non-Escherichia genus and divide by total number of non-escherichia strains (9)
    mutate(ess_non_eco = sum(Essentiality[genus != "Citrobacter"] == "1")/9) %>%
  # get the ones in which ess_eco > ess_non_eco
  filter(ess_eco > ess_non_eco) %>%
  # get only the names of the clusters.
    select(Cluster) %>%
    # remove duplicates
    unique() %>%
    # make a vector
    unlist() %>% length()

print(giant_table_citrobacter)



# ok now import the tabel connecting the clusters to the genes
rep_genes_clusters <- read.table("./results/all-hieranoid/cluster_info.tsv", header = F, sep = "\t", col.names = "gene" ) %>%
  # add enumerated clusters per row (Clust1 for row 1 and so on)
    mutate(Cluster = paste("Clust", 0:(nrow(.)-1), sep = ""))

rgene_raw <- rep_genes_clusters

# add clusters containing all genes to rep_genes_clusters to rep_genes_clusters
# read clusters
all_clusts <- read.table("./results/all-hieranoid/clusters.txt", header = F, sep = ",") %>%
  mutate(clusters = gsub("\t", ",", V1)) %>% select(clusters) %>% unlist() %>% as.character()

rep_genes_clusters$whole_cluster <- sapply(paste0(rep_genes_clusters$gene, ","), function(x) {
  all_clusts[grepl(x, all_clusts)]
  #all_clusts[grepl(x, all_clusts)]
})

rep_genes_clusters$whole_cluster <- sapply(rep_genes_clusters$whole_cluster,function (x) x[[1]])

# remove duplicate rows that have the same cluster, and remove the second occuring row
rep_genes_clusters <- rep_genes_clusters[!duplicated(rep_genes_clusters$whole_cluster),]


# get all locus tags hidden in the all_clusts data. Remove the number with which they end. (any length of numbers before
# first non number)
# split them by , to extract all lts
lts <- unlist(strsplit(all_clusts, ",")) %>% unlist() %>% gsub("[0-9]*$", "", .) %>% unique()
lts <- gsub("_.*", "", lts) %>% unique()
lts <- lts[!grepl("DEG", lts)]
lts <- gsub("^SEN.*", "SEN", lts) %>% unique()
lts <- lts[!grepl("ordLocTag|AE0LocTag", lts)]
lts <- lts[lts!="ENC"]

# get mappings
names_genomes  <- gsub("Locus: ", "", colnames(table_final[,18:64]))
# add locus tag shortform for them :  [1] "Klebsiella pneumoniae Ecl8"                                              "Klebsiella pneumoniae RH201207"
#  [3] "Escherichia coli ST131 EC958"                                            "Escherichia coli UPEC ST131 NCTC13441"
#  [5] "Escherichia coli BW25113 (Keio)"                                         "Escherichia coli BW25113"
#  [7] "Citrobacter rodentium ICC168"                                            "Salmonella Typhi Ty2"
#  [9] "Salmonella Typhimurium A130"                                             "Salmonella Typhimurium D23580"
# [11] "Salmonella Typhimurium SL3261"                                           "Salmonella Typhimurium SL1344"
# [13] "Salmonella Enteritidis P125109"                                          "Secondary endosymbiont of Ctenarytaina eucalypti"
# [15] "Candidatus Moranella endobia PCIT"                                       "Secondary endosymbiont of Heteropsylla cubana"
# [17] "Candidatus Baumannia cicadellinicola strain BGSS"                        "Candidatus Baumannia cicadellinicola strain B-GSS"
# [19] "Baumannia cicadellinicola str. Hc (Homalodisca coagulata)"               "Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757"
# [21] "Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675"  "Candidatus Blochmannia vafer str. BVAF"
# [23] "Candidatus Blochmannia floridanus"                                       "Candidatus Blochmannia pennsylvanicus str. BPEN"
# [25] "Candidatus Blochmannia chromaiodes str. 640"                             "Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis"
# [27] "Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans" "Buchnera aphidicola str. Sg (Schizaphis graminum)"
# [29] "Buchnera aphidicola str. G002 (Myzus persicae)"                          "Buchnera aphidicola str. USDA (Myzus persicae)"
# [31] "Buchnera aphidicola str. F009 (Myzus persicae)"                          "Buchnera aphidicola str. W106 (Myzus persicae)"
# [33] "Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)"                      "Buchnera aphidicola str. APS (Acyrthosiphon pisum)"
# [35] "Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)"                     "Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)"
# [37] "Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)"                    "Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)"
# [39] "Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)"                     "Buchnera aphidicola str. 5A (Acyrthosiphon pisum)"
# [41] "Buchnera aphidicola str. Ua (Uroleucon ambrosiae)"                       "Buchnera aphidicola str. Bp (Baizongia pistaciae)"
# [43] "Buchnera aphidicola BCc"                                                 "Buchnera aphidicola (Cinara tujafilina)"
# [45] "Sodalis glossinidius str. morsitans"                                     "Sodalis praecaptivus strain HS1"
# [47] "Candidatus Sodalis pierantonius str. SOPE"

lt_genomes <- c("BN373", "ERS227112", "EC958", "NCTC13441", "b", "BW25113", "ROD", "t", "STM", "STMMW", "SL3261",
                "SL1344", "SEN", "A359", "MEPCIT", "A35E", "IM45", "AB162", "BCI", "BOBLI757", "BTURN675", "BVAF",
                "Bfl", "BPEN", "BCHRO640", "BA000021", "WIGMOR", "BUsg", "BUMPG002", "BUMPUSDA", "BUMPF009", "BUMPW106",
                "BAKON", "BA000003", "BUAPTUC7", "CWO", "CWQ", "CWU", "CWS", "BUAP5A", "BUAMB", "bbp", "BCc", "BCTU",
                "SG", "Sant", "SOPEG")
names(names_genomes) <- lt_genomes

# create new, columns for all the lt and extract the respective lt from the whole_cluster column
for (lt in lts)
{
  # get genome name
  genome <- names_genomes[lt]
  str_ltag <- paste0("Locus: ", genome)
  # from the whole_cluster column, extract the lt if it exists, else NA
  re <- paste0(".*(^|,)(", lt, "[^,]+),.*")
  added_col <-  sapply(rep_genes_clusters$whole_cluster, function(x) {
    tags <- ifelse(grepl(re, x), gsub(re, "\\2", x), NA)
    })

    # add the new column to the rep_genes_clusters df
  rep_genes_clusters[[str_ltag]]  <- added_col
}

# manually adjust the one for STM
lt <- "STM"
genome <- names_genomes[lt]
str_ltag <- paste0("Locus: ", genome)
str_TraDIS <- paste0("TraDIS Essentiality: ", genome)
re <- paste0(".*(^|,)(", lt, "_[^,]+),.*")
added_col <-  sapply(rep_genes_clusters$whole_cluster, function(x) {
  tags <- ifelse(grepl(re, x), gsub(re, "\\2", x), NA)
})
rep_genes_clusters[[str_ltag]]  <- added_col

# adjust keio bw (ecogene)
rep_genes_clusters$`Locus: Escherichia coli BW25113 (Keio)` <- ifelse(grepl("^b[0-9]+$", rep_genes_clusters$gene),
                                                                      rep_genes_clusters$gene, NA)

# find locuses that contain STMM
rep_genes_clusters$`Locus: Salmonella Typhimurium A130`[grepl("STMM", rep_genes_clusters$`Locus: Salmonella Typhimurium A130`)]
# from Locus: Salmonella Typhimurium A130 replace lts that contain STMM with NA
rep_genes_clusters$`Locus: Salmonella Typhimurium A130` <- gsub("STMMw.*", NA, rep_genes_clusters$`Locus: Salmonella Typhimurium A130`)


# ok, now we need to get the tradis data, saved in txt files for each genome in a folder.
for (genome in names(names_genomes)[1:13][!names(names_genomes)[1:13]=="b"])
{
  # read in the tradis data
  tradis <- read.table(paste0("./results/biases/dbscan/", genome, ".txt"), header = F, sep = "\t")
  # rename cols
  colnames(tradis) <- c("locus_tag", "Insertions", "Essentiality", "score_loge")
  # from the loge score , convert it to log2
  tradis$score_log2 <- log2(exp(tradis$score_loge))
  # add the log2 score to a new column in the rep_genes_clusters df, starting with "TraDIS: "
  genome_full_name <- names_genomes[genome]
  rep_genes_clusters[[paste0("TraDIS Essentiality: ", genome_full_name)]] <- ifelse(rep_genes_clusters[[paste0("Locus: ", genome_full_name)]] %in% tradis$locus_tag,
                                                                         tradis$score_log2[match(rep_genes_clusters[[paste0("Locus: ", genome_full_name)]], tradis$locus_tag)], NA)
  # print # of genes < = 0 in the resp. genome, and the total number of genes in the genome
  # use the curated rep_genes_clusters for this
  tot_genes <- sum(na.omit(rep_genes_clusters[[paste0("Locus: ", genome_full_name)]]) != "NA")
  ess_genes <- sum(na.omit(rep_genes_clusters[[paste0("TraDIS Essentiality: ", genome_full_name)]]) <= 0)
  print(names_genomes[genome])
  print(c(ess_genes, tot_genes))
}

sum(na.omit(rep_genes_clusters$`TraDIS Essentiality: Klebsiella pneumoniae Ecl8`)<0)
# remove row containing gene ROD__00001 in rep_genes_clusters$gene
rep_genes_clusters <- rep_genes_clusters[rep_genes_clusters$gene != "ROD__00001",]

# Ok now get the locus tags only of the endosymbionts
names_endosymbionts <- paste0 ("Locus: ",names_genomes[14:47])


rep_genes_clusters_endosymbionts <- rep_genes_clusters[names_endosymbionts]
# find out how many rows have 0 NAs:
sum(rowSums(is.na(rep_genes_clusters_endosymbionts)) == 0)

# add it as a column to the table:
rep_genes_clusters["Conserved in Symbiont Genomes"] <- ifelse(rowSums(is.na(rep_genes_clusters_endosymbionts)) == 0, 1, 0)
# replace nas in rep_genes_clusters["Conserved in Symbiont Genomes"] with 0
rep_genes_clusters["Conserved in Symbiont Genomes"] <- ifelse(is.na(rep_genes_clusters[["Conserved in Symbiont Genomes"]]), 0, rep_genes_clusters[["Conserved in Symbiont Genomes"]])



# add column "EcoGene Essentiality: Escherichia coli BW25113". take values from final_table
rep_genes_clusters["EcoGene Essentiality: Escherichia coli BW25113"] <- ifelse(!is.na(rep_genes_clusters$`Locus: Escherichia coli BW25113`),
                                                                            table_final$`EcoGene Essentiality: Escherichia coli BW25113`[match(rep_genes_clusters$`Locus: Escherichia coli BW25113`, table_final$`Locus: Escherichia coli BW25113`)], 0)
# replace NAs with 0
rep_genes_clusters["EcoGene Essentiality: Escherichia coli BW25113"] <- ifelse(is.na(rep_genes_clusters[[ "EcoGene Essentiality: Escherichia coli BW25113"]]), 0, rep_genes_clusters[[ "EcoGene Essentiality: Escherichia coli BW25113"]])
# add column of "Keio gene name" in rep_genes_clusters in the same way
rep_genes_clusters["Keio gene name"] <- ifelse(!is.na(rep_genes_clusters$`Locus: Escherichia coli BW25113`),
                                               table_final$`Keio gene name`[match(rep_genes_clusters$`Locus: Escherichia coli BW25113`, table_final$`Locus: Escherichia coli BW25113`)], NA)

# add column of Core Essentiality in rep_genes_clusters in the same way
rep_genes_clusters["Core Essentiality"] <- ifelse(!is.na(rep_genes_clusters$`Locus: Escherichia coli BW25113`),
                                                  table_final$`Core Essentiality`[match(rep_genes_clusters$`Locus: Escherichia coli BW25113`, table_final$`Locus: Escherichia coli BW25113`)], 0)
# replace NAs with 0
rep_genes_clusters["Core Essentiality"] <- ifelse(is.na(rep_genes_clusters$`Core Essentiality`), 0, rep_genes_clusters$`Core Essentiality`)

# add column of " Ancestral Essentiality" in rep_genes_clusters in the same way
rep_genes_clusters["Ancestral Essentiality"] <- ifelse(!is.na(rep_genes_clusters$`Locus: Escherichia coli BW25113`),
                                                      table_final$`Ancestral Essentiality`[match(rep_genes_clusters$`Locus: Escherichia coli BW25113`, table_final$`Locus: Escherichia coli BW25113`)], 0)
# replace NAs with 0
rep_genes_clusters["Ancestral Essentiality"] <- ifelse(is.na(rep_genes_clusters$`Ancestral Essentiality`), 0, rep_genes_clusters$`Ancestral Essentiality`)

# add 1 ancestral gene in row containing SL1344_1967 in "Locus: Salmonella Typhimurium SL1344"
rep_genes_clusters["Ancestral Essentiality"] <- ifelse(grepl("SL1344_1967", rep_genes_clusters$`Locus: Salmonella Typhimurium SL1344`), 1, rep_genes_clusters$'Ancestral Essentiality')



# remove rows that have NA in all columns from column 4 to 62
rep_genes_clusters <- rep_genes_clusters[rowSums(is.na(rep_genes_clusters[,4:62])) != 59,]


# print("get all genes for PastML - gene (Cluster) as column and genome as rows")
# tab_all_genes_pastml <- rep_genes_clusters[,grepl("Locus:|Cluster", colnames(rep_genes_clusters))][,c(1,6:18)]
# # make all NAs a 0 and all non-NAs a 1 in all columns except the first
# tab_all_genes_pastml[,2:14] <- ifelse(is.na(tab_all_genes_pastml[,2:14]), 0, 1)
# # make the first column the rownames
# rownames(tab_all_genes_pastml) <- tab_all_genes_pastml$Cluster
# # remove the first column
# tab_all_genes_pastml <- tab_all_genes_pastml[,-1]
# # transpose
# tab_all_genes_pastml <- t(tab_all_genes_pastml)
# # rename the rownames by the genome names
# rownames(tab_all_genes_pastml) <- gsub("Locus: ", "", rownames(tab_all_genes_pastml))
# # make a df from the matrix
# tab_all_genes_pastml <- as.data.frame(tab_all_genes_pastml)
# # add rownames as column
# tab_all_genes_pastml <- tab_all_genes_pastml %>% rownames_to_column("Genome")
# # get number of clusters that are fully conserved in all genomes
# fully_conserved_clusters <- sum(colSums(tab_all_genes_pastml[-1]) == 13)
# # its 2455!
# # remove the clusters that are fully conserved in all genomes
# tab_all_genes_pastml <- tab_all_genes_pastml[,c(T,colSums(tab_all_genes_pastml[-1])!=13)]
# # check the nr of genes that are not present in any of the genomes
# genes_not_present <- sum(colSums(tab_all_genes_pastml[-1]) == 0)
# # its 3476. remove them as well
# tab_all_genes_pastml <- tab_all_genes_pastml[,c(T,colSums(tab_all_genes_pastml[-1])!=0)]
#
# # ok now save the table as tsv
# write.table(tab_all_genes_pastml, "./data/pastml_reconstruction/table_all_genes.tsv", sep = "\t", quote = F,
#             row.names = F, col.names = T)

# ok now import the table with the essentiality of the genes in the genomes
pastml_result_all <- read.table("./data/pastml_reconstruction/pastml_all_genomes_all_genes.txt", header = T, sep = "\t")

# get rowsums of pastml result, excluding the first column and na.omit
pastml_result_all_rowsums <- rowSums(pastml_result_all[,-1], na.rm = T)
names(pastml_result_all_rowsums) <- pastml_result_all$node
# get # of fully conserved genes from
fully_conserved_clusters
# get all colnames that have a 1 in the first row.
clusters_conserved_all <- colnames(pastml_result_all)[pastml_result_all[1,] == 1]
pastml_result_all_rowsums + fully_conserved_clusters



# add gene to giant_table_all

giant_table_all <- giant_table_all %>% left_join(rgene_raw, by = "Cluster") %>%
  # make gene the second column
    select(Cluster, gene, everything())


tab_only_ess <- giant_table_all[rowSums(giant_table_all[,-c(1,2)]) >= 0,]
# make first column rownames
rownames(tab_only_ess) <- tab_only_ess$gene
# remove first column
tab_only_ess <- tab_only_ess[,-c(1,2)]
# transpose
tab_only_ess <- t(tab_only_ess)
# rename rownames bythe initial codes
rownames(tab_only_ess) <- c("Klebsiella pneumoniae Ecl8", "Klebsiella pneumoniae RH201207", "Citrobacter rodentium ICC168",
                            "Salmonella typhimurium SL1344", "Salmonella typhimurium SL3261",
                            "Salmonella Typhimurium D23580", "Salmonella Typhimurium A130",
                            "Salmonella Enteritidis P125109", "Salmonella Typhi Ty2", "Escherichia coli ST131 EC958",
                            "Escherichia coli UPEC ST131", "Escherichia coli BW25113", "Escherichia coli K-12 MG1655")
# make df from matrix
tab_only_ess <- as.data.frame(tab_only_ess)

# add rownames as column
tab_only_ess <- tab_only_ess %>% rownames_to_column("Genome")

# get nr of genes fully conserved essentiality (colsums 13)
fully_cons_essentiality <- sum(colSums(tab_only_ess[-1]) == 13)

genes_fully_ess <- colnames(tab_only_ess)[-1][colSums(tab_only_ess[-1]) == 13]

# get # of genes that occur in all of our strains:
genes_fully_conserved <- giant_table_all$gene[rowSums(giant_table_all[,-c(1,2)]) == 13]

# OK now we got the genes that are ancestrally esse

# save as tsv. take care that the rownames are not saved as a column
write.table(tab_only_ess[,c(T,colSums(tab_only_ess[-1])!=13)], "./data/pastml_reconstruction/table_only_ess_genes.tsv", sep = "\t", quote = F,
            row.names = F, col.names = T)





# read pastml result
pastml_result <- read.table("./data/pastml_reconstruction/pastml_ess_genes_nonconserved.txt", header = T, sep = "\t")

#get rowsums of pastml result, excluding the first column and na.omit
pastml_result_rowsums <- rowSums(pastml_result[,-1], na.rm = T)
names(pastml_result_rowsums) <- pastml_result$node
pastml_result_rowsums + fully_cons_essentiality

# get all colnames that have a 1 in the first row.
clusters_conserved <- colnames(pastml_result)[pastml_result[1,] == 1]

# add pastml ancestral genes to rep_genes_clusters
rep_genes_clusters["PastML Ancestral Essentiality"] <- ifelse((rep_genes_clusters$Cluster %in% clusters_conserved)|
                                                                (rep_genes_clusters$`Core Essentiality`==1), 1, 0)

# find keio gene names not ancestral in pastml but in other ancestral
keio_not_anc <- rep_genes_clusters$`Keio gene name`[rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`PastML Ancestral Essentiality`==0]

# get pastml numbers of all genes (pastml_all_genomes_nonconserved.txt)
pastml_all <- read.table("./data/pastml_reconstruction/pastml_all_genomes_nonconserved.txt", header = T, sep = "\t")

# get rowsums of pastml result, excluding the first column and na.omit
pastml_all_rowsums <- rowSums(pastml_all[,-1], na.rm = T)
names(pastml_all_rowsums) <- pastml_all$node
# get # of fully conserved genes from
pastml_all_rowsums + sum(fully_cons_essentiality)

# get all colnames that have a 1 in the first row.
clusters_conserved_all <- colnames(pastml_all)[pastml_all[1,] == 1]

# add pastml ancestral genes to rep_genes_clusters
rep_genes_clusters["PastML Ancestral genes"] <- ifelse(rep_genes_clusters$Cluster %in% clusters_conserved_all, 1, 0)


# save as xlsx and make title row bold
library(xlsx)
saved_xlsx <- rep_genes_clusters[,c("Cluster", "PastML Ancestral Essentiality" , colnames(table_final))]

# add endosymbiont
write.table(saved_xlsx, "./data/final_table_JJ.tsv", row.names = F, sep = "\t", quote = F)


# map these clusters back to the genes
genes_conserved <- rgene_raw$gene[rgene_raw$Cluster %in% clusters_conserved]

# get ancestral essential genes (genes_conserved and genes_fully_ess)
anc_ess_genes <- c(genes_conserved, genes_fully_ess)

# ok now get the ecogene essential genes:
ecogene_ess_genes <- read.table("./results/ecogene-k12.txt", header = F, sep = "\t")$V1

# okay now create an euler diagram with the three sets of genes
library(eulerr)
# get number of genes in each set
n_all_three <- sum(genes_fully_ess %in% ecogene_ess_genes & genes_fully_ess %in% anc_ess_genes)
n_ecogene_anc <- sum(genes_fully_ess %in% ecogene_ess_genes & !(genes_fully_ess %in% anc_ess_genes))
n_ecogene <- sum(!(genes_fully_ess %in% ecogene_ess_genes) & genes_fully_ess %in% anc_ess_genes)
n_anc <- sum(!(genes_fully_ess %in% ecogene_ess_genes) & !(genes_fully_ess %in% anc_ess_genes))
n_ecogene_only <- sum(! (ecogene_ess_genes %in% anc_ess_genes))
n_anc_only <- sum(! (anc_ess_genes %in% ecogene_ess_genes))
n_anc_and_ecogene_but_not_fully <- sum(anc_ess_genes %in% ecogene_ess_genes) - n_all_three

library(viridis)
cols <- viridis(8, option = "viridis")

# create euler diagram
fit <- euler(c("Ancestral" = n_anc_only, "Ecogene" = n_ecogene_only, "Ancestral&Ecogene" = n_anc_and_ecogene_but_not_fully,
               "Ancestral&Core&Ecogene" = n_all_three,  "Core" = n_ecogene))
plot(fit, quantities = TRUE, labels = c("A", "E",  "ACE", "AE"),
     fills = list(fill = c("steelblue", cols[8], cols[1], cols[6]), alpha = 0.7))

svg("./figures/euler_venn_entero.svg", width = 5, height = 5)
plot(fit, quantities = TRUE, labels = c("A", "E",  "ACE", "AE"),
     fills = list(fill = c("steelblue", cols[8], cols[1], cols[6]), alpha = 0.7))
dev.off()


# do same with ancestral, core and endosymbionts. get endosymbionts genes:
endosymbiont_ess_genes <- table_final$`Locus: Escherichia coli BW25113 (Keio)`[table_final$`Conserved in Symbiont Genomes`==1]

# get number of genes in each set
n_all_three <- sum(genes_fully_ess %in% endosymbiont_ess_genes & genes_fully_ess %in% anc_ess_genes)
n_endosymbiont_anc <- sum(endosymbiont_ess_genes %in% anc_ess_genes & !(endosymbiont_ess_genes %in% genes_fully_ess))
n_endosymbiont <- sum(!(endosymbiont_ess_genes %in% genes_fully_ess) & !(endosymbiont_ess_genes %in% anc_ess_genes))
n_anc <- sum(!(anc_ess_genes %in% endosymbiont_ess_genes) & !(anc_ess_genes %in% genes_fully_ess))
n_core <- sum(!(genes_fully_ess %in% endosymbiont_ess_genes) & !(genes_fully_ess %in% intersect(anc_ess_genes, endosymbiont_ess_genes)))

# create euler diagram
fit <- euler(c("Ancestral" = n_anc,  "Endosymbiont" = n_endosymbiont, "Core&Ancestral" = n_core,
               "Ancestral&Core&Endosymbiont" = n_all_three, "Ancestral&Endosymbiont" = n_endosymbiont_anc))
plot(fit, quantities = TRUE, labels = NULL,
        fills = list(fill = c("steelblue", "orange2", cols[1]), alpha = 0.7))

svg("./figures/euler_venn_endosymbiont.svg", width = 5, height = 5)
plot(fit, quantities = TRUE, labels = NULL,
     fills = list(fill = c("steelblue", "orange2", cols[1]), alpha = 0.7))
dev.off()

# do same for old (fitchs) ancestral genes
n_all_three <- sum(rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==1 &
                     rep_genes_clusters$`EcoGene Essentiality: Escherichia coli BW25113`==1)
n_ecogene_anc <- sum(rep_genes_clusters$`EcoGene Essentiality: Escherichia coli BW25113`==1 &
                       rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==0)
n_ecogene <- sum(rep_genes_clusters$`EcoGene Essentiality: Escherichia coli BW25113`==1 &
                      rep_genes_clusters$`Ancestral Essentiality`==0 & rep_genes_clusters$`Core Essentiality`==0)
n_anc <- sum(rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==0 &
                  rep_genes_clusters$`EcoGene Essentiality: Escherichia coli BW25113`==0)

fit <- euler(c("Ancestral" = n_anc, "Ecogene" = n_ecogene, "Ancestral&Ecogene" = n_ecogene_anc,
               "Ancestral&Core&Ecogene" = n_all_three))
plot(fit, quantities = TRUE, labels = c("A", "E",  "AE", "ACE"),
        fills = list(fill = c("steelblue", cols[8], cols[1], cols[6]), alpha = 0.7))

svg("./figures/euler_venn_old.svg", width = 5, height = 5)
plot(fit, quantities = TRUE, labels = c("A", "E",  "AE", "ACE"),
     fills = list(fill = c("steelblue", cols[8], cols[1], cols[6]), alpha = 0.7))
dev.off()

# same for old euler with endosymbionts
n_all_three <- sum(rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==1 &
                     rep_genes_clusters$`Conserved in Symbiont Genomes`==1)
n_endosymbiont_anc <- sum(rep_genes_clusters$`Conserved in Symbiont Genomes`==1 &
                       rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==0)
n_endosymbiont <- sum(rep_genes_clusters$`Conserved in Symbiont Genomes`==1 &
                        rep_genes_clusters$`Ancestral Essentiality`==0 & rep_genes_clusters$`Core Essentiality`==0)
n_anc <- sum(rep_genes_clusters$`Ancestral Essentiality`==1 & rep_genes_clusters$`Core Essentiality`==0 &
                    rep_genes_clusters$`Conserved in Symbiont Genomes`==0)
n_core_anc <- sum(rep_genes_clusters$`Core Essentiality`==1 & rep_genes_clusters$`Ancestral Essentiality`==1 &
                     rep_genes_clusters$`Conserved in Symbiont Genomes`==0)

fit <- euler(c("Ancestral" = n_anc, "Endosymbiont" = n_endosymbiont, "Core&Ancestral" = n_core_anc,
                "Ancestral&Core&Endosymbiont" = n_all_three, "Ancestral&Endosymbiont" = n_endosymbiont_anc))
plot(fit, quantities = TRUE, labels = NULL,
        fills = list(fill = c("steelblue", "orange2", cols[1]), alpha = 0.7))

svg("./figures/euler_venn_old_endosymbiont.svg", width = 5, height = 5)
plot(fit, quantities = TRUE, labels = NULL,
     fills = list(fill = c("steelblue", "orange2", cols[1]), alpha = 0.7))
dev.off()






# get old table of ancestreal essentiality!
old_anc <- read.csv("./results/core-ancestral/ancestral.csv", header = T, sep = ",")

sum(old_anc$gene_tag %in% genes_conserved)
# just 2 are not ther. now check which genes are neither in the conserved nor in the fully conserved genes
old_anc$gene_tag[!old_anc$gene_tag %in% genes_conserved & !old_anc$gene_tag %in% genes_fully_ess]


# from table_final, only select TRADIS and Keio columns, and add an index column as first column 1
# Warning: `funs()` was deprecated in dplyr 0.8.0.
# ℹ Please use a list of either functions or lambdas
table_final_filtered <- table_final %>% select(contains("TraDIS"), contains("Locus: Escherichia coli BW25113 (Keio)")) %>%
  # add 0 if NA and 1 if not NA. do not use funs function, rather lambda or function
    mutate_all(function(x) ifelse(x=="NA", 0, 1)) %>%
    # add index column
    mutate(index = 1:nrow(.)) %>%
    # make index column first column
    select(index, everything())

# get table_final_filtered but make NA string an actual NA and make numeric columns
table_final_filtered_no_na <- table_final %>% select(contains("TraDIS")) %>%
    # add 0 if NA and 1 if not NA. do not use funs function, rather lambda or function
        mutate_all(function(x) ifelse(x=="NA", NA, as.numeric(x))) %>%
        # add index column
        mutate(index = 1:nrow(.)) %>%
        # make index column first column
        select(index, everything())




#  change colnames to the genome ids. order is different now:
#  [1] "index"
#  [2] "TraDIS Essentiality: Klebsiella pneumoniae Ecl8"
#  [3] "TraDIS Essentiality: Klebsiella pneumoniae RH201207"
#  [4] "TraDIS Essentiality: Escherichia coli ST131 EC958"
#  [5] "TraDIS Essentiality: Escherichia coli UPEC ST131 NCTC13441"
#  [6] "TraDIS Essentiality: Escherichia coli BW25113"
#  [7] "TraDIS Essentiality: Citrobacter rodentium ICC168"
#  [8] "TraDIS Essentiality: Salmonella Typhi Ty2"
#  [9] "TraDIS Essentiality: Salmonella Typhimurium A130"
# [10] "TraDIS Essentiality: Salmonella Typhimurium D23580"
# [11] "TraDIS Essentiality: Salmonella Typhimurium SL3261"
# [12] "TraDIS Essentiality: Salmonella Typhimurium SL1344"
# [13] "TraDIS Essentiality: Salmonella Enteritidis P125109"
# [14] "Locus: Escherichia coli BW25113 (Keio)"
# but change to ids instead of names
colnames(table_final_filtered) <- c("index","BN373", "ERS227112", "EC958", "NCTC13441", "BW25113", "ROD", "t", "SEN",
                                    "STM", "SL3261", "SL1344", "STMMW", "b")

# transpose the table, with genome and indexes as colnames
table_final_filtered <- table_final_filtered %>% select(-index) %>% t() %>% as.data.frame() %>% rownames_to_column("Genome")

# get only colss which only contain 1s
table_final_filtered_conserved <- table_final_filtered %>% select_if(function(x) all(x==1))

# get only colss which not only contain 1s
table_final_filtered_not_conserved <- table_final_filtered %>% select_if(function(x) !all(x==1))

# get # of genes that are conserved fully (i.e. only ones in column)


# save as tsv. take care that the rownames are not saved as a column
write.table(table_final_filtered, "./data/pastml_reconstruction/all_genes_final_filtered.tsv", sep = "\t", quote = F,
            row.names = F, col.names = T)


# go through the table ./resullts/all_hieranoid/clusters.txt . lines have tab separated vaues with different numbers of
# columns. go through the lines, and save all locus tags ()
# go through lines in file:
cluster_list <- list()
i <- 1
for (line in readLines("./results/all-hieranoid/clusters.txt"))
{
  cluster_name <- paste0("cluster_", i)
  # split line by tab to get all genes in cluster
    genes <- str_split(line, "\t")[[1]]
  # remove ones that start with "exDEG"
    genes <- genes[!grepl("^exDEG", genes)]
  # extract '([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+' from each gene to get the locus tag
    genes <- gsub("([a-zA-Z0-9]+)_[a-zA-Z0-9]+", "\\1", genes)
  # include all SEN genes
    genes <- gsub(".*SEN.*", "SEN", genes)
    genes <- gsub("(DEG\\d\\d\\d\\d).+", "\\1", genes)
  # make t from genes starting with t\\d
    genes <- gsub("(t|b)\\d+", "\\1", genes)
  # make SG from genes starting with SG\\d
    genes <- gsub("(SG|SEN|Bfl)\\d+", "\\1", genes)
  # remove _\\d
    genes <- gsub("_\\d+", "", genes)
    # remove _
    genes <- gsub("_", "", genes)
  # remove empty strings
    genes <- genes[genes != ""]
  print(genes)
  cluster_list[[cluster_name]] <- genes
  i <- i + 1
}

# now create a table with the cluster names as rownames and the locus tags as columns. 1 if the locus tag is in the
# cluster, 0 if not
# first, get all unique locus tags
locus_tags <- unique(unlist(cluster_list))
locus_tags

# make a matrix with the locus tags as rownames and the cluster names as colnames
cluster_table <- matrix(0, nrow = length(locus_tags), ncol = length(cluster_list))
rownames(cluster_table) <- locus_tags
colnames(cluster_table) <- names(cluster_list)

# now go through the cluster_list and add 1 to the matrix if the locus tag is in the cluster
for (cluster in names(cluster_list))
{
  for (locus_tag in cluster_list[[cluster]])
  {
    cluster_table[locus_tag, cluster] <- 1
  }
}

# select only the rows from the genome ids previously selected
cluster_table <- cluster_table[rownames(cluster_table) %in% table_final_filtered$Genome,]

# make a df from the matrix and add the locus tags as first column
cluster_table <- as.data.frame(cluster_table)
cluster_table <- cluster_table %>% rownames_to_column("Genome")

x <- rowSums(cluster_table[-1])
names(x) <- cluster_table$Genome
x

# rename the rownames to full names (Salmonella enteritidis P125109, Salmonella typhimurium SL1344, instead of
#  [1] "SEN"       "SL1344"    "SL3261"    "STMMW"     "STM"       "t"
#  [7] "ROD"       "BW25113"   "b"         "NCTC13441" "EC958"     "ERS227112"
# [13] "BN373"
cluster_table$Genome <- c("Salmonella enteritidis P125109", "Salmonella typhimurium SL1344", "Salmonella typhimurium SL3261",
                          "Salmonella typhimurium D23580", "Salmonella typhimurium A130", "Salmonella typhi Ty2",
                        "Citrobacter rodentium ICC168", "Escherichia coli BW25113", "Escherichia coli K-12 MG1655",
                          "Escherichia coli UPEC ST131", "Escherichia coli ST131 EC958", "Klebsiella pneumoniae RH201207",
                            "Klebsiella pneumoniae Ecl8")

# get nr of clusters that only have 1s (13 colsums)
nr_conserved_genes <- sum(colSums(cluster_table[,-1]) == 13)

# tak these from the df
cluster_table_not_conserved <- cluster_table[,c(T,colSums(cluster_table[,-1]) != 13)]
# also take out the ones which have colsums of 0
cluster_table_not_conserved <- cluster_table_not_conserved[,c(T,colSums(cluster_table_not_conserved[,-1]) != 0)]

#save as tsv. take care that the rownames are not saved as a column
write.table(cluster_table_not_conserved, "./data/pastml_reconstruction/all_genomes_nonconserved.tsv", sep = "\t", quote = F,
            row.names = F, col.names = T)

# now run pastml.

# read in results of pastml:
pastml_result_tot <- read.table("./data/pastml_reconstruction/pastml_all_genomes_nonconserved.txt", header = T, sep = "\t")

#get rowsums of pastml result, excluding the first column and na.omit
pastml_result_tot_rowsums <- rowSums(pastml_result_tot[,-1], na.rm = T)
names(pastml_result_tot_rowsums) <- pastml_result_tot$node
pastml_result_tot_rowsums + nr_conserved_genes

# # create tibble with pastml_result_tot_rowsums,pastml_result_rowsums and the genome names as columns
# pastml_results <- tibble(genome = names(pastml_result_tot_rowsums),
#                          total_genes = pastml_result_tot_rowsums + nr_conserved_genes,
#                          ess_genes = pastml_result_rowsums+ fully_cons_essentiality)
#
# df_total <- data.frame(value = pastml_result_tot_rowsums + nr_conserved_genes, node = names(pastml_result_tot_rowsums) )
# df_essgenes <- data.frame(value = pastml_result_rowsums+ fully_cons_essentiality, node = names(pastml_result_rowsums) )
#
# View(df_total)
# View(df_essgenes)

rpoe <- c("rpoE", "degS", "rseP")
# do same for genes that are almost genus specific (trpS,ipA,rsgA,ftsX,fepG,wzyE,ubiH,ftsE,fepC,ubiF)
alm_genus_spec <- c("trpS","fepC","fepG", "ftsX","ftsE","rsgA", "ubiF","lipA", "ubiH", "wzyE")

# get a heatmap, showing all genes that are ancestrally essential but not core essential, and are essential in ecogene.
heatmap_table <- rep_genes_clusters %>%
  filter((`Ancestral Essentiality` == 1 & `Core Essentiality` == 0 & `EcoGene Essentiality: Escherichia coli BW25113` == 1)| (`Keio gene name` %in% rpoe) | (`Keio gene name` %in% alm_genus_spec)) %>%
  mutate("TraDIS Essentiality: Escherichia coli BW25113 (EcoGene)" = ifelse(`EcoGene Essentiality: Escherichia coli BW25113` == 1, -6.5, 6.5)) %>%
    #select(colnames(rep_genes_clusters)[grepl("TraDIS Ess|Keio gene name", colnames(rep_genes_clusters))]) %>%
  select(c(paste0("TraDIS Essentiality: ", c("Citrobacter rodentium ICC168", "Salmonella Typhi Ty2",
                                                   "Salmonella Enteritidis P125109", "Salmonella Typhimurium SL1344",
                                                   "Salmonella Typhimurium SL3261", "Salmonella Typhimurium D23580",
                                                    "Salmonella Typhimurium A130", "Escherichia coli UPEC ST131 NCTC13441",
                                                    "Escherichia coli ST131 EC958", "Escherichia coli BW25113",
                                                    "Escherichia coli BW25113 (EcoGene)",
                                             "Klebsiella pneumoniae RH201207",
                                                    "Klebsiella pneumoniae Ecl8")), "Keio gene name")) %>%
  # select only genes with clear evidence for non-essentiality (essentiality score > 1) in at least two strains
  filter(rowSums(.[,1:14] > 1) >= 2 | `Keio gene name` %in% c(rpoe, alm_genus_spec)) %>%
  column_to_rownames("Keio gene name") %>%
  t() %>% as.data.frame() %>%
  # change row names gsub
  rownames_to_column("Genome") %>% mutate(Genome = gsub("TraDIS Essentiality: ", "", Genome)) %>%
  mutate(Genome = gsub("Klebsiella", "K.", Genome)) %>%
    mutate(Genome = gsub("Escherichia coli", "E. coli", Genome)) %>%
  mutate(Genome = gsub("Salmonella", "S.", Genome)) %>%
    mutate(Genome = gsub("Citrobacter", "C.", Genome)) %>%
  mutate(Genome = gsub("UPEC ", "", Genome)) %>%
  column_to_rownames("Genome")

# make rpoe columns last 3
heatmap_table <- cbind(heatmap_table[,!colnames(heatmap_table) %in% c(rpoe, alm_genus_spec)],
                        heatmap_table[,colnames(heatmap_table) %in% rpoe], heatmap_table[,colnames(heatmap_table) %in% alm_genus_spec])

# same with pastml results
heatmap_table_pastml <- rep_genes_clusters %>%
  filter((`PastML Ancestral Essentiality` == 1 & `Core Essentiality` == 0 & `EcoGene Essentiality: Escherichia coli BW25113` == 1)| (`Keio gene name` %in% rpoe) | (`Keio gene name` %in% alm_genus_spec)) %>%
    mutate("TraDIS Essentiality: Escherichia coli BW25113 (EcoGene)" = ifelse(`EcoGene Essentiality: Escherichia coli BW25113` == 1, -6.5, 6.5)) %>%
  select(c(paste0("TraDIS Essentiality: ", c("Citrobacter rodentium ICC168", "Salmonella Typhi Ty2",
                                                   "Salmonella Enteritidis P125109", "Salmonella Typhimurium SL1344",
                                                   "Salmonella Typhimurium SL3261", "Salmonella Typhimurium D23580",
                                                    "Salmonella Typhimurium A130", "Escherichia coli UPEC ST131 NCTC13441",
                                                    "Escherichia coli ST131 EC958", "Escherichia coli BW25113",
                                                    "Escherichia coli BW25113 (EcoGene)",
                                             "Klebsiella pneumoniae RH201207",
                                                    "Klebsiella pneumoniae Ecl8")), "Keio gene name")) %>%
    # select only genes with clear evidence for non-essentiality (essentiality score > 1) in at least two strains
    filter(rowSums(.[,1:14] > 1) >= 2 | `Keio gene name` %in% c(rpoe, alm_genus_spec)) %>%
    column_to_rownames("Keio gene name") %>%
    t() %>% as.data.frame() %>%
    # change row names gsub
    rownames_to_column("Genome") %>% mutate(Genome = gsub("TraDIS Essentiality: ", "", Genome)) %>%
    mutate(Genome = gsub("Klebsiella", "K.", Genome)) %>%
    mutate(Genome = gsub("Escherichia coli", "E. coli", Genome)) %>%
    mutate(Genome = gsub("Salmonella", "S.", Genome)) %>%
    mutate(Genome = gsub("Citrobacter", "C.", Genome)) %>%
    mutate(Genome = gsub("UPEC ", "", Genome)) %>%
    column_to_rownames("Genome")




# make rpoe columns last 3
heatmap_table_pastml <- cbind(heatmap_table_pastml[,!colnames(heatmap_table_pastml) %in% c(rpoe, alm_genus_spec)],
                        heatmap_table_pastml[,colnames(heatmap_table_pastml) %in% rpoe], heatmap_table_pastml[,colnames(heatmap_table_pastml) %in% alm_genus_spec])





library(ComplexHeatmap)
library(circlize)
# create complex heatmaps

# add coloring
col_fun <- colorRamp2(c(-6.5, -5,0, 5, 6.5), c("steelblue", "steelblue", "white", "darkorange", "darkorange"))

ht_anc_fitch <- Heatmap(heatmap_table, cluster_rows = F, cluster_columns = T,
               col = col_fun,
               show_heatmap_legend = F,
               column_split = factor(c(rep("Ancestrally essential but not\ncore essential genes", 24), rep("σE\ngenes", 3),
                                rep("Genus specific\ngenes", 10)), levels = c("Ancestrally essential but not\ncore essential genes",
                                "σE\ngenes", "Genus specific\ngenes")),

               row_title_side = "right", row_title_rot = 0,
                        rect_gp = gpar(col = "white", lwd = 1.5),
               border = T, border_gp = gpar(col = "black", lwd = 1.5),
                        column_gap = unit(0.3, "cm"),
               column_names_gp = gpar(fontsize = 11, fontface = "bold.italic"),
               column_title_gp = gpar(fontsize = 12),
               row_names_gp = gpar(fontsize = 10),
               row_title = NULL,
               width = unit(20, "cm"), height = unit(7, "cm"))

ht_anc_fitch

lgd <- Legend(col_fun = col_fun, title = expression("Essentiality score"), #direction = "horizontal",
             title_gp = gpar(fontsize = 15), labels = c("", "-5", " 0"," 5", ""), legend_height = unit(7, "cm"),
              grid_width = unit(0.5, "cm"),
              labels_gp = gpar(fontsize = 12),
             at = c(-6.5,-5, 0, 5,6.5), border = "black",
             title_position = "leftcenter-rot")

svg("./figures/heatmap_anc_fitch_new.svg", width = 12, height = 5)
ht_anc_fitch
draw(lgd, x = unit(2, "cm"), y = unit(5.92, "cm"))
dev.off()

# do same with pastml results
ht_anc_pastml <- Heatmap(heatmap_table_pastml, cluster_rows = F, cluster_columns = T,
               col = col_fun,
               show_heatmap_legend = F,
               column_split = factor(c(rep("Ancestrally essential but not\ncore essential genes", 16), rep("σE\ngenes", 3),
                                rep("Genus specific\ngenes", 10)), levels = c("Ancestrally essential but not\ncore essential genes",
                                "σE\ngenes", "Genus specific\ngenes")),

               row_title_side = "right", row_title_rot = 0,
                        rect_gp = gpar(col = "white", lwd = 1.5),
               border = T, border_gp = gpar(col = "black", lwd = 1.5),
                        column_gap = unit(0.3, "cm"),
               column_names_gp = gpar(fontsize = 11, fontface = "bold.italic"),
               column_title_gp = gpar(fontsize = 12),
               row_names_gp = gpar(fontsize = 10),
               row_title = NULL,
               width = unit(17, "cm"), height = unit(7, "cm"))

ht_anc_pastml

svg("./figures/heatmap_anc_pastml_update.svg", width = 12, height = 5)
ht_anc_pastml
draw(lgd, x = unit(3, "cm"), y = unit(5.92, "cm"))
dev.off()







