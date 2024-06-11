# This script will do some additional statistics taken from the output of the main analysis.
# It wil run hypergeometric tests

library(dplyr)
library(tidyr)
library(KEGGREST)
library(xlsx)
library(tibble)
library(readxl)


setwd("./EnTrI")

# Load the data
final_table <- read_xlsx("./data/final_table_JJ.xlsx", sheet = 1)
old_final_table <- read_xlsx("./data/final_table.xlsx", sheet = 1)

# make all values of TraDIS -6.49212768400034 to -6.5
final_table_s <- final_table %>% mutate(across(starts_with("TraDIS"), ~ifelse(. == -6.49212768400034, -6.5, .)))
# save as xlsx with first row as column names
write.table(final_table_s, file = "./data/final_table_JJ_new.tsv", sep = "\t", row.names = FALSE)




# compare different ancestral essentiality:
g_not_in_pastml <- final_table$`Keio gene name`[final_table$`Ancestral Essentiality` == 1 & final_table$`PastML Ancestral Essentiality`==0]
g_only_in_pastml <- final_table$`Keio gene name`[final_table$`Ancestral Essentiality` == 0 & final_table$`PastML Ancestral Essentiality`==1]
# print without "" and combine into one string
paste(g_not_in_pastml, collapse = ", ")

# make an euler plot showing all the different ancestrality genes
library(eulerr)
n_only_pastml <- length(g_only_in_pastml)
n_only_ancestral <- length(g_not_in_pastml)
n_both <- length(final_table$`Keio gene name`[final_table$`Ancestral Essentiality` == 1 & final_table$`PastML Ancestral Essentiality`==1])
# create a tibble with the data
fit <- euler(c("PastML" = n_only_pastml, "Fitch" = n_only_ancestral, "PastML&Fitch" = n_both))
# plot the euler plot
svg("./figures/euler_plot_pastml_fitch.svg")
plot(fit, quantities = TRUE, fills = c("yellow", "steelblue", "lightyellow"))
dev.off()


# extract the tradis data and the locus tags
tradis_data <- final_table %>% select("Cluster","Keio gene name",colnames(final_table)[grepl("TraDIS", colnames(final_table))])
old_tradis_data <- old_final_table %>% select("Keio gene name",colnames(old_final_table)[grepl("TraDIS", colnames(old_final_table))])

tradis_data_bin <- final_table %>% filter(`EcoGene Essentiality: Escherichia coli BW25113` == 1 & `Ancestral Essentiality` == 1 & `Core Essentiality` == 0) %>%
    select("Cluster","Keio gene name",colnames(final_table)[grepl("TraDIS", colnames(final_table))])

old_tradis_data_bin <- old_final_table %>% filter(`EcoGene Essentiality: Escherichia coli BW25113` == 1 & `Ancestral Essentiality` == 1 & `Core Essentiality` == 0) %>%
    select("Keio gene name",colnames(old_final_table)[grepl("TraDIS", colnames(old_final_table))])

# get tradis data with clear evidence for non-essentiality (essentiality score > 1) in at least two strains
tradis_data_bin %>% filter(rowSums(across(starts_with("TraDIS"), ~. > 1)) >= 1) %>% View()
# get 24 genes from table:
twenty_four <- tradis_data_bin %>% filter(rowSums(across(starts_with("TraDIS"), ~.=="NA")) == 0) %>%
   # filter data with tradis values >= 2
     filter(rowSums(across(starts_with("TraDIS"), ~. > 1)) > 0) %>% select("Keio gene name") %>% unlist() %>%
  as.character()

na_gnes <- tradis_data_bin %>% filter(rowSums(across(starts_with("TraDIS"), ~.=="NA")) >0 ) %>%
  select("Keio gene name") %>% unlist() %>% as.character()


# filter tradis data with na values in tradis data
tradis_data_28 <- tradis_data_bin %>%
  # filter data with tradis values >= 2
    filter(!rowSums(across(starts_with("TraDIS"), ~. >= 1)) > 0) %>%
  # filter data with >2 tradis values >=0
    filter(!rowSums(across(starts_with("TraDIS"), ~. >= 0)) > 2) %>% select("Keio gene name") %>% unlist() %>% as.character()
# get nr of nas in old_tradis_data_bin
nr_nas <- old_tradis_data_bin %>% filter(rowSums(across(starts_with("TraDIS"), ~.=="NA")) > 0) %>% nrow()

old_tradis_data_28 <- old_tradis_data_bin %>%
  # make all tradis cols number
    mutate(across(starts_with("TraDIS"), as.numeric)) %>%
  # filter data with tradis values >= 2
    filter(!rowSums(across(starts_with("TraDIS"), ~. >= 1)) > 0) %>%
  # filter data with >2 tradis values >=0
    filter(!rowSums(across(starts_with("TraDIS"), ~. >= 0)) > 2) %>% select("Keio gene name") %>% unlist() %>% as.character()

tradis_data_28[!tradis_data_28 %in% old_tradis_data_28]
old_tradis_data_28[!old_tradis_data_28 %in% tradis_data_28]

# find genes that are essential in all strains except 1 with a value of betwen 0 and +1
# get the number of strains
tradis_data_bin_27 <- tradis_data_bin %>%
  mutate(across(starts_with("TraDIS"), ~ifelse(. >= 0 & . < 1, 1, ifelse(. < 0 , -1, NA)))) %>%
  # remove all rows containing NAs
    #filter_all(all_vars(!is.na(.))) %>%
  # remove all rows with only -1s across tradis data
  #filter(rowSums(across(starts_with("TraDIS"), ~. == -1)) < 12) %>%
  # remove all with more than 1 1
    filter(rowSums(across(starts_with("TraDIS"), ~. == 1)) < 3)

old_tradis_data_bin_27 <- old_tradis_data_bin %>% mutate(across(starts_with("TraDIS"), as.numeric)) %>%
  mutate(across(starts_with("TraDIS"), ~ifelse(. >= 0 & . < 1.1, 1, ifelse(. < 0 , -1, NA)))) %>%
  # remove all rows containing NAs
  #filter_all(all_vars(!is.na(.))) %>%
  # remove all rows with only -1s across tradis data
  #filter(rowSums(across(starts_with("TraDIS"), ~. == -1)) < 12) %>%
  # remove all with more than 1 1
  filter(rowSums(across(starts_with("TraDIS"), ~. == 1)) < 3)

dim(tradis_data_bin_27)
dim(old_tradis_data_bin_27)


poss_ast <- tradis_data_bin %>% select("Keio gene name") %>% unique() %>% unlist() %>% as.character()
# find wheher pos_ast gesen are essential in ecogene
poss_ast_ecogene <- final_table$`Keio gene name`[final_table$`Keio gene name` %in% poss_ast & final_table$`EcoGene Essentiality: Escherichia coli BW25113` == 1]

# get # of genes with at least one strain has a value of > 1 in the tradis data
tradis_24 <- tradis_data_bin %>% filter(rowSums(across(starts_with("TraDIS"), ~. > 1)) > 0) %>% select("Keio gene name") %>% unique() %>% unlist() %>% as.character()

tradis_data_bin$`Keio gene name`[!tradis_data_bin$`Keio gene name` %in% tradis_24 & !tradis_data_bin$`Keio gene name` %in% tradis_data_28]

# get locus tag data
locus_tags_data <- final_table %>% select("Cluster","Ancestral Essentiality",
                                          "EcoGene Essentiality: Escherichia coli BW25113", "Core Essentiality",
                                          "Conserved in Symbiont Genomes" ,
                                          colnames(final_table)[grepl("Locus:", colnames(final_table))]) %>% select(1:18)

# filter out all rows with the string "NA" in any of the columns
locus_tags_data <- locus_tags_data %>% filter_all(all_vars(!grepl("NA", .)))

# ok now load the keggrest data
keggnames <- keggList("pathway", "eco")
keggnames <- gsub(" - Escherichia coli K-12 MG1655", "", keggnames)
link_kegg <- keggLink("pathway", "eco")
names(link_kegg) <- gsub("eco:", "", names(link_kegg))
link_kegg <- gsub("path:", "", link_kegg)
# create a tibble with the kegg data
kegg_data <- tibble(locus_tag = names(link_kegg), kegg_name = link_kegg,
                    pathway = keggnames[link_kegg])

locus_tags_all <- locus_tags_data$`Locus: Escherichia coli BW25113 (Keio)`
locus_tags_ancestral <- locus_tags_data$`Locus: Escherichia coli BW25113 (Keio)`[locus_tags_data$`Ancestral Essentiality` == 1]
locus_tags_core <- locus_tags_data$`Locus: Escherichia coli BW25113 (Keio)`[locus_tags_data$`Core Essentiality` == 1]
locus_tags_EcoGene <- locus_tags_data$`Locus: Escherichia coli BW25113 (Keio)`[locus_tags_data$`EcoGene Essentiality: Escherichia coli BW25113` == 1]
locus_tags_Endosymb <- locus_tags_data$`Locus: Escherichia coli BW25113 (Keio)`[locus_tags_data$`Conserved in Symbiont Genomes` == 1]

# go through kegg pathways and do hypergeometric tests
# for each pathway, get the number of genes in the pathway
# and the number of genes in the pathway that are Ancestrally essential.
N <- length(locus_tags_all)
s_a <- length(locus_tags_ancestral)
s_c <- length(locus_tags_core)
s_e <- length(locus_tags_EcoGene)
s_endo <- length(locus_tags_Endosymb)


# go through pathways and chaeck for hypergeometric p-values
pathway_pvalues_a <- c()
pathway_pvalues_c <- c()
pathway_pvalues_e <- c()
pathway_pvalues_endo <- c()

# go through pathways and chaeck for hypergeometric p-values. We go through each pathway and do following steps:
# 1. get the genes in the pathway
# 2. get the genes in the pathway that are Ancestrally essential
# 3. get the genes in the pathway that are Core essential
# 4. get the genes in the pathway that are EcoGene essential
# 5. get the genes in the pathway that are conserved in endosymbionts
# 6. calculate the hypergeometric p-value for each of the four groups of genes, by comparing the number of genes in the pathway
# that are essential in the group to randomly drawn genes from the group, considering the total number of genes in the genome,
# and the number of genes in the pathway and the number of essential genes in the group.
for (pathway_id in names(keggnames)) {
    pathway_name <- keggnames[pathway_id]
    genes_in_pathway <- kegg_data %>% filter(pathway == pathway_name) %>% pull(locus_tag)
    genes_in_pathway_ancestral <- genes_in_pathway [genes_in_pathway %in% locus_tags_ancestral]
    genes_in_pathway_core <- genes_in_pathway [genes_in_pathway %in% locus_tags_core]
    genes_in_pathway_EcoGene <- genes_in_pathway [genes_in_pathway %in% locus_tags_EcoGene]
    genes_in_pathway_Endosymb <- genes_in_pathway [genes_in_pathway %in% locus_tags_Endosymb]

    M <- length(genes_in_pathway)

    k_a <- length(genes_in_pathway_ancestral)
    k_c <- length(genes_in_pathway_core)
    k_e <- length(genes_in_pathway_EcoGene)
    k_endo <- length(genes_in_pathway_Endosymb)

    print(paste0("Pathway: ", pathway_name,", N:", N, " M: ", M, " s: ", s_a, " k: ", k_a))

    pvalue_a <- phyper(k_a-1, s_a, N-s_a, M, lower.tail = FALSE)
    pvalue_c <- phyper(k_c-1, s_c, N-s_c, M, lower.tail = FALSE)
    pvalue_e <- phyper(k_e-1, s_e, N-s_e, M, lower.tail = FALSE)
    pvalue_endo <- phyper(k_endo-1, s_endo, N-s_endo, M, lower.tail = FALSE)


    pathway_pvalues_a[pathway_name] <- pvalue_a #* length(keggnames)
    pathway_pvalues_c[pathway_name] <- pvalue_c #* length(keggnames)
    pathway_pvalues_e[pathway_name] <- pvalue_e #* length(keggnames)
    pathway_pvalues_endo[pathway_name] <- pvalue_endo # length(keggnames)
    print(pvalue_a)
    print(pvalue_c)
    print(pvalue_e)
    print(pvalue_endo)
    # print their gene names
    #print(final_table$`Keio gene name`[final_table$`Locus: Escherichia coli BW25113 (Keio)` %in% genes_in_pathway])
}
# do fdr correction for the p-values
pathway_pvalues_a <- p.adjust(pathway_pvalues_a, method = "fdr")
pathway_pvalues_c <- p.adjust(pathway_pvalues_c, method = "fdr")
pathway_pvalues_e <- p.adjust(pathway_pvalues_e, method = "fdr")
pathway_pvalues_endo <- p.adjust(pathway_pvalues_endo, method = "fdr")

names(pathway_pvalues_a)[pathway_pvalues_a < 0.01][order(pathway_pvalues_a[pathway_pvalues_a < 0.01])]
names(pathway_pvalues_c)[pathway_pvalues_c < 0.01]
names(pathway_pvalues_e)[pathway_pvalues_e < 0.01]
names(pathway_pvalues_endo)[pathway_pvalues_endo < 0.01]


# create a plot of the p-values (barplot with -log10(p-value) on the x axis). use ggplot2. show only pathways with p-value <0.05
pathway_pvalues <- tibble(pathway = factor(names(pathway_pvalues_a), levels = names(pathway_pvalues_a)[order(pathway_pvalues_a, decreasing = TRUE)]),
                                                        Ancestral = pathway_pvalues_a,
                          Core = pathway_pvalues_c, EcoGene = pathway_pvalues_e, Endosymbiont = pathway_pvalues_endo) %>%
  # preserve only rows with p-value < 0.05 in one of the three columns
    filter(Ancestral < 0.01 | Core < 0.01 | EcoGene < 0.01 | Endosymbiont < 0.01) %>%
  # make longer (expand)
    pivot_longer(cols = c("Ancestral", "Core", "EcoGene", "Endosymbiont"), names_to = "Essentiality", values_to = "pathway_pvalues") %>%
  # make 1 out of pvalues that are >1
    mutate(pathway_pvalues = ifelse(pathway_pvalues > 1, 0.99, pathway_pvalues)) %>%
  # make Essentiality a factor and order it
    mutate(Essentiality = factor(Essentiality, levels = c("Endosymbiont", "EcoGene", "Core", "Ancestral")))



library(ggplot2)
library(viridis)

keggplot_all <- ggplot(pathway_pvalues, aes(x = pathway, y = -log10(pathway_pvalues), fill = Essentiality)) +
  # add geom bar but make bars next to one another
  geom_bar(stat = "identity", position = "dodge")  +
  scale_fill_manual(values = c( "darkorange", viridis(100)[99],viridis(100)[3], "steelblue")) +
  ylab("-log10 P-value (adjusted)") +
  xlab("") +
  theme_classic() +
  # transpose x and y axis
    coord_flip() +
  # add line to show significance threshold
  geom_hline(yintercept = -log10(0.01), linetype = "dotted", color = "black", size = 0.5) +
  theme(axis.text.x = element_text(size = 10, family="Sans", color = "black"),
        axis.text.y = element_text(size = 10, family="Sans", color = "black"),
        axis.title = element_text(size = 10, family="Sans"),
        # add legend text
        legend.text = element_text(size = 8, family="Sans"),
        # remove legend title
        legend.title = element_blank(),
        # make legend within plot, not at bottom  (position 0.9)
        legend.position = c(0.8, 0.5))

keggplot_all

ggsave("./figures/keggplot_all.svg", keggplot_all, width = 6, height = 4)


# get allgenes conserved in endosymbionts but not ancestrally essential in tradis data
# get endosymbiont essential genes
endosymb_ess_not_anc <- final_table$`Locus: Escherichia coli BW25113 (Keio)`[final_table$`Ancestral Essentiality` == 0 & final_table$`Conserved in Symbiont Genomes` == 1][order(final_table$`Keio gene name`[final_table$`Ancestral Essentiality` == 0 & final_table$`Conserved in Symbiont Genomes` == 1])]

# do the same kegg analysis for these genes
N <- length(locus_tags_all)
s <- length(endosymb_ess_not_anc)

# go through pathways and chaeck for hypergeometric p-values
pathway_endo_pvalues <- c()
for (pathway_id in names(keggnames)) {
  pathway_name <- keggnames[pathway_id]
  genes_in_pathway <- kegg_data %>% filter(pathway == pathway_name) %>% pull(locus_tag)
  genes_in_pathway_endo <- genes_in_pathway [genes_in_pathway %in% endosymb_ess_not_anc]
  M <- length(genes_in_pathway)
  k <- length(genes_in_pathway_endo)
  print(paste0("Pathway: ", pathway_name,", N:", N, " M: ", M, " s: ", s, " k: ", k))
  pvalue <- phyper(k-1, s, N-s, M, lower.tail = FALSE)
  pathway_endo_pvalues[pathway_name] <- pvalue *length(keggnames)
  print(pvalue)
  # print their gene names
  print(final_table$`Keio gene name`[final_table$`Locus: Escherichia coli BW25113 (Keio)` %in% genes_in_pathway])
}

pathway_endo_pvalues[pathway_endo_pvalues < 0.05]



# now let's find genus-specific essential genes with a permutation test
# get gene-specific ess. genes:
tradis_table_all_long <- tradis_data %>% select(-"Keio gene name") %>%
  pivot_longer(cols = -Cluster, names_to = "Strain", values_to = "Essentiality") %>%
  mutate(genus = gsub(".*Salmonella.*", "Salmonella", Strain)) %>%
    mutate(genus = gsub(".*Escherichia.*", "Escherichia", genus)) %>%
    mutate(genus = gsub(".*Klebsiella.*", "Klebsiella", genus)) %>%
    mutate(genus = gsub(".*Citrobacter.*", "Citrobacter", genus)) %>%
  mutate(Strain = gsub("TraDIS Essentiality: ", "", Strain)) %>%
  mutate(Essentiality = ifelse(Essentiality <= 0 , 1, 0))

# get the genus-specific essential genes
genus_specific_escherichia <- tradis_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
    mutate(ess_eco = sum(Essentiality[genus == "Escherichia"] == 1)) %>%
  mutate(ess_others = sum(Essentiality[genus != "Escherichia"] == 1)) %>%
  # one row per cluster
    distinct(Cluster, .keep_all = TRUE) %>%
  # get the ones that have 4 essential genes in escherichia and 0 in the others
    filter(ess_eco == 4 & ess_others == 0)

# do same for salmonella
genus_specific_salmonella <- tradis_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
  mutate(ess_sal = sum(Essentiality[genus == "Salmonella"] == 1)) %>%
  mutate(ess_others = sum(Essentiality[genus != "Salmonella"] == 1)) %>%
  # one row per cluster
  distinct(Cluster, .keep_all = TRUE) %>%
  # get the ones that have 4 essential genes in escherichia and 0 in the others
  filter(ess_sal == 5 & ess_others == 0)

# do same for klebsiella
genus_specific_klebsiella <- tradis_table_all_long %>% group_by(Cluster) %>%
  # get the number of essential genes in the Escherichia genus and divide by total number of escherichia strains (4)
  mutate(ess_kleb = sum(Essentiality[genus == "Klebsiella"] == 1)) %>%
  mutate(ess_others = sum(Essentiality[genus != "Klebsiella"] == 1)) %>%
  # one row per cluster
  distinct(Cluster, .keep_all = TRUE) %>%
  # get the ones that have 4 essential genes in escherichia and 0 in the others
  filter(ess_kleb == 2 & ess_others == 0)


genus_spec <- genus_specific_escherichia %>% filter(ess_eco == 3) %>% select(Cluster)


tradis_table_all_long <- tradis_data %>%
  # remove rows with NAs
  filter_all(all_vars(!grepl("NA", .))) %>%
  select(-"Keio gene name") %>%
  # select only ones with no NAs in the TraDIS data
  pivot_longer(cols = -Cluster, names_to = "Strain", values_to = "Essentiality") %>%
  mutate(genus = gsub(".*Salmonella.*", "Salmonella", Strain)) %>%
    mutate(genus = gsub(".*Escherichia.*", "Escherichia", genus)) %>%
    mutate(genus = gsub(".*Klebsiella.*", "Klebsiella", genus)) %>%
    mutate(genus = gsub(".*Citrobacter.*", "Citrobacter", genus)) %>%
  mutate(Strain = gsub("TraDIS Essentiality: ", "", Strain)) %>%
  mutate(Essentiality = ifelse(Essentiality <= 0 , 1, 0)) %>%
  # remove clusters with 0 essential genes in all strains
  group_by(Cluster) %>%
    filter(sum(Essentiality) > 0) %>%
  # add keio gene names again
    left_join(final_table[,c("Cluster","Keio gene name")], by = "Cluster") %>%
    ungroup()



# now walk through each cluster, shuffle the genus column 10 times and count the number of times there is evidence for
# genus-specificity for the gene.

e_coli_genus_spec_sign <- c()
salm_genus_spec_sign <- c()

for (cluster in unique(tradis_table_all_long$Cluster)) {
  cluster_data <- tradis_table_all_long %>% filter(Cluster == cluster)
  # add 1 row with cols "Clust3", "Escherichia coli K12 Keio", [essentiality], "Escherichia", [gene]
  # get keio essentiality
  keio_ess <- final_table[final_table$Cluster == cluster,
                            c("Keio gene name", "EcoGene Essentiality: Escherichia coli BW25113")]
  # add row
  cluster_data <- rbind(cluster_data, c(cluster, "Escherichia coli K12 Keio",
                                        keio_ess$`EcoGene Essentiality: Escherichia coli BW25113`, "Escherichia",
                                        keio_ess$`Keio gene name`))
  # check whether cluster is genus specific:
  cluster_data$Essentiality <- as.numeric(cluster_data$Essentiality)
  is_genus_specific_escher <- sum(cluster_data$genus == "Escherichia" & cluster_data$Essentiality == 1) /
    sum(cluster_data$genus == "Escherichia") - sum(cluster_data$genus != "Escherichia" & cluster_data$Essentiality == 1)/
    sum(cluster_data$genus != "Escherichia")

  genename <- keio_ess$`Keio gene name`[1]

  if (is_genus_specific_escher > 0) {
    print(paste0("Escherichia coli: ",genename))
      genus_specificity <- c()
      for (i in 1:10000) {
        cluster_data <- cluster_data %>% mutate(genus = sample(genus))
        is_genus_specific_p <- sum(cluster_data$genus == "Escherichia" & cluster_data$Essentiality == 1) /
          sum(cluster_data$genus == "Escherichia") - sum(cluster_data$genus != "Escherichia" & cluster_data$Essentiality == 1)/
          sum(cluster_data$genus != "Escherichia")
        genus_specificity <- c(genus_specificity, is_genus_specific_p)
      }
    # plot the distribution of genus-specificity
    pval <- sum(genus_specificity >= is_genus_specific_escher) / length(genus_specificity)
    print(pval)
    e_coli_genus_spec_sign[[genename]] <- pval
    if (pval < 0.05) {
      print(paste0("Genus specific gene: ", genename, "; p-value: ", pval))

    }
  }
  # same for salmonella
  is_genus_specific_salmonella <- sum(cluster_data$genus == "Salmonella" & cluster_data$Essentiality == 1) /
      sum(cluster_data$genus == "Salmonella") - sum(cluster_data$genus != "Salmonella" & cluster_data$Essentiality == 1)/
      sum(cluster_data$genus != "Salmonella")



  if (is_genus_specific_salmonella > 0) {
    print(paste0("Salmonella: ",genename))
    genus_specificity <- c()
    for (i in 1:10000) {
      cluster_data <- cluster_data %>% mutate(genus = sample(genus))
      is_genus_specific_p <- sum(cluster_data$genus == "Salmonella" & cluster_data$Essentiality == 1) /
        sum(cluster_data$genus == "Salmonella") - sum(cluster_data$genus != "Salmonella" & cluster_data$Essentiality == 1)/
        sum(cluster_data$genus != "Salmonella")
      genus_specificity <- c(genus_specificity, is_genus_specific_p)
    }
    # plot the distribution of genus-specificity

    pval <- sum(genus_specificity >= is_genus_specific_salmonella) / length(genus_specificity)
    print(pval)
    salm_genus_spec_sign[[genename]] <- pval
    if (pval < 0.05) {
      print(paste0("Genus specific gene: ", genename, "; p-value: ", pval))
    }
  }

}
# save the p-values of Escherichia coli and Salmonella in one file
write.table(e_coli_genus_spec_sign, file = "./data/e_coli_genus_spec_sign.txt")
write.table(salm_genus_spec_sign, file = "./data/salm_genus_spec_sign.txt")


# import the p-values
e_coli_genus_spec_sign <- read.table("./data/e_coli_genus_spec_sign.txt")
salm_genus_spec_sign <- read.table("./data/salm_genus_spec_sign.txt")

all_pvals_ecoli <- unlist(e_coli_genus_spec_sign)
all_pvals_salmonella <- unlist(salm_genus_spec_sign)

# correct the p-values for multiple testing (p.adjust)
e_coli_genus_spec_sign <- p.adjust(all_pvals_ecoli, method = "fdr")
salm_genus_spec_sign <- p.adjust(all_pvals_salmonella, method = "fdr")













