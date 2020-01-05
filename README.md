Here the steps for different analyses in the manuscript are explained.
# Cluster homologous proteins
* cd bin
* homclust.py ../data/fasta-protein/chromosome/Citrobacter_rodentium_ICC168_FN543502_v1.fasta ../data/fasta-protein/chromosome/seqdb.fasta ../results/homclust -ns 14 -s 5
* The output is saved in results/homclust
# Generate phylogenetic tree
* cd results/phylosift
* find ../../data/fasta-genome/chromosome/ -maxdepth 1 -name "*.fa" -exec ~/program-bank/phylosift_v1.0.1/phylosift search --isolate --besthit {} \;
* find ../../data/fasta-genome/chromosome/ -maxdepth 1 -name "*.fa" -exec ~/program-bank/phylosift_v1.0.1/phylosift align --isolate --besthit {} \;
* for f in PS_temp/*/alignDir/concat.updated.1.fasta; do cat "$f" >> protein_alignment.fa; done
* raxmlHPC -s protein_alignment.fa -n phylosift-aa.raxmlbootstrap -m PROTGAMMALG4M -p 1234 -f a -x 1234 -# 100
# Cluster orthologous proteins
* Run hieranoid with BLAST as its similarity search tool on the proteins in data/fasta-protein/chromosome using the tree in results/phylosift/RAxML_bestTree.phylosift-aa.raxmlbootstrap rooted at midpoint. The output of the program is in results/hieranoid
* Convert the output to a file that is easy to read for other programs in this analysis by running edit-hieranoid.py
# Calculate insertion index
* Run calculate-insertion-index.R to calculate insertion index for trimmed genes on both 5' and 3' sides.
* Run calculate-insertion-index_not-trimmed.R to calculate insertion index for untrimmed genes.
* The outputs are saved in results/insertion-indices
# Benchmark the methods of essentiality evaluation
* ecogene-k12.txt in results directory contains list of E.coli K-12 essential genes from ecogenes.
* Run find-k12-counterparts.py
* The results of Monte Carlo method (Turner, K.H., Wessel, A.K., Palmer, G.C., Murray, J.L., and Whiteley, M. (2015). Essential genome of Pseudomonas aeruginosa in cystic fibrosis sputum. PNAS 112, 4110–4115.) are saved in results/Tn-seq
* Run benchmark-essentiality-calls.R. The results are saved in figures/essential-call-comparison-*.pdf
# STUDY BIASES
## Nucleotide composition bias
* Run make-logos-top-100.py to prepare data for making logos from 100 insertion sites with the highest number of insertions.
* Generate logos:
	* results/logos
	* weblogo -F pdf -A dna -f 100logos.txt -o ../../figures/100logo-prob.pdf -s large -U probability
	* weblogo -F pdf -A dna -f 100logos.txt -o ../../figures/100logo-bits.pdf -s large --composition "{'A':23, 'C':27, 'G':27, 'T':23}"
## Position bias within genes
* Run insertion-position-bias.py followed by insertion-position-bias.R. The output figure is saved in figures/insertion-position-bias.pdf
## G-C and distance from origin biases
* Run check-biases.py, the output is saved in results/biases/check-biases
* Run check-biases_not-trimmed.py to check the biases on the genes that are not trimmed on the two ends. The results are saved in results/biases/check-biases_not-trimmed
## Correct biases
* Run check-biases-ii.R to correct the biases. The code saves the corrected insertion indices in results/biases/dbscan and visualises different biases in figures/biases.pdf
* Run compare-bias-predictions.R to compare the impact of bias corrections. The output figure is saved in figures/compare-bias-predictions.pdf
# DBSCAN better clusters essential and non-essential genes compared to gamma fit
* Run compare-dbscan-gamma.R. The resulting figure is saved in figures/gamma-vs-dbscan.pdf
# Find ideal insertion density
* Run density-SL1344.R. The output is saved in false-positive-rate_density.pdf.
# Genus-specific, single-copy, and multi-copy genes
* Run this command to add insertion indices to homologous clusters (the output is saved in results/merge-clust-plot): python merge-clust-plot.py ../results/homclust/EFam-clusters/ ../results/biases/dbscan/ ../results/ecogene-k12.txt ../results/merge-clust-plot
* Run clust2plot.R to plot the genus-specific, single-copy, and multi-copy insertion indices in figures/cluster-essentiality.pdf
# Define core essential and ancestrally essential genes
* To find core essential genes run define-core-accessory-hieranoid-fasta.py. The outputs are saved in results/define-core-accessory-hieranoid-fasta-80
* To find ancestrally essential genes run define-core-accessory-hieranoid-fitch.py. The outputs are saved in results/define-core-accessory-hieranoid-fitch*
* Annotate the phylogenetic tree using results/define-core-accessory-hieranoid-fitch/info.txt
* Run tree-data.R which plots results/fitch.pdf
# Define core essential genes for DEG dataset
* Cluster the genes in DEG dataset using the steps explained in "Generate phylogenetic tree" and "Cluster orthologous proteins". The results are saved in results/deg/clusters.txt
* Run define-core-accessory-hieranoid-deg.py. The output is saved in results/deg/define-core-accessory-hieranoid
# Define core genes for endosymbionts
* Cluster the genes in DEG dataset using the steps explained in "Generate phylogenetic tree" and "Cluster orthologous proteins". The results are saved in results/endosymbionts/clusters.txt
* Run define-core-accessory-hieranoid-endosymbionts.py. The output is saved in results/endosymbionts/define-core-accessory-hieranoid
# Conservation vs. essentiality
* Use the clusters in results/define-core-accessory-hieranoid-fitch-cores as the set of ancestral genes
* Use the essentiality data from results/define-core-accessory-hieranoid-fitch-core-essentials and generate results/walking-hypergeometric-test/essentiality.txt
* Define the conservation level of genes as described in the methods section of the paper
* Save the conservation measure for each gene in results/walking-hypergeometric-test/conservation.out
* Run "Rscript gsea.R -d ../results/walking-hypergeometric-test/essentiality.txt -l ../results/walking-hypergeometric-test/conservation.out -o ../figures". The results are saved in figures/essential folder
