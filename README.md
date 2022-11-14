Here the steps for different analyses in the manuscript are explained.
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
# Ancestrally essential genes
* To find ancestrally essential genes run define-core-accessory-hieranoid-fitch.py. The outputs are saved in results/define-core-accessory-hieranoid-fitch*
# Core essential genes for TraDIS, DEG and endosymbionts
## Generate phylogenetic tree
* cd EnTrI/results/all-hieranoid/
* find fasta-genome/ -maxdepth 1 -name "*.fa*" -exec ~/program-bank/phylosift_v1.0.1/phylosift search --isolate --besthit {} \;
* find fasta-genome/ -maxdepth 1 -name "*.fa*" -exec ~/program-bank/phylosift_v1.0.1/phylosift align --isolate --besthit {} \;
* for f in PS_temp/*/alignDir/concat.updated.1.fasta; do cat "$f" >> protein_alignment.fa; done
* raxmlHPC -s protein_alignment.fa -n phylosift-aa.raxmlbootstrap -m PROTGAMMALG4M -p 1234 -f a -x 1234 -# 100
## Cluster orthologous proteins
* Run hieranoid with BLAST as its similarity search tool on the proteins in data/fasta-protein/chromosome using the tree in results/phylosift/RAxML_bestTree.phylosift-aa.raxmlbootstrap rooted at midpoint. The output of the program is in results/all-hieranoid/hieranoid-result.txt
* Convert the output to a file that is easy to read for other programs in this analysis by running edit-hieranoid.py
## Generate tables and figures summarising essentiality
* Run annotate-clust.py. The output is saved in results/all-hieranoid/cluster-representatives.txt
* Run eggnog-mapper on cluster-representatives.txt to find gene names for clusters.
* Download KEGG annotations and save them in results/KEGG
* Run giant-tab-1.py
* Run giant-tab-2.R
* Run venn-ancestral-core-endosymbiont.R


