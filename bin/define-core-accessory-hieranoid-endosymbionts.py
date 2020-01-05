#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match, findall
from collections import defaultdict


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict


def makedir(dirname):
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

seqdb = '/home/fatemeh/EnTrI/results/endosymbionts/fasta-proteins/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/endosymbionts/clusters.txt'
outdir = '/home/fatemeh/EnTrI/results/endosymbionts/define-core-accessory-hieranoid'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
species_names = ['A359', 'bbp', 'BOBLI757', 'BUMPF009', 'CWO', 'Sant', 'A35E', 'BCc', 'BPEN', 'BUMPG002', 'CWQ', 'SG',
                 'AB162', 'BCHRO640', 'BTURN675', 'BUMPUSDA', 'CWS', 'SOPEG', 'BA000003', 'BCI', 'BUAMB', 'BUMPW106',
                 'CWU', 'WIGMOR', 'BA000021', 'BCTU', 'BUAP5A', 'BUsg', 'IM45', 'BAKON', 'Bfl', 'BUAPTUC7', 'BVAF',
                 'MEPCIT']


esscoredir = outdir + '/core-essential-genomes/'
makedir(esscoredir)
sesscoredir = outdir + '/core-sometimes-essential-genomes/'
makedir(sesscoredir)

# esscoregenes = []
# sesscoregenes = []

num_species = len(species_names)
gene_dict = {species_names[i]: 0 for i in range(num_species)}
with open(clusters) as from_file:
    counter = 1
    for line in from_file:
        for key in gene_dict.keys():
            gene_dict[key] = 0
        list_of_genes = []

        genes = line.split()
        for g in genes:
            s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
            if s in gene_dict.keys():
                list_of_genes.append(g)
                gene_dict[s] = 1
        # if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
        if sum(gene_dict.values()) >= len(gene_dict.keys()) * 80 / 100:
            # esscoregenes += list_of_genes
            with open(esscoredir + 'clust' + str(counter), 'w') as tofile:
                for gene in list_of_genes:
                    tofile.write('>' + gene + '\n')
                    tofile.write(str(sequences[gene].seq) + '\n')
        else:
            with open(sesscoredir + 'clust' + str(counter), 'w') as tofile:
                for gene in list_of_genes:
                    tofile.write('>' + gene + '\n')
                    tofile.write(str(sequences[gene].seq) + '\n')
        counter += 1
            # sesscoregenes += list_of_genes
# esscoregenes = list(set(esscoregenes))
# esscoregenes.sort()
# sesscoregenes = list(set(sesscoregenes))
# sesscoregenes.sort()

# prev_name = ''
# for gene in esscoregenes:
#     name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
#     if name != prev_name:
#         if 'esscoregenesfile' in locals():
#             esscoregenesfile.close()
#         esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
#     esscoregenesfile.write(sequences[gene].format('fasta'))
#     prev_name = name
#
# prev_name = ''
# for gene in sesscoregenes:
#     name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
#     if name != prev_name:
#         if 'sesscoregenesfile' in locals():
#             sesscoregenesfile.close()
#         sesscoregenesfile = open(sesscoredir + '/' + name + '.fasta', 'w')
#     sesscoregenesfile.write(sequences[gene].format('fasta'))
#     prev_name = name