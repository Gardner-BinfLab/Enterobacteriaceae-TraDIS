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


def read_gene_essentiality(indir):
    list_of_files = listdir(indir)
    iidict = {}
    for filename in list_of_files:
        with open(indir+'/'+filename) as from_file:
            for line in from_file:
                cells = line.split()
                iidict[cells[0]] = cells[2]
    return iidict


def read_k12(inpath, iidict):
    with open(inpath) as from_file:
        for line in from_file:
            cells = line.split()
            iidict[cells[0]] = 'essential'
    return iidict

seqdb = '/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/hieranoid/clusters.txt'
insertion_indices = '/home/fatemeh/EnTrI/results/biases/dbscan'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fasta-80'
makedir(outdir)
sequences = read_fasta_sequences(seqdb)
gene_essentiality = read_gene_essentiality(insertion_indices)
gene_essentiality = read_k12(k12path, gene_essentiality)
species_names = defaultdict()
species_names = {"all":["BN373", "ERS227112", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "b", "SL3261","BW25113", "EC958"],"typhimurium":["STM", "SL1344", "STMMW", "SL3261"],
    "salmonella":["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"], "ecoli":["NCTC13441", "b","BW25113", "EC958"],
    "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"], "salmonellacitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "ROD"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "NCTC13441", "b", "ROD","BW25113", "EC958"],
    "salmonellaty2":["t"], "salmonellap125109":["SEN"], "salmonellasl1344":["SL1344"], "salmonellasl3261":["SL3261"],
    "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"], "ecolik12":["b"],
    "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"],"ecoliBW25113":["BW25113"], "ecoliEC958":["EC958"]}

with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fasta-80/info.txt', 'w') as infofile:
    infofile.write('speciesname\tcoreessential\tcore\n')

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes'
    makedir(esscoredir)
    sesscoredir = speciesdir + '/core-sometimes-essential-genomes'
    makedir(sesscoredir)
    nesscoredir = speciesdir + '/core-never-essential-genomes'
    makedir(nesscoredir)
    essaccessorydir = speciesdir + '/accessory-essential-genomes'
    makedir(essaccessorydir)
    sessaccessorydir = speciesdir + '/accessory-sometimes-essential-genomes'
    makedir(sessaccessorydir)
    nessaccessorydir = speciesdir + '/accessory-never-essential-genomes'
    makedir(nessaccessorydir)

    esscoregenes = []
    sesscoregenes = []
    nesscoregenes = []
    essaccessorygenes = []
    sessaccessorygenes = []
    nessaccessorygenes = []

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: 0 for i in range(num_species)}
    essentiality_dict = {species_names[item][i]: 0 for i in range(num_species)}
    # for filename in list_of_files:
    with open (clusters) as from_file:
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = 0
                essentiality_dict[key] = 0
            list_of_genes = []
            list_of_essential_genes = []

            genes = line.split()
            # findall_result = findall('(([a-zA-Z0-9]+?)_[a-zA-Z0-9]+):', line)
            # temp_result = findall('(([a-zA-Z]+?)\d+):', line)
            # for matches in temp_result:
            #     if matches[1] != 'n':
            #         findall_result.append(matches)
            for g in genes:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    gene_dict[s] = 1
                    if g in gene_essentiality and gene_essentiality[g] == 'essential':
                        list_of_essential_genes.append(g)
                        essentiality_dict[s] = 1
            #if gene_dict[min(gene_dict, key=gene_dict.get)] == 1:
            if(sum(gene_dict.values()) >= 80/100 * len(gene_dict.keys())):
                #if essentiality_dict[min(essentiality_dict, key=essentiality_dict.get)] == 1:
                #if len(list_of_essential_genes) == len(list_of_genes):
                if len(list_of_essential_genes) >= len(list_of_genes) * 80 / 100:
                    esscoregenes += list_of_genes
                elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
                    nesscoregenes += list_of_genes
                else:
                    sesscoregenes += list_of_genes
            else:
                #if not sum([gene_dict[item] - essentiality_dict[item] for item in essentiality_dict.keys()]):
                if len(list_of_essential_genes) >= len(list_of_genes) * 80 / 100:
                    essaccessorygenes += list_of_genes
                elif essentiality_dict[max(essentiality_dict, key=essentiality_dict.get)] == 0:
                    nessaccessorygenes += list_of_genes
                else:
                    sessaccessorygenes += list_of_genes
    esscoregenes = list(set(esscoregenes))
    esscoregenes.sort()
    sesscoregenes = list(set(sesscoregenes))
    sesscoregenes.sort()
    nesscoregenes = list(set(nesscoregenes))
    nesscoregenes.sort()
    essaccessorygenes = list(set(essaccessorygenes))
    essaccessorygenes.sort()
    sessaccessorygenes = list(set(sessaccessorygenes))
    sessaccessorygenes.sort()
    nessaccessorygenes = list(set(nessaccessorygenes))
    nessaccessorygenes.sort()

    with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fasta-80/info.txt', 'a') as infofile:
        infofile.write(str(item) + '\t' + str(len(esscoregenes)/len(species_names[item])) + '\t' +
                       str((len(esscoregenes)+len(sesscoregenes)+len(nesscoregenes))/len(species_names[item])) + '\n')

    prev_name = ''
    for gene in esscoregenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'esscoregenesfile' in locals():
                esscoregenesfile.close()
            esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
        esscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in sesscoregenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'sesscoregenesfile' in locals():
                sesscoregenesfile.close()
            sesscoregenesfile = open(sesscoredir + '/' + name + '.fasta', 'w')
        sesscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in nesscoregenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'nesscoregenesfile' in locals():
                nesscoregenesfile.close()
            nesscoregenesfile = open(nesscoredir + '/' + name + '.fasta', 'w')
        nesscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in essaccessorygenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'essaccessorygenesfile' in locals():
                essaccessorygenesfile.close()
            essaccessorygenesfile = open(essaccessorydir + '/' + name + '.fasta', 'w')
        essaccessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in sessaccessorygenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'sessaccessorygenesfile' in locals():
                sessaccessorygenesfile.close()
            sessaccessorygenesfile = open(sessaccessorydir + '/' + name + '.fasta', 'w')
        sessaccessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in nessaccessorygenes:
        name = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gene).group(1).strip('_')
        if name != prev_name:
            if 'nessaccessorygenesfile' in locals():
                nessaccessorygenesfile.close()
            nessaccessorygenesfile = open(nessaccessorydir + '/' + name + '.fasta', 'w')
        nessaccessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name
