#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, mkdir
from Bio import SeqIO
from shutil import rmtree
from re import match


class Stack:
    def __init__(self):
        self.items = []

    def isEmpty(self):
        return self.items == []

    def push(self, item):
        self.items.append(item)

    def pop(self):
        return self.items.pop()

    def size(self):
        return len(self.items)

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
                if cells[2] == 'essential':
                        iidict[cells[0]] = {1}
                else:
                    iidict[cells[0]] = {0}
    return iidict

def read_k12(inpath, iidict):
    with open(inpath) as from_file:
        for line in from_file:
            cells = line.split()
            iidict[cells[0]] = {1}
    return iidict

seqdb = '/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta'
clusters = '/home/fatemeh/EnTrI/results/hieranoid/clusters.txt'
insertion_indices = '/home/fatemeh/EnTrI/results/biases/dbscan'
k12path = '/home/fatemeh/EnTrI/results/ecogene-k12.txt'
outdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch'
coredir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-cores'
coreessdir = '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch-core-essentials'
makedir(outdir)
makedir(coredir)
makedir(coreessdir)
speciestreedir = '/home/fatemeh/EnTrI/bin/speciestrees-no-k12'
sequences = read_fasta_sequences(seqdb)
gene_essentiality = read_gene_essentiality(insertion_indices)
gene_essentiality = read_k12(k12path, gene_essentiality)
species_names = {"all":["BN373", "ERS227112", "NCTC13441", "ROD", "SEN", "SL1344", "STM",
    "STMMW", "t", "SL3261","BW25113", "EC958", "b"],"typhimurium":["STM", "SL1344", "STMMW", "SL3261"], "salmonella":["SEN", "SL1344", "STM", "STMMW", "t", "SL3261"],
    "ecoli":["NCTC13441","BW25113", "EC958", "b"], "klebsiella":["ERS227112", "BN373"], "citrobacter":["ROD"],
    "salmonellacitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "ROD"],
    "salmonellaecolicitrobacter":["SEN", "SL1344", "SL3261", "STM", "STMMW", "t", "NCTC13441", "ROD","BW25113", "EC958", "b"],
    "salmonellaty2":["t"], "salmonellap125109":["SEN"],
    "salmonellasl1344":["SL1344"], "salmonellasl3261":["SL3261"], "salmonellaa130":["STM"], "salmonellad23580":["STMMW"], "ecolist131":["NCTC13441"],
    "klebsiellarh201207":["ERS227112"], "klebsiellaecl8":["BN373"],"ecoliBW25113":["BW25113"], "ecoliEC958":["EC958"],
    "ecolik12":["b"]}
counter = 0

with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/info.txt', 'w') as infofile:
    infofile.write('speciesname\tcoreessential\tcore\n')

for item in species_names.keys():
    speciesdir = outdir + '/' + item
    makedir(speciesdir)
    esscoredir = speciesdir + '/core-essential-genomes'
    makedir(esscoredir)
    nesscoredir = speciesdir + '/core-never-essential-genomes'
    makedir(nesscoredir)
    aesscoredir = speciesdir + '/core-sometimes-essential-genomes'
    makedir(aesscoredir)
    accessorydir = speciesdir + '/accessory-genomes'
    makedir(accessorydir)

    esscoregenes = []
    nesscoregenes = []
    aesscoregenes = []
    accessorygenes = []

    lenesscoregenes = 0
    lennescoregenes = 0

    num_species = len(species_names[item])
    gene_dict = {species_names[item][i]: {0} for i in range(num_species)}
    essentiality_dict = {species_names[item][i]: {0} for i in range(num_species)}
    ii = []
    presence = []

    with open(clusters) as from_file:
        for line in from_file:
            for key in gene_dict.keys():
                gene_dict[key] = {0}
                essentiality_dict[key] = {0}
            list_of_genes = []
            stack_of_species = Stack()
            stack_of_genes = Stack()

            genes = line.split()
            for g in genes:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    gene_dict[s] = {1}
                    if g in gene_essentiality.keys():
                        essentiality_dict[s] = gene_essentiality[g]
            if num_species > 1:
                with open(speciestreedir + '/' + item + '.tre') as tree_file:
                    treeline = tree_file.readline()
                    sp_name = ''
                    i = 0
                    char = treeline[i]
                    while char != ';':
                        if char == '(':
                            i += 1
                            char = treeline[i]
                        elif char == ',':
                            if sp_name != '':
                                stack_of_species.push(essentiality_dict[sp_name])
                                stack_of_genes.push(gene_dict[sp_name])
                                sp_name = ''
                            i += 1
                            char = treeline[i]
                        elif char == ')':
                            if sp_name != '':
                                stack_of_species.push(essentiality_dict[sp_name])
                                stack_of_genes.push(gene_dict[sp_name])
                                sp_name = ''

                            val1 = stack_of_species.pop()
                            val2 = stack_of_species.pop()
                            ave = val1 & val2
                            if len(ave) == 0:
                                ave = val1 | val2
                            stack_of_species.push(ave)

                            val1 = stack_of_genes.pop()
                            val2 = stack_of_genes.pop()
                            ave = val1 & val2
                            if len(ave) == 0:
                                ave = val1 | val2
                            stack_of_genes.push(ave)

                            i += 1
                            char = treeline[i]
                        else:
                            sp_name += char
                            i += 1
                            char = treeline[i]

                    ii.append(stack_of_species.pop())
                    presence.append(stack_of_genes.pop())
            else:
                ii.append(list(essentiality_dict.values())[0])
                presence.append(list(gene_dict.values())[0])

    gene_dict = {species_names[item][i]: {0} for i in range(num_species)}
    essentiality_dict = {species_names[item][i]: {0} for i in range(num_species)}
    with open(clusters) as from_file:
        index = 0
        for line in from_file:
            gene_name = ''
            for key in gene_dict.keys():
                gene_dict[key] = {0}
                essentiality_dict[key] = {0}
            list_of_genes = []
            genes = line.split()
            for g in genes:
                s = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', g).group(1).strip('_')
                if s in gene_dict.keys():
                    list_of_genes.append(g)
                    gene_dict[s] = {1}
                    if g in gene_essentiality.keys():
                        essentiality_dict[s] = gene_essentiality[g]
                if s == 'b':
                    l = sequences[g].description
                    gene_name = match('(?:[^\[]*\[){4}([^\]]*)', l).group(1)
            union = set()
            for key in essentiality_dict.keys():
                union = union | essentiality_dict[key]
            if num_species > 1 and presence[index] == {1}:
                # if ii[index] == {0}:
                #     nesscoregenes += list_of_genes
                #     lennescoregenes += 1
                if item =='all': # this if statement makes clusters of core genes
                    counter += 1
                    with open (coredir + '/clust'+ str(counter), 'w') as corefile:
                        for clustitem in list_of_genes:
                            corefile.write('>'+clustitem+'\n')
                            corefile.write(str(sequences[clustitem].seq)+'\n')

                if ii[index] == {1}:
                    esscoregenes += list_of_genes
                    lenesscoregenes += 1
                    if item == 'all':
                        with open(
                                '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/always-ess.txt',
                                'a') as efile:
                            efile.write(l + '\n')
                        with open(coreessdir+ '/clust'+ str(counter), 'w') as coreessfile:
                            for clustitem in list_of_genes:
                                coreessfile.write('>'+clustitem+'\n')
                                coreessfile.write(str(sequences[clustitem].seq)+'\n')
                else:
                    lennescoregenes += 1
                    if union != {0}:
                        aesscoregenes += list_of_genes
                        if item == 'all':
                            with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/sometimes-ess.txt',
                                      'a') as sefile:
                                sefile.write(l + '\n')
                    else:
                        nesscoregenes += list_of_genes
                        if item == 'all':
                            with open(
                                    '/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/never-ess.txt',
                                    'a') as nefile:
                                nefile.write(l + '\n')

                index += 1

            elif presence[index] == {1} and num_species == 1:
                if ii[index] == {0}:
                    nesscoregenes += list_of_genes
                    lennescoregenes += 1
                elif ii[index] == {1}:
                    esscoregenes += list_of_genes
                    lenesscoregenes += 1
                # else:
                #     aesscoregenes += list_of_genes
                #     lennescoregenes += 1
                index += 1
            else:
                accessorygenes += list_of_genes
                index += 1

    esscoregenes = list(set(esscoregenes))
    esscoregenes.sort()
    nesscoregenes = list(set(nesscoregenes))
    nesscoregenes.sort()
    accessorygenes = list(set(accessorygenes))
    accessorygenes.sort()
    aesscoregenes = list(set(aesscoregenes))
    aesscoregenes.sort()

    with open('/home/fatemeh/EnTrI/results/define-core-accessory-hieranoid-fitch/info.txt', 'a') as infofile:
        infofile.write(str(item) + '\t' + str(lenesscoregenes) + '\t' +
                       str(lennescoregenes+lenesscoregenes) + '\n')

    prev_name = ''
    for gene in esscoregenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'esscoregenesfile' in locals():
                esscoregenesfile.close()
            esscoregenesfile = open(esscoredir + '/' + name + '.fasta', 'w')
        esscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in nesscoregenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'nesscoregenesfile' in locals():
                nesscoregenesfile.close()
            nesscoregenesfile = open(nesscoredir + '/' + name + '.fasta', 'w')
        nesscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in aesscoregenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'aesscoregenesfile' in locals():
                aesscoregenesfile.close()
            aesscoregenesfile = open(aesscoredir + '/' + name + '.fasta', 'w')
        aesscoregenesfile.write(sequences[gene].format('fasta'))
        prev_name = name

    prev_name = ''
    for gene in accessorygenes:
        match_result = match('([a-zA-Z0-9]+?)_\S+', gene)
        if not match_result:
            match_result = match('([a-zA-Z]+?)\d+\S*', gene)
        name = match_result.group(1)
        if name != prev_name:
            if 'accessorygenesfile' in locals():
                accessorygenesfile.close()
            accessorygenesfile = open(accessorydir + '/' + name + '.fasta', 'w')
        accessorygenesfile.write(sequences[gene].format('fasta'))
        prev_name = name