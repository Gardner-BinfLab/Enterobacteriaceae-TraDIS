from Bio import SeqIO
from re import findall, match
from os import system, path, mkdir
from shutil import rmtree, copyfile


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


# seqdb = '/home/fatemeh/EnTrI/data/fasta-protein/chromosome/seqdb.fasta'
# hieranoid = '/home/fatemeh/EnTrI/results/hieranoid/hieranoid-result.txt'
# outdir = '/home/fatemeh/EnTrI/results/hieranoid/edited-hieranoid-result'
# clusterspath = '/home/fatemeh/EnTrI/results/hieranoid/clusters.txt'

seqdb = '/home/fatemeh/EnTrI/results/deg/fasta-proteins/seqdb.fasta'
hieranoid = '/home/fatemeh/EnTrI/results/deg/hieranoid-result.txt'
outdir = '/home/fatemeh/EnTrI/results/deg/edited-hieranoid-result'
clusterspath = '/home/fatemeh/EnTrI/results/deg/clusters.txt'

makedir(outdir)
sequences = read_fasta_sequences(seqdb)
clusters = list()
with open(hieranoid) as from_file:
    for line in from_file:
        genes = findall('[\(,]([a-zA-Z0-9_]+):', line)
        if genes[0].startswith('DEG'):
            species = [gene[0:7] for gene in genes]
        else:
            species = [s.strip('_') for s in findall('\'([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+\'', str(genes))]
        unique_species = list(set(species))
        while len(unique_species) < len(species):
            uniques_index = []
            newgenes = []
            newspecies = []
            for i in range(len(species)):
                j = 0
                while j in range(len(species)):
                    if i != j and species[i] == species[j]:
                        break
                    else:
                        j += 1
                if j == len(species):
                    uniques_index.append(i)
            if len(uniques_index) == 0:
                uniques_index = [0]
            with open(outdir + '/uniques.fasta', 'w') as uniquesfile:
                with open(outdir + '/others.fasta', 'w') as othersfile:
                    for item in range(len(genes)):
                        if item in uniques_index:
                            uniquesfile.write(sequences[genes[item]].format('fasta'))
                            newgenes.append(genes[item])
                            newspecies.append(species[item])
                        else:
                            othersfile.write(sequences[genes[item]].format('fasta'))
            if len(uniques_index) == 1:
                copyfile(outdir + '/uniques.fasta', outdir + '/mafft.fasta')
            else:
                try:
                    system('mafft --text --quiet {0} > {1}'.format(outdir + '/uniques.fasta',
                                                                   outdir + '/mafft.fasta'))
                except:
                    raise SystemExit
            # try:
            #     system('fasttree {0} > {1}'.format(outdir + '/mafft.fasta',
            #                                        outdir + '/fasttree.tre'))
            # except:
            #     raise SystemExit

            try:
                system('hmmbuild -o /dev/null {0} {1}'.format(outdir + '/hmmbuild.res',
                                                              outdir + '/mafft.fasta'))
            except:
                raise SystemExit
            try:
                system('hmmsearch -o /dev/null --noali --domtblout {0} {1} {2}'.
                       format(outdir + '/hmmsearch.res', outdir + '/hmmbuild.res', outdir + '/others.fasta'))
            except:
                raise SystemExit
            with open(outdir + '/hmmsearch.res') as hmmfile:
                for row in hmmfile:
                    if not row.startswith('#'):
                        cells = row.split()
                        gn = cells[0]
                        sp = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', gn).group(1).strip('_')
                        if sp not in newspecies:
                            newspecies.append(sp)
                            newgenes.append(gn)
            clusters.append(newgenes)
            [genes.remove(g) for g in newgenes]
            species = [s.strip('_') for s in findall('\'([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+\'', str(genes))]
            unique_species = list(set(species))
        clusters.append(genes)
list_of_genes = [item for sublist in clusters for item in sublist]
for item in sequences.keys():
    if item not in list_of_genes:
        clusters.append([item])
with open(clusterspath, 'w') as clustersfile:
    for item in clusters:
        for ii in item:
            clustersfile.write(ii+'\t')
        clustersfile.write('\n')
