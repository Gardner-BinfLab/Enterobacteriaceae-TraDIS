from Bio import SeqIO
from re import findall, match
from os import system, path, mkdir
from shutil import rmtree, copyfile


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict

clustpath = '../results/all-hieranoid/clusters.txt'
annotpath = '../results/all-hieranoid/cluster-representatives.txt'
seqpath = '../results/all-hieranoid/seqdb.fasta'
# large_clustpath = '../results/all-hieranoid/large-clusters.txt'

sequences = read_fasta_sequences(seqpath)

with open(clustpath, 'r') as fromfile:
    with open(annotpath, 'w') as annotfile:
        # with open(large_clustpath, 'w') as clustfile:
        for line in fromfile:
            cells = line.split()
            cells = [item for item in cells if not item.startswith('exDEG')]
            if len(cells) > 0:
                # for item in cells[:-1]:
                #     clustfile.write(item + '\t')
                # clustfile.write(cells[-1] + '\n')
                rep = findall('\'(b[0-9]+)\'', str(cells))
                if not rep:
                    rep = findall('\'(DEG[0-9]+)\'', str(cells))
                if not rep:
                    rep = [cells[0]]
                rep = rep[0]
                annotfile.write(sequences[rep].format('fasta'))
