from collections import defaultdict
from re import match, findall
from os import listdir, path, mkdir
from shutil import rmtree

def makedir(dirname):
    # print ('Making directory \'{0}\'..'.format(dirname))
    if path.exists(dirname):
        input_var = 'i'
        while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
            input_var = input('Directory {0} already exists. Replace? [Y/n] '.format(dirname))
        if input_var == 'Y' or input_var == 'y':
            rmtree(dirname)
        else:
            raise SystemExit
    mkdir(dirname)

seqdb = '../data/fasta-dna/chromosome/seqdb.fasta'
# iis = '../results/maximise_MCC/pca'
iis = '../results/insertion-indices/gamma'
results = '../results/biases/check-biases'
makedir(results)
genome_length = {"SL1344": 4878012, "STMMW": 4879400, "SEN": 4685848, "t": 4791961, "STM": 4895639, "ETEC": 5153435,
                 "b": 4641652, "CS17": 4994793, "NCTC13441": 5174631, "ROD": 5346659, "BN373": 5324709,
                 "ERS227112": 5869288, "ENC": 4908759, "SL3261": 4878012, "EC958": 5109767,
                 "BW25113": 4631469}
dnaa ={"ROD":4262871, "CS17":4234263, "ENC":414484, "ETEC":4305897, "NCTC13441":4952702, "ERS227112":453004,
       "BN373": 5024509, "SEN":3919680, "STM":4019091, "SL1344":4066338, "STMMW":4067900, "t":3790618,
       "SL3261":4066338, "b":3883729, "BW25113":3875686, "EC958":4240812}

list_of_files = listdir(iis)
for filename in list_of_files:
    iis_dict = defaultdict(tuple)
    with open('{0}/{1}'.format(iis, filename), 'r') as iifile:
        for line in iifile:
            cells = line.split()
            iis_dict[cells[0]] = (float(cells[1]), cells[2])

    gene_name = ''
    with open(results + '/' + filename, 'w') as tofile:
        with open(seqdb, 'r') as sequencefile:
            for line in sequencefile:
                if line.startswith('>'):
                    if gene_name != '':
                        dist = min(abs(start[0]-dnaa[strain_name]), genome_length[strain_name]-abs(start[0]-dnaa[strain_name])) + 1
                        tofile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n'.format(gene_name, iis_dict[gene_name][0],
                                                                             iis_dict[gene_name][1],
                                                                             float(dist) / genome_length[strain_name],
                                                                             float(gc) / gene_length,
                                                                             float(gene_length), float(start[0])))
                        # tofile.write('{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n'.format(gene_name, iis_dict[gene_name][0], iis_dict[gene_name][1], float(dist)/genome_length[strain_name], float(gc)/gene_length, float(gene_length)))
                    gc = 0
                    match_result = match('>\s*((\S+?)_+\S+)\s+\[\S+/((\d+\-\d+\s)+)\(', line)
                    if match_result is None:
                        match_result = match('>\s*(([a-zA-Z]+)\d+[a-zA-Z]?)\s+\[\S+/((\d+\-\d+\s)+)\(', line)
                    if match_result is not None:
                        gene_name = match_result.group(1)
                        strain_name = match_result.group(2)
                        starts_ends = findall('\d+\-\d+\s', match_result.group(3))
                        start = []
                        end = []
                        for item in starts_ends:
                            match_result = match('(\d+)-(\d+)', item)
                            start.append(int(match_result.group(1)))
                            end.append(int(match_result.group(2)))
                    else:
                        gene_name = ''
                    if gene_name not in iis_dict.keys():
                        gene_name = ''
                    else:
                        # gene_insertions = 0
                        gene_length = 0
                        for i in range(0,len(start)):
                            # gene_insertions += sum([1 for x in iis_dict[strain_name][start[i]-1:end[i]] if x > 0])
                            gene_length += end[i] - start[i] + 1
                        # ii = (float(gene_insertions)/gene_length)/(float(genome_insertions[strain_name])/genome_length[strain_name])
                else:
                    gc += line.count('g') + line.count('c')