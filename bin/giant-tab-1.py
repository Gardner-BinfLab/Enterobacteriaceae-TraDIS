from os import listdir, path, mkdir
from shutil import rmtree
from re import match
import pandas as pd

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

def read_k12(inpath):
    all = []
    with open(inpath) as from_file:
        for line in from_file:
            cell = line.split()
            all.append(cell[0])
    return all

def read_ancestral(inpath):
    genes = []
    with open(inpath, 'r') as fromfile:
        for line in fromfile:
            if line.startswith('>'):
                genes.append(line.split(' ')[0][1:])
    return(genes)

outdir = '../results/giant-tab/'
makedir(outdir)
outpath = outdir + 'giant-tab.tsv'

colnames = ['DEG: Synechococcus elongatus PCC 7942', 'DEG: Mycobacterium tuberculosis H37Rv',
            'DEG: Porphyromonas gingivalis ATCC 33277', 'DEG: Bacteroides thetaiotaomicron VPI-5482',
            'DEG: Bacteroides fragilis 638R', 'DEG: Francisella novicida U112', 'DEG: Acinetobacter baumannii ATCC 17978',
            'DEG: Pseudomonas aeruginosa PAO1', 'DEG: Shewanella oneidensis MR-1', 'TraDIS Essentiality: Klebsiella pneumoniae Ecl8',
            'TraDIS Presence: Klebsiella pneumoniae Ecl8', 'TraDIS Essentiality: Klebsiella pneumoniae RH201207',
            'TraDIS Presence: Klebsiella pneumoniae RH201207', 'TraDIS Essentiality: Escherichia coli ST131 EC958',
            'TraDIS Presence: Escherichia coli ST131 EC958', 'TraDIS Essentiality: Escherichia coli UPEC ST131 NCTC13441',
            'TraDIS Presence: Escherichia coli UPEC ST131 NCTC13441', 'EcoGene Essentiality: Escherichia coli BW25113',
            'EcoGene Presence: Escherichia coli BW25113',
            'TraDIS Essentiality: Escherichia coli BW25113', 'TraDIS Presence: Escherichia coli BW25113',
            'TraDIS Essentiality: Citrobacter rodentium ICC168', 'TraDIS Presence: Citrobacter rodentium ICC168',
            'TraDIS Essentiality: Salmonella Typhi Ty2', 'TraDIS Presence: Salmonella Typhi Ty2',
            'DEG: Salmonella enterica serovar Typhi', 'TraDIS Essentiality: Salmonella Typhimurium A130',
            'TraDIS Presence: Salmonella Typhimurium A130', 'TraDIS Essentiality: Salmonella Typhimurium D23580',
            'TraDIS Presence: Salmonella Typhimurium D23580', 'TraDIS Essentiality: Salmonella Typhimurium SL3261',
            'TraDIS Presence: Salmonella Typhimurium SL3261', 'TraDIS Essentiality: Salmonella Typhimurium SL1344',
            'TraDIS Presence: Salmonella Typhimurium SL1344', 'TraDIS Essentiality: Salmonella Enteritidis P125109',
            'TraDIS Presence: Salmonella Enteritidis P125109', 'Symbiont: Secondary endosymbiont of Ctenarytaina eucalypti',
            'Symbiont: Candidatus Moranella endobia PCIT', 'Symbiont: Secondary endosymbiont of Heteropsylla cubana',
            'Symbiont: Candidatus Baumannia cicadellinicola strain BGSS', 'Symbiont: Candidatus Baumannia cicadellinicola strain B-GSS',
            'Symbiont: Baumannia cicadellinicola str. Hc (Homalodisca coagulata)',
            'Symbiont: Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757',
            'Symbiont: Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675',
            'Symbiont: Candidatus Blochmannia vafer str. BVAF', 'Symbiont: Candidatus Blochmannia floridanus',
            'Symbiont: Candidatus Blochmannia pennsylvanicus str. BPEN', 'Symbiont: Candidatus Blochmannia chromaiodes str. 640',
            'Symbiont: Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
            'Symbiont: Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans',
            'Symbiont: Buchnera aphidicola str. Sg (Schizaphis graminum)', 'Symbiont: Buchnera aphidicola str. G002 (Myzus persicae)',
            'Symbiont: Buchnera aphidicola str. USDA (Myzus persicae)', 'Symbiont: Buchnera aphidicola str. F009 (Myzus persicae)',
            'Symbiont: Buchnera aphidicola str. W106 (Myzus persicae)', 'Symbiont: Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)',
            'Symbiont: Buchnera aphidicola str. APS (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)', 'Symbiont: Buchnera aphidicola str. 5A (Acyrthosiphon pisum)',
            'Symbiont: Buchnera aphidicola str. Ua (Uroleucon ambrosiae)', 'Symbiont: Buchnera aphidicola str. Bp (Baizongia pistaciae)',
            'Symbiont: Buchnera aphidicola BCc', 'Symbiont: Buchnera aphidicola (Cinara tujafilina)',
            'Symbiont: Sodalis glossinidius str. morsitans', 'Symbiont: Sodalis praecaptivus strain HS1',
            'Symbiont: Candidatus Sodalis pierantonius str. SOPE', 'DEG: Haemophilus influenzae Rd KW20', 'DEG: Vibrio cholerae N16961',
            'DEG: Burkholderia pseudomallei K96243', 'DEG: Burkholderia thailandensis E264', 'DEG: Rhodopseudomonas palustris CGA009',
            'DEG: Agrobacterium fabrum str. C58', 'DEG: Caulobacter crescentus', 'DEG: Brevundimonas subvibrioides ATCC 15264',
            'DEG: Helicobacter pylori 26695', 'DEG: Mycoplasma genitalium G37', 'DEG: Mycoplasma pulmonis UAB CTIP',
            'DEG: Streptococcus agalactiae A909', 'DEG: Streptococcus pyogenes NZ131', 'DEG: Streptococcus pyogenes MGAS5448',
            'DEG: Staphylococcus aureus NCTC 8325', 'DEG: Staphylococcus aureus N315', 'DEG: Bacillus subtilis 168',
            'Locus: Klebsiella pneumoniae Ecl8',
            'Locus: Klebsiella pneumoniae RH201207', 'Locus: Escherichia coli ST131 EC958',
            'Locus: Escherichia coli UPEC ST131 NCTC13441',
            'Locus: Escherichia coli BW25113 (Keio)', 'Locus: Escherichia coli BW25113',
            'Locus: Citrobacter rodentium ICC168',
            'Locus: Salmonella Typhi Ty2',
            'Locus: Salmonella Typhimurium A130',
            'Locus: Salmonella Typhimurium D23580', 'Locus: Salmonella Typhimurium SL3261',
            'Locus: Salmonella Typhimurium SL1344',
            'Locus: Salmonella Enteritidis P125109', 'Locus: Secondary endosymbiont of Ctenarytaina eucalypti',
            'Locus: Candidatus Moranella endobia PCIT', 'Locus: Secondary endosymbiont of Heteropsylla cubana',
            'Locus: Candidatus Baumannia cicadellinicola strain BGSS', 'Locus: Candidatus Baumannia cicadellinicola strain B-GSS',
            'Locus: Baumannia cicadellinicola str. Hc (Homalodisca coagulata)',
            'Locus: Blochmannia endosymbiont of Camponotus (Colobopsis) obliquus strain 757',
            'Locus: Blochmannia endosymbiont of Polyrhachis (Hedomyrma) turneri strain 675',
            'Locus: Candidatus Blochmannia vafer str. BVAF', 'Locus: Candidatus Blochmannia floridanus',
            'Locus: Candidatus Blochmannia pennsylvanicus str. BPEN', 'Locus: Candidatus Blochmannia chromaiodes str. 640',
            'Locus: Wigglesworthia glossinidia endosymbiont of Glossina brevipalpis',
            'Locus: Wigglesworthia glossinidia endosymbiont of Glossina morsitans morsitans',
            'Locus: Buchnera aphidicola str. Sg (Schizaphis graminum)', 'Locus: Buchnera aphidicola str. G002 (Myzus persicae)',
            'Locus: Buchnera aphidicola str. USDA (Myzus persicae)', 'Locus: Buchnera aphidicola str. F009 (Myzus persicae)',
            'Locus: Buchnera aphidicola str. W106 (Myzus persicae)', 'Locus: Buchnera aphidicola str. Ak (Acyrthosiphon kondoi)',
            'Locus: Buchnera aphidicola str. APS (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. LL01 (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. TLW03 (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. JF98 (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. JF99 (Acyrthosiphon pisum)', 'Locus: Buchnera aphidicola str. 5A (Acyrthosiphon pisum)',
            'Locus: Buchnera aphidicola str. Ua (Uroleucon ambrosiae)', 'Locus: Buchnera aphidicola str. Bp (Baizongia pistaciae)',
            'Locus: Buchnera aphidicola BCc', 'Locus: Buchnera aphidicola (Cinara tujafilina)',
            'Locus: Sodalis glossinidius str. morsitans', 'Locus: Sodalis praecaptivus strain HS1',
            'Locus: Candidatus Sodalis pierantonius str. SOPE'
            ]

tags = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008', 'DEG1006',
        'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002', 'DEG1001']
# genome_dict = dict(zip(tags, colnames))
tradis = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'STM', 'STMMW', 'SL3261',
                'SL1344', 'SEN']
enterobacter = ['BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
                'SL3261', 'SL1344', 'SEN']
endosymbiont = ['A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
                'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON',
                'BA000003', 'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG',
                'Sant', 'SOPEG']
gammaprot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003']
prot = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373',
        'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344',
        'SEN', 'A359', 'MEPCIT', 'A35E', 'IM45', 'AB162', 'BCI', 'BOBLI757', 'BTURN675', 'BVAF', 'Bfl', 'BPEN',
        'BCHRO640', 'BA000021', 'WIGMOR', 'BUsg', 'BUMPG002', 'BUMPUSDA', 'BUMPF009', 'BUMPW106', 'BAKON', 'BA000003',
        'BUAPTUC7', 'CWO', 'CWQ', 'CWU', 'CWS', 'BUAP5A', 'BUAMB', 'bbp', 'BCc', 'BCTU', 'SG', 'Sant', 'SOPEG',
        'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
gammaprotexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b',
                   'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003']
protexsymb = ['DEG1012', 'DEG1043', 'DEG1036', 'DEG1029', 'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113',
              'ROD', 't', 'DEG1016', 'STM', 'STMMW', 'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035',
              'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020', 'DEG1046', 'DEG1008']
bactexsymb = ['DEG1040', 'DEG1027', 'DEG1022', 'DEG1023', 'DEG1034', 'DEG1012', 'DEG1043', 'DEG1036', 'DEG1029',
              'BN373', 'ERS227112', 'EC958', 'NCTC13441', 'b', 'BW25113', 'ROD', 't', 'DEG1016', 'STM', 'STMMW',
              'SL3261', 'SL1344', 'SEN', 'DEG1005', 'DEG1003', 'DEG1035', 'DEG1024', 'DEG1041', 'DEG1045', 'DEG1020',
              'DEG1046', 'DEG1008', 'DEG1006', 'DEG1014', 'DEG1042', 'DEG1038', 'DEG1037', 'DEG1017', 'DEG1002',
              'DEG1001']
lentradis = len(tradis)
lenenterobacter = len(enterobacter)
lenendosymbiont = len(endosymbiont)
lengammaprot = len(gammaprot)
lenprot = len(prot)
lengammaprotexsymb = len(gammaprotexsymb)
lenprotexsymb = len(protexsymb)
lenbact = len(tags)
lenbactexsymb = len(bactexsymb)

essentiality = '../results/biases/dbscan/'
ii_dict = dict()
list_of_files = listdir(essentiality)
for filename in list_of_files:
    with open(essentiality + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            ii_dict[cells[0]] = cells[3]

annotations = '../results/all-hieranoid/cluster-representatives.emapper.annotations'
annot_dict = {}
with open(annotations, 'r') as fromfile:
    for line in fromfile:
        if not line.startswith('#'):
            cells = line.split('\t')
            annot_dict[cells[0]] = cells[5]
list_annot_keys = list(annot_dict.keys())

keioannot = '../data/fasta-protein/chromosome/U00096.fasta'
keio_dict = {}
with open(keioannot, 'r') as fromfile:
    for line in fromfile:
        if line.startswith('>'):
            loc = line.split(' ')[0][1:]
            name = line.split('] [')[3]
            keio_dict[loc] = name


k12path = '../results/ecogene-k12.txt'
k12genes = read_k12(k12path)

kegg = '../results/KEGG/escherichia_coli_K-12_MG1655.dat'
kegg_dict = {}
acc_dict = {}
with open(kegg, 'r') as fromfile:
    for line in fromfile:
        cells = line.split('\t')
        cells[0] = cells[0][1:-1]
        cells[1] = cells[1][1:-1]
        cells[2] = cells[2][1:-2]
        if cells[0] != 'gene_id' and cells[0] not in kegg_dict.keys():
            cells[2] = match('(.*) - Escherichia coli K-12 MG1655', cells[2]).group(1)
            kegg_dict[cells[0]] = cells[2]
            acc_dict[cells[0]] = cells[1]
        elif cells[0] != 'gene_id':
            cells[2] = match('(.*) - Escherichia coli K-12 MG1655', cells[2]).group(1)
            kegg_dict[cells[0]] += ' / ' + cells[2]
            acc_dict[cells[0]] += ' / ' + cells[1]

anc_all = read_ancestral('../results/define-core-accessory-hieranoid-fitch/all/core-essential-genomes/BN373.fasta')
anc_kle = read_ancestral('../results/define-core-accessory-hieranoid-fitch/klebsiella/core-essential-genomes/BN373.fasta')
anc_eco = read_ancestral('../results/define-core-accessory-hieranoid-fitch/ecoli/core-essential-genomes/b.fasta')
anc_sal = read_ancestral('../results/define-core-accessory-hieranoid-fitch/salmonella/core-essential-genomes/SL1344.fasta')
anc_salcit = read_ancestral('../results/define-core-accessory-hieranoid-fitch/salmonellacitrobacter/core-essential-genomes/SL1344.fasta')
anc_salecocit = read_ancestral('../results/define-core-accessory-hieranoid-fitch/salmonellaecolicitrobacter/core-essential-genomes/b.fasta')
# This is used to add ancestrally essential genes that are not annotated or are present in less than three genomes:
anc = set(anc_all) | set(anc_kle) | set(anc_eco) | set(anc_sal) | set(anc_salcit) | set(anc_salecocit)

clusters = '../results/all-hieranoid/clusters.txt'
#prev_name = ''
with open(outpath, 'w') as tofile:
    tofile.write('Cluster\tGene_EGGNOG\tGene_Keio\tKEGG_acc\tKEGG_pathway\t')
    for item in colnames:
        tofile.write(item + '\t')
    tofile.write('Enterobacteriaceae (excluding Salmonella enterica serovar Typhi) %essential\t' +
                 'Enterobacteriaceae %essential\tEndosymbiont %conserved\tGammaproteobacteria %essential\t' +
                 'Gammaproteobacteria (excluding symbionts) %essential\tProteobacteria %essential\t' +
                 'Proteobacteria (excluding symbionts) %essential\tBacteria %essential\t' +
                 'Bacteria (excluding symbionts) %essential\tEnterobacteriaceae: Ancestral essentiality\t' +
                 'E. coli + Salmonella + Citrobacter: Ancestral essentiality\t' +
                 'Salmonella + Citrobacter: Ancestral essentiality\t' +
                 'Klebsiella: Ancestral essentiality\tE. coli: Ancestral essentiality\t' +
                 'Salmonella: Ancestral essentiality\n'
                 )
    with open(clusters, 'r') as fromfile:
        for line in fromfile:
            pathway = '-'
            acc = '-'
            genes = line.split()
            genes = [item for item in genes if not item.startswith('exDEG')]
            genes = [item for item in genes if not item.startswith('ENC_')]
            if len(genes) > 0:
                tofile.write(','.join(genes) + '\t')
                intersect = list(set(genes) & set(list_annot_keys))
                if len(intersect) != 1:
                    intersect = '-'
                else:
                    intersect = intersect[0]
                clustdict = {}
                for item in tags:
                    if item in tradis and item != 'b':
                        clustdict[item] = '-'
                    else:
                        clustdict[item] = '0'
                locusdict = {el: '-' for el in tradis + endosymbiont}

                numtradis = 0
                numenterobacter = 0
                numendosymbiont = 0
                numgammaprot = 0
                numprot = 0
                numgammaprotexsymb = 0
                numprotexsymb = 0
                numbact = 0
                numbactexsymb = 0
                if intersect != '-':
                    tofile.write(annot_dict[intersect] + '\t')
                    if intersect in keio_dict.keys():
                        tofile.write(keio_dict[intersect])
                    else:
                        tofile.write('-')
                else:
                    tofile.write('-\t-')
                tofile.write('\t')
                keio_pres = '0'
                for gene in genes:
                    if gene.startswith('DEG'):
                        species = gene[0:7]
                    else:
                        species = match('([a-zA-Z0-9]+_|[a-zA-Z]+)[a-zA-Z0-9]+', str(gene)).group(1).strip('_')
                    if species in tradis:
                        locusdict[species] = gene
                        if species != 'b':
                            if gene in ii_dict.keys():
                                clustdict[species] = str(ii_dict[gene])
                                if float(ii_dict[gene]) <= 0:
                                    numtradis += 1
                                    if species in enterobacter:
                                        numenterobacter += 1
                                    if species in endosymbiont:
                                        numendosymbiont += 1
                                    if species in gammaprot:
                                        numgammaprot += 1
                                    if species in gammaprotexsymb:
                                        numgammaprotexsymb += 1
                                    if species in prot:
                                        numprot += 1
                                    if species in protexsymb:
                                        numprotexsymb += 1
                                    if species in tags:
                                        numbact += 1
                                    if species in bactexsymb:
                                        numbactexsymb += 1
                        else:
                            keio_pres = '1'
                            if gene in k12genes:
                                clustdict[species] = '1'
                                numtradis += 1
                                if species in enterobacter:
                                    numenterobacter += 1
                                if species in endosymbiont:
                                    numendosymbiont += 1
                                if species in gammaprot:
                                    numgammaprot += 1
                                if species in gammaprotexsymb:
                                    numgammaprotexsymb += 1
                                if species in prot:
                                    numprot += 1
                                if species in protexsymb:
                                    numprotexsymb += 1
                                if species in tags:
                                    numbact += 1
                                if species in bactexsymb:
                                    numbactexsymb += 1
                            # else:
                            #     clustdict[species] = 'N'
                            if gene in kegg_dict.keys():
                                pathway = kegg_dict[gene]
                                acc = acc_dict[gene]
                    else:
                        if species in endosymbiont:
                            locusdict[species] = gene
                        clustdict[species] = '1'

                        if species in enterobacter:
                            numenterobacter += 1
                        if species in endosymbiont:
                            numendosymbiont += 1
                        if species in gammaprot:
                            numgammaprot += 1
                        if species in gammaprotexsymb:
                            numgammaprotexsymb += 1
                        if species in prot:
                            numprot += 1
                        if species in protexsymb:
                            numprotexsymb += 1
                        if species in tags:
                            numbact += 1
                        if species in bactexsymb:
                            numbactexsymb += 1

                tofile.write(acc + '\t' + pathway + '\t')
                for item in tags:
                    tofile.write(clustdict[item] + '\t')
                    if item in tradis and item != 'b':
                        presence = '1'
                        if clustdict[item] == '-':
                            presence = '0'
                        tofile.write(presence + '\t')
                    elif item == 'b':
                        tofile.write(keio_pres + '\t')
                for item in tradis:
                    tofile.write(locusdict[item] + '\t')
                for item in endosymbiont:
                    tofile.write(locusdict[item] + '\t')
                percenttradis = round(numtradis * 100 / lentradis, 1)
                percenterobacter = round(numenterobacter * 100 / lenenterobacter, 1)
                percendosymbiont = round(numendosymbiont * 100 / lenendosymbiont, 1)
                percgammaprot = round(numgammaprot * 100 / lengammaprot, 1)
                percgammaprotexsymb = round(numgammaprotexsymb * 100 / lengammaprotexsymb, 1)
                percprot = round(numprot * 100 / lenprot, 1)
                percprotexsymb = round(numprotexsymb * 100 / lenprotexsymb, 1)
                percbact = round(numbact * 100 / lenbact, 1)
                percbactexsymb = round(numbactexsymb * 100 / lenbactexsymb, 1)
                tofile.write(str(percenttradis) + '\t' + str(percenterobacter) + '\t' + str(percendosymbiont) + '\t'
                             + str(percgammaprot) + '\t' + str(percgammaprotexsymb) + '\t' + str(percprot) + '\t' +
                             str(percprotexsymb) + '\t' + str(percbact) + '\t' + str(percbactexsymb) + '\t')
                anc_all_yn = int(len(set(genes) & set(anc_all)) > 0)
                anc_salecocit_yn = int(len(set(genes) & set(anc_salecocit)) > 0)
                anc_salcit_yn = int(len(set(genes) & set(anc_salcit)) > 0)
                anc_kle_yn = int(len(set(genes) & set(anc_kle)) > 0)
                anc_eco_yn = int(len(set(genes) & set(anc_eco)) > 0)
                anc_sal_yn = int(len(set(genes) & set(anc_sal)) > 0)
                tofile.write(str(anc_all_yn) + '\t' + str(anc_salecocit_yn) + '\t' + str(anc_salcit_yn) + '\t'
                             + str(anc_kle_yn) + '\t' + str(anc_eco_yn) + '\t' + str(anc_sal_yn) + '\n')

                # prev_name = annot_dict[intersect]

tbl = pd.read_csv(outpath, sep='\t', header=0)
tbl = tbl.iloc[tbl.Gene_EGGNOG.str.lower().argsort()]
tbl.reset_index(inplace=True)
tbl.drop('index', axis=1, inplace=True)
tbl.to_csv(outpath, sep='\t', index=False, na_rep='-')
