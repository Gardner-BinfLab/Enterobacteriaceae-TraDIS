# USAGE: homclust.py [-options] <seq_fasta_file> <seqdb_fasta_file> <output_directory>

# Needs Python 3, HMMER 3, and mafft installed.

# Clusters homologous proteins in different species. The first input file is a protein fasta file from one
# species and the second input file is a file resulted from merging all fasta files (including the first input file)
# for all given species.

# For merging fasta files, one can use:
# cat <file1> <file2> ... <filen> > seqdb.fasta

#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from os import path, mkdir, system, listdir, remove, getcwd
from shutil import rmtree, copyfile, move
from Bio import SeqIO
from collections import defaultdict
from re import search
from argparse import ArgumentParser


def check_directory_name(dirname):
    if dirname.endswith('/'):
        dirname = dirname[:-1]
    return dirname


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


def removedir(dirname):
    # print ('Removing directory \'{0}\'..'.format(dirname))
    rmtree(dirname)


def get_dir_name_exten(filepath):
    filedir = path.dirname(filepath)
    filename_extension = path.basename(filepath)
    filename, fileexten = path.splitext(filename_extension)
    return filedir, filename, fileexten


def merge_dir_name_exten(filedir, filename, fileexten):
    filepath = filedir + '/' + filename + fileexten
    return filepath


def run_jackhmmer(seq1, seq2, output_dir, evalt, ievalt, msv, num_iter, bias, num_cpus):
    print ('Running Jackhmmer on \'{0}\' vs. \'{1}\'..'.format(seq1, seq2))
    seq1_dir, seq1_name, seq1_exten = get_dir_name_exten(seq1)
    seq2_dir, seq2_name, seq2_exten = get_dir_name_exten(seq2)
    output_name = '{0}-{1}-jackhmmer'.format(seq1_name, seq2_name)
    output_exten = '.domtblout'
    output_path = merge_dir_name_exten(output_dir, output_name, output_exten)
    if bias:
        bias_term = ''
    else:
        bias_term = ' --nobias'
    try:
        system('jackhmmer -N {0} -E {1} --incE {2} --F1 {3} --noali{4} --cpu {5} -o /dev/null --domtblout {6} {7} {8}'
               .format(num_iter, evalt, ievalt, msv, bias_term, num_cpus, output_path, seq1, seq2))
    except:
        raise SystemExit
    return output_path


# def read_fasta_ids(filepath):
#     # print ('Reading IDs from fasta file \'{0}\'..'.format(filepath))
#     with open(filepath, 'rU') as fasta_file:
#         record_list = list(SeqIO.parse(fasta_file, 'fasta'))
#     id_list = [record.id for record in record_list]
#     return id_list


def read_jackhmmer_ids(filepath):
    # print ('Reading IDs from Jackhmmer output \'{0}\'..'.format(filepath))
    id_list = []
    with open(filepath, 'r') as jackhmmer_file:
        for line in jackhmmer_file:
            cells = line.split()
            if not cells[0].startswith('#'):
                id_list.append(cells[0])
    return id_list


def list_orfans(sequence, jackhmmer_path, orfans_dir):
    print ('Listing orfans..')
    orfans_name = 'orfans'
    orfans_exten = '.fasta'
    orfans_path = merge_dir_name_exten(orfans_dir, orfans_name, orfans_exten)
    jackhmmer_ids = read_jackhmmer_ids(jackhmmer_path)
    with open(orfans_path, 'w') as orfans_file:
        with open(sequence) as seq_file:
            for record in SeqIO.parse(seq_file, 'fasta'):
                fastaid = record.id
                if not any(fastaid in s for s in jackhmmer_ids):
                    SeqIO.write(record, orfans_file, 'fasta')
    return orfans_path


def merge_jackhmmers(path_1, path_2, output_dir):
    output_name = 'merged-jackhmmer'
    output_exten = '.domtblout'
    output_path = merge_dir_name_exten(output_dir, output_name, output_exten)
    copyfile(path_1, output_path)
    with open(output_path, 'a') as to_file:
        with open(path_2) as from_file:
            for line in from_file:
                if not line.startswith('#'):
                    to_file.write(line)
    return output_path


def cluster_jackhmmer_output(input_dir, output_dir, fastadb, id_thresh):
    # FRAGMENTS CLASS STORES THE INFORMATION FOR EACH PREDICTED PROTEIN
    class Fragments(object):
        def __init__(self, name, domain_start, domain_end, protein_length):
            self.name = name
            self.domain_start = domain_start
            self.domain_end = domain_end
            self.domain_length = 0
            self.protein_length = protein_length
            self.write_to_file = ''

        def calculate_domain_length(self):
            self.domain_length += self.domain_end - self.domain_start + 1

        def make_write_to_file(self):
            self.write_to_file = self.write_to_file + \
                '{0}\t{1}\t{2}\n'.format(self.name, self.domain_start, self.domain_end)

    print ('Clustering the outputs of Jackhmmer..')
    sequence_dict = read_fasta_sequences(fastadb)
    query = '0'
    with open(input_dir, 'r') as jackhmmer_file:
        for line in jackhmmer_file:
            cells = line.split()
            if cells[0].startswith('#'):  # IGNORE COMMENT LINES
                continue
            if cells[3] != query:  # A NEW CLUSTER HAS STARTED
                if query != '0':  # THE CLUSTER IS NOT THE FIRST CLUSTER IN THE FILE
                    for key in frag_dict.keys():  # FOR EVERY MEMBER OF THE CLUSTER DO A QUALITY CONTROL AND IF THEY
                                                # PASS, ADD THEM TO THE FILE MADE FOR THAT SPECIAL CLUSTER
                        if frag_dict[key].domain_length >= id_thresh * frag_dict[key].protein_length and \
                           frag_dict[query].domain_length >= id_thresh * frag_dict[query].protein_length and \
                           frag_dict[key].domain_length >= id_thresh * frag_dict[query].domain_length and \
                           frag_dict[query].domain_length >= id_thresh * frag_dict[key].domain_length:
                            with open('{0}/{1}.txt'.format(output_dir, query), 'a') as cluster_file:
                                cluster_file.write(frag_dict[key].write_to_file)
                        # else:
                        #     with open('{0}/{1}.txt'.format(output_dir, key), 'a') as cluster_file:
                        #         cluster_file.write(frag_dict[key].write_to_file)
                query = cells[3]  # SET THE QUERY ID FOR THE NEW CLUSTER
                frag_dict = dict()  # EMPTY FRAG_DICT TO FILL IT WITH THE MEMBERS OF THE NEW CLUSTER
            domid = cells[0]
            domstart = int(cells[19])
            domend = int(cells[20])
            protlen = len(sequence_dict[domid])
            if domid not in frag_dict:  # MAKE AN OBJECT FOR DOMAIN_ID'S THAT ARE NOT REPEATED (DUE TO GAPS)
                frag_dict[domid] = Fragments(domid, domstart, domend, protlen)  # HAVE USED A DICTIONARY OF
                # CLASS OBJECTS TO BE ABLE TO NAME THE OBJECTS BY THEIR SEQUENCE ID AND ALSO TO CHECK EASILY IF A
                # SEQUENCE ID EXISTS IN THE DICTIONARY
            else:
                frag_dict[domid].domain_start = domstart
                frag_dict[domid].domain_end = domend
                frag_dict[domid].protein_length = protlen
            frag_dict[domid].calculate_domain_length()
            frag_dict[domid].make_write_to_file()
        for key in frag_dict.keys():  # QC AND MAKE A FILE FOR THE LAST CLUSTER
            if frag_dict[key].domain_length >= id_thresh * frag_dict[key].protein_length and \
               frag_dict[query].domain_length >= id_thresh * frag_dict[query].protein_length and \
               frag_dict[key].domain_length >= id_thresh * frag_dict[query].domain_length and \
               frag_dict[query].domain_length >= id_thresh * frag_dict[key].domain_length:
                with open('{0}/{1}.txt'.format(output_dir, query), 'a') as cluster_file:
                    cluster_file.write(frag_dict[key].write_to_file)
            # else:
            #     with open('{0}/{1}.txt'.format(output_dir, key), 'a') as cluster_file:
            #         cluster_file.write(frag_dict[key].write_to_file)


def sequence_level_intersect(clust1, clust2, threshold):
    overlaps = 0
    for counter1 in range(0, len(clust1.seqid)):
        for counter2 in range(0, len(clust2.seqid)):
            if clust1.seqid[counter1] == clust2.seqid[counter2]:
                intersect = float(min(clust1.end[counter1], clust2.end[counter2]) -
                                  max(clust1.start[counter1], clust2.start[counter2]) + 1) /\
                    float(min(clust1.end[counter1] - clust1.start[counter1],
                              clust2.end[counter2] - clust2.start[counter2]) + 1)
                if intersect >= threshold:
                    overlaps += 1
                    break
    final_intersect1 = float(overlaps) / len(clust1.seqid)
    overlaps = 0
    for counter2 in range(0, len(clust2.seqid)):
        for counter1 in range(0, len(clust1.seqid)):
            if clust1.seqid[counter1] == clust2.seqid[counter2]:
                intersect = float(min(clust1.end[counter1], clust2.end[counter2]) -
                                  max(clust1.start[counter1], clust2.start[counter2]) + 1) /\
                    float(min(clust1.end[counter1] - clust1.start[counter1],
                              clust2.end[counter2] - clust2.start[counter2]) + 1)
                if intersect >= threshold:
                    overlaps += 1
                    break
    final_intersect2 = float(overlaps) / len(clust2.seqid)
    final_intersect = max(final_intersect1, final_intersect2)
    return final_intersect


def merge_clusters(input_dir, output_dir, seq_id_thresh, clu_id_thresh):
    class Clusters(object):
        def __init__(self):
            self.seqid = []
            self.start = []
            self.end = []

        def clustappend(self, sequence_id, sequence_start, sequence_end):
            self.seqid.append(sequence_id)
            self.start.append(int(sequence_start))
            self.end.append(int(sequence_end))

    print ('Merging similar clusters..')
    adj_list = defaultdict(lambda: defaultdict(lambda: 0))   # ADJACENCY LIST FOR THE GRAPH (GRAPH NODES = FILES, GRAPH
    # EDGES = MORE THAN CLU_ID_THRESH% IDENTITY)
    sequence_occurrence = defaultdict(list)
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as cluster_file:
            for line in cluster_file:
                cells = line.split()
                sequence_occurrence[cells[0]].append(filename)
    for key in sequence_occurrence:
        for i in range(0, len(sequence_occurrence[key])-1):
            with open('{0}/{1}'.format(input_dir, sequence_occurrence[key][i])) as from_file:
                fromclust = Clusters()
                for line in from_file:
                    fromcells = line.split()
                    fromclust.clustappend(fromcells[0], fromcells[1], fromcells[2])
            for j in range(i+1, len(sequence_occurrence[key])):
                if sequence_occurrence[key][i] in adj_list and sequence_occurrence[key][j] in \
                        adj_list[sequence_occurrence[key][i]]:
                    continue
                with open('{0}/{1}'.format(input_dir, sequence_occurrence[key][j])) as to_file:
                    toclust = Clusters()
                    for line in to_file:
                        tocells = line.split()
                        toclust.clustappend(tocells[0], tocells[1], tocells[2])
                intersect = float(len(list(set(fromclust.seqid) & set(toclust.seqid))))
                if intersect >= clu_id_thresh * float(len(set(fromclust.seqid))) or intersect >= clu_id_thresh *\
                        float(len(set(toclust.seqid))):
                    seq_intersect = sequence_level_intersect(fromclust, toclust, seq_id_thresh)
                    if seq_intersect >= clu_id_thresh:
                        # ADD NODES THAT ARE CONNECTED BY AN EDGE TO THE ADJACENCY LIST
                        adj_list[sequence_occurrence[key][i]][sequence_occurrence[key][j]] = 1
                        adj_list[sequence_occurrence[key][j]][sequence_occurrence[key][i]] = 1
                    else:
                        adj_list[sequence_occurrence[key][i]][sequence_occurrence[key][j]] = 0
                        adj_list[sequence_occurrence[key][j]][sequence_occurrence[key][i]] = 0
                else:
                    adj_list[sequence_occurrence[key][i]][sequence_occurrence[key][j]] = 0
                    adj_list[sequence_occurrence[key][j]][sequence_occurrence[key][i]] = 0
    # FIND ALL NODES THAT ARE CONNECTED TOGETHER AND MERGE THEIR CORRESPONDING FILES
    visited_files = set()
    for filename in list_of_files:
        if filename not in visited_files:
            union = [filename]
            visited_files.add(filename)
            i = 0
            while i < len(union):
                for item in adj_list[union[i]]:
                    if adj_list[union[i]][item] and item not in visited_files:
                        visited_files.add(item)
                        union.append(item)
                i += 1
            with open('{0}/{1}'.format(output_dir, filename), 'w') as outfile:
                for unionname in union:
                    with open('{0}/{1}'.format(input_dir, unionname)) as infile:
                        for line in infile:
                            outfile.write(line)

def calculate_median(input_list):
    input_list = sorted(input_list)
    if len(input_list) < 1:
        return None
    if len(input_list) %2 == 1:
        return input_list[int((len(input_list)+1)/2)-1]
    else:
        return float(sum(input_list[int(len(input_list)/2)-1:int(len(input_list)/2)+1]))/2.0


def uniquify_clusters(input_dir, output_dir, threshold):
    print ('Uniquifying overlapping sequences within clusters..')
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        sequences = defaultdict(list)  # A DICTIONARY OF LIST OF TUPLES. EACH TUPLE CONTAINS THE START AND THE END OF
        # THE FRAGMENT. THE DICTIONARY CONTAINS ID OF THE SEQUENCE AND A LIST OF TUPLES FOR STARTS AND ENDS OF FRAGMENTS
        new_sequences = defaultdict(list)
        lengths = []
        with open('{0}/{1}'.format(input_dir, filename)) as from_file:
            for line in from_file:
                cells = line.split()
                sequences[cells[0]].append((int(cells[1]), int(cells[2])))
        for key in sequences:
            # SORT START-END TUPLES. THEN MERGE ALL OVERLAPPING TUPLES.
            sequences[key].sort()
            counter = 0
            while counter < len(sequences[key]):
                minimum = sequences[key][counter][0]
                maximum = sequences[key][counter][1]
                while counter < len(sequences[key]) - 1 and maximum > sequences[key][counter+1][0]:
                    maximum = max(maximum, sequences[key][counter+1][1])
                    counter += 1
                new_sequences[key].append((minimum, maximum))
                counter += 1
            summation1 = sum([pair[0] for pair in new_sequences[key]])
            summation2 = sum([pair[1] for pair in new_sequences[key]])
            lengths.append(summation2 - summation1 + len(new_sequences[key]))
        median = calculate_median(lengths)
        for key in new_sequences:
            summation1 = sum([pair[0] for pair in new_sequences[key]])
            summation2 = sum([pair[1] for pair in new_sequences[key]])
            if summation2 - summation1 + len(new_sequences[key]) >= median * threshold:
                with open('{0}/{1}'.format(output_dir, filename), 'a') as to_file:
                    for pair in new_sequences[key]:
                        to_file.write('{0}\t{1}\t{2}\n'.format(key, pair[0], pair[1]))



def collect_monsters(input_dir, seq_source, output_dir, output_thr):
    print('Separating monster families..')
    out_clusters = '{0}/clusters'.format(output_dir)
    makedir(out_clusters)
    out_fastas = '{0}/fastas'.format(output_dir)
    makedir(out_fastas)
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as from_file:
            num_lines = sum(1 for line in from_file)
        if num_lines >= output_thr:
            move('{0}/{1}'.format(input_dir, filename), '{0}/monster-{1}'.format(out_clusters, filename))
    make_fasta_clusters(out_clusters, out_fastas, seq_source, 'cluster')
    return out_clusters, out_fastas


def collect_singletons(input_dir, output_dir, output_thr, sequence_source):
    print('Collecting sequences in short clusters..')
    visited = set()
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        with open('{0}/{1}'.format(input_dir, filename)) as from_file:
            num_lines = sum(1 for line in from_file)
        if num_lines <= output_thr:
            remove('{0}/{1}'.format(input_dir, filename))
        else:
            with open('{0}/{1}'.format(input_dir, filename)) as from_file:
                for line in from_file:
                    cells = line.split()
                    visited.add(cells[0])
    singletons_f = merge_dir_name_exten(output_dir, 'singletons', '.fasta')
    with open(sequence_source) as source_file:
        with open(singletons_f, 'a') as singletons_file:
            for record in SeqIO.parse(source_file, 'fasta'):
                if record.id not in visited:
                    record.id = '{0}_{1}-{2}'.format(record.id, 1, len(record))
                    SeqIO.write(record, singletons_file, 'fasta')
    return singletons_f


def read_fasta_sequences(filepath):
    with open(filepath, 'rU') as fasta_file:
        fasta_dict = SeqIO.to_dict(SeqIO.parse(fasta_file, 'fasta'))
    return fasta_dict


def select_cluster_representative(input_dir, input_name):
    rep_file = '{0}/representative.fasta'.format(input_dir)
    flag = 0
    with open('{0}/{1}'.format(input_dir, input_name)) as from_file:
        with open(rep_file, 'w') as to_file:
            for line in from_file:
                if line.startswith('>'):
                    flag += 1
                if flag < 2:
                    to_file.write(line)
                else:
                    break
    with open(rep_file, 'rU') as fasta_file:
        sequence = list(SeqIO.parse(fasta_file, 'fasta'))
    rep_id = sequence[0].id
    return rep_file, rep_id


def copy_to_clusters(input_dir, output_dir):
    list_of_files = listdir(input_dir)
    for filename in list_of_files:
        move('{0}/{1}'.format(input_dir, filename), '{0}/{1}'.format(output_dir, filename))


def break_monster_families(output_dir, input_clusters, unique_input_clusters, input_fastas, eval_thresh, ieval_thresh,
                           msv_thresh, num_iter, bias, num_cpus, ovlp_thresh):
    list_of_files = listdir(input_fastas)
    for filename in list_of_files:
        fasta_path = '{0}/{1}'.format(input_fastas, filename)
        cluster_path = '{0}/{1}.txt'.format(input_clusters, filename[0:-6])
        sequences = read_fasta_sequences(fasta_path)
        counter = 0
        while sequences:
            with open(fasta_path, 'w') as new_cluster_file:
                for key in sequences:
                    garbage = SeqIO.write(sequences[key], new_cluster_file, 'fasta')
            rep_path, rep_id = select_cluster_representative(input_fastas, filename)
            jackhmmer_path = run_jackhmmer(rep_path, fasta_path, input_fastas, eval_thresh, ieval_thresh, msv_thresh,
                                           num_iter, bias, num_cpus)
            remove_list = read_jackhmmer_ids(jackhmmer_path)
            remove(jackhmmer_path)
            if rep_id in remove_list:
                with open('{0}/m-{1}-{2}.fasta'.format(input_fastas, filename[8:-6], counter), 'w') as fasta_file:
                    with open('{0}/m-{1}-{2}.txt'.format(input_clusters, filename[8:-6], counter), 'w') as cluster_file:
                        for key in remove_list:
                            if key in sequences:
                                fasta_file.write(sequences[key].format('fasta'))
                                cluster_line = search('(\w+)_(\d+)-(\d+)', sequences[key].id)
                                cluster_file.write('{0}\t{1}\t{2}\n'.format(cluster_line.group(1),
                                                                            cluster_line.group(2), cluster_line.group(3)
                                                                            ))
                                del sequences[key]
            else:
                with open('{0}/m-{1}-{2}.fasta'.format(input_fastas, filename[8:-6], counter), 'w') as fasta_file:
                    with open('{0}/m-{1}-{2}.txt'.format(input_clusters, filename[8:-6], counter), 'w') as cluster_file:
                        fasta_file.write(sequences[rep_id].format('fasta'))
                        cluster_line = search('(\w+)_(\d+)-(\d+)', sequences[rep_id].id)
                        cluster_file.write('{0}\t{1}\t{2}\n'.format(cluster_line.group(1), cluster_line.group(2),
                                                                    cluster_line.group(3)))
                        del sequences[rep_id]
            counter += 1
        remove(cluster_path)
        remove(fasta_path)
    uniquify_clusters(input_clusters, unique_input_clusters, ovlp_thresh)
    copy_to_clusters(unique_input_clusters, output_dir)


def code_clusters(input_dir, output_dir):
    class Proteins(object):
        def __init__(self, protein_id):
            self.protein_id = protein_id
            self.contig_id = 0

        def increase_contig_id(self):
            self.contig_id += 1
            return self.contig_id

    print ('Making EFam IDs for clusters..')
    list_of_files = listdir(input_dir)
    for clustid in range(0, len(list_of_files)):
        prot_dict = dict()
        protid = 0
        with open('{0}/{1}'.format(input_dir, list_of_files[clustid])) as infile:
            with open('{0}/EFam-{1:06d}.txt'.format(output_dir, clustid), 'w') as outfile:
                for line in infile:
                    cells = line.split()
                    prev_id = cells[0]
                    if prev_id not in prot_dict.keys():
                        prot_dict[prev_id] = Proteins(protid)
                        protid += 1
                    outfile.write('EFam-{0:06d}_{1:04d}_{2:02d}\t{3}\t{4}\t{5}\n'.format(clustid,
                                                                                         prot_dict[prev_id].protein_id,
                                                                                         prot_dict[prev_id].contig_id,
                                                                                         cells[0], cells[1], cells[2]))
                    prot_dict[prev_id].increase_contig_id()


def make_fasta_clusters(input_dir, output_dir, sequence_source, filetype):
    print ('Fetching fasta sequences for members of clusters..')
    # INDEX SEQUENCE_SOURCE
    try:
        system('esl-sfetch --index {0}'.format(sequence_source))
    except:
        raise SystemExit
    list_of_files = listdir(input_dir)
    for file_name in list_of_files:
        with open('{0}/{1}'.format(input_dir, file_name), 'r') as infile:
            for line in infile:
                cells = line.split()
                # FETCH THE SEQUENCE WITH GIVEN COORDINATES AND WRITE IT TO A FILE
                if filetype == 'EFam':
                    try:
                        system('esl-sfetch -c {0}-{1} {2} {3} | cat >> {4}/{5}.fasta'
                               .format(cells[2], cells[3], sequence_source, cells[1], output_dir, file_name[0:-4]))
                    except:
                        raise SystemExit
                elif filetype == 'cluster':
                    try:
                        system('esl-sfetch -c {0}-{1} {2} {3} | cat >> {4}/{5}.fasta'
                               .format(cells[1], cells[2], sequence_source, cells[0], output_dir, file_name[0:-4]))
                    except:
                        raise SystemExit
        try:
            system('sed -i \'s/\//_/g\' {0}/{1}.fasta'.format(output_dir, file_name[0:-4]))
        except:
            raise SystemExit


def multiple_sequence_alignment(input_dir, output_dir):
    print ('Performing multiple sequence alignments..')
    list_of_files = listdir(input_dir)
    output_exten = '.msa'
    for input_name in list_of_files:
        input_path = '{0}/{1}'.format(input_dir, input_name)
        output_name = input_name[0:-6]
        output_path = merge_dir_name_exten(output_dir, output_name, output_exten)
        count = 0
        for record in SeqIO.parse(input_path, 'fasta'):
            count += 1
            if count > 1:
                try:
                    system('mafft --text --quiet {0} > {1}'.format(input_path, output_path))
                except:
                    raise SystemExit
                break
        if count == 1:
            copyfile(input_path, output_path)


def build_hmms(input_dir, output_dir):
    print ('Building HMM profiles..')
    list_of_files = listdir(input_dir)
    output_exten = '.hmm'
    for input_name in list_of_files:
        input_path = '{0}/{1}'.format(input_dir, input_name)
        output_name = input_name[0:-4]
        output_path = merge_dir_name_exten(output_dir, output_name, output_exten)
        try:
            system('hmmbuild -o /dev/null {0} {1}'.format(output_path, input_path))
        except:
            raise SystemExit


def merge_files(input_dir, output_path):
    list_of_files = listdir(input_dir)
    with open(output_path, 'w') as to_file:
        for filename in list_of_files:
            with open('{0}/{1}'.format(input_dir, filename)) as from_file:
                for line in from_file:
                    to_file.write(line)


def run_hmmsearch(fastas, hmms, output_dir, num_cpus):
    output_path = merge_dir_name_exten(output_dir, 'hmmsearch', '.domtblout')
    try:
        system('hmmsearch -o /dev/null --noali --cpu {0} --domtblout {1} {2} {3}'.format(num_cpus, output_path, hmms,
                                                                                         fastas))
    except:
        raise SystemExit
    return output_path


def merge_hmmsearch_n_clusters(hmmsearch, out_dir, thresh):
    visited = set()
    with open(hmmsearch) as hmmsearch_file:
        for line in hmmsearch_file:
            cells = line.split()
            if not cells[0].startswith('#'):  # IGNORE COMMENT LINES
                query = cells[3]
                if not cells[0] in visited:
                    sequences = defaultdict(list)
                    lengths = []
                    with open('{0}/{1}.txt'.format(out_dir, query)) as cluster_file:
                        for new_line in cluster_file:
                            new_cells = new_line.split()
                            sequences[cells[0]].append((int(new_cells[1]), int(new_cells[2])))
                    for key in sequences:
                        summation1 = sum([pair[0] for pair in sequences[key]])
                        summation2 = sum([pair[1] for pair in sequences[key]])
                        lengths.append(summation2 - summation1 + len(sequences[key]))
                    median = calculate_median(lengths)
                    if int(cells[18]) - int(cells[17]) + 1 >= median * thresh:
                        with open('{0}/{1}.txt'.format(out_dir, query), 'a') as cluster_file:
                            seqid = search('(\w+)_\d+-\d+', cells[0])
                            cluster_file.write('{0}\t{1}\t{2}\n'.format(seqid.group(1), cells[17], cells[18]))
                        visited.add(cells[0])
    return list(visited)


def remove_hmmsearch_hits(input_dir, input_db, deletion_list):
    output_db = merge_dir_name_exten(input_dir, 'unique-singletons', '.fasta')
    with open(input_db, 'rU') as from_file:
        with open(output_db, 'w') as to_file:
            for record in SeqIO.parse(from_file, 'fasta'):
                if record.id not in deletion_list:
                    SeqIO.write(record, to_file, 'fasta')


def make_cluster_for_fastas(input_fasta, output_cluster):
    list_of_files = listdir(input_fasta)
    for filename in list_of_files:
        with open('{0}/{1}.txt'.format(output_cluster, filename[0:-6]), 'w') as to_file:
            with open('{0}/{1}'.format(input_fasta, filename)) as from_file:
                for record in SeqIO.parse(from_file, 'fasta'):
                    to_file.write('{0}\t{1}\t{2}\n'.format(record.id, 1, len(record)))


def merge_singletons(singleton_dir, fastadb, uniseq_dir, fasta_source, eval_thresh, ieval_thresh, msv_thresh, num_iter,
                     num_cpus, ovlp_thresh):
    uniseq_fastas = '{0}/uniseq-fastas'.format(singleton_dir)
    makedir(uniseq_fastas)
    uniseq_msas = '{0}/uniseq-msas'.format(singleton_dir)
    makedir(uniseq_msas)
    uniseq_hmms = '{0}/uniseq-hmms'.format(singleton_dir)
    makedir(uniseq_hmms)
    monster_fasta = '{0}/monster-fasta'.format(singleton_dir)
    makedir(monster_fasta)
    monster_cluster = '{0}/monster-cluster'.format(singleton_dir)
    makedir(monster_cluster)
    unique_monster_cluster = '{0}/unique-monster-cluster'.format(singleton_dir)
    makedir(unique_monster_cluster)
    make_fasta_clusters(uniseq_dir, uniseq_fastas, fasta_source, 'cluster')
    multiple_sequence_alignment(uniseq_fastas, uniseq_msas)
    build_hmms(uniseq_msas, uniseq_hmms)
    hmmsdb = merge_dir_name_exten(singleton_dir, 'uniseq-hmms', '.hmm')
    merge_files(uniseq_hmms, hmmsdb)
    hmmsearch_result = run_hmmsearch(fastadb, hmmsdb, singleton_dir, num_cpus)
    remove_list = merge_hmmsearch_n_clusters(hmmsearch_result, uniseq_dir, ovlp_thresh)
    remove_hmmsearch_hits(monster_fasta, fastadb, remove_list)
    make_cluster_for_fastas(monster_fasta, monster_cluster)
    break_monster_families(uniseq_dir, monster_cluster, unique_monster_cluster, monster_fasta, eval_thresh,
                           ieval_thresh, msv_thresh, num_iter, 0, num_cpus, ovlp_thresh)

# READ INPUTS AND SET PARAMETERS
et = 10**-10
iet = 10**-3
mt = 10**-3
it = 80
ni = 5
parser = ArgumentParser(description='Clusters homologous proteins. Needs python 3 or higher and HMMER 3.')
parser.add_argument('-e', '--eva', help='E-value threshold for Jackhmmer', default=et)
parser.add_argument('-i', '--iev', help='Inclusion E-value threshold for Jackhmmer', default=iet)
parser.add_argument('-m', '--msv', help='MSV threshold for Jackhmmer', default=mt)
parser.add_argument('-o', '--ovlp', help='Sequence overlap for merging clusters', default=it)
parser.add_argument('-c', '--cov', help='Coverage for sequences within a cluster', default=it)
parser.add_argument('-n', '--num', help='Number of Jackhmmer iterations', default=ni)
parser.add_argument('-p', '--cpu', help='Number of CPUs', default=1)
parser.add_argument('-j1', '--jck1', help='The results of first Jackhmmer (representative sequences file vs. seqdb)')
parser.add_argument('-j2', '--jck2', help='The results of second Jackhmmer (orfans vs. seqdb)')
parser.add_argument('-j', '--jck', help='The results of both Jackhmmer runs (merged-jackhmmer)')
parser.add_argument('--nojck2', help='Does not run the second jackhmmer', action="store_true")
parser.add_argument('-ns', '--nspec', help='Number of species under study')
parser.add_argument('-s', '--sing', help='Threshold for defining singletons')
parser.add_argument('seq', help='Fasta file for protein sequences in one species')
parser.add_argument('seqdb', help='Fasta file for protein sequences in all species')
parser.add_argument('outdir', help='Directory for outputs')
args = parser.parse_args()
eval_threshold = float(args.eva)
inclusion_eval_threshold = float(args.iev)
msv_threshold = float(args.msv)
num_iterations = int(args.num)
seq_identity_threshold = float(args.ovlp)/100
clu_identity_threshold = float(args.cov)/100
cpu = int(args.cpu)
jackhmmer_1 = args.jck1
jackhmmer_2 = args.jck2
jackhmmer = args.jck
number_of_species = int(args.nspec)
monster_thresh = 3 * number_of_species
singleton_thresh = int(args.sing)
seq = args.seq
seqdb = args.seqdb
outdir = args.outdir
outdir = check_directory_name(outdir)
# MAKE DIRECTORIES FOR OUTPUTS
print ('Making directories for results..')
makedir(outdir)
temp = '{0}/temp'.format(outdir)
makedir(temp)
clusters = '{0}/clusters'.format(outdir)
makedir(clusters)
uniseq_clusters = '{0}/unique-sequence-clusters'.format(outdir)
makedir(uniseq_clusters)
efamclusters = '{0}/EFam-clusters'.format(outdir)
makedir(efamclusters)
efamfastas = '{0}/EFam-fastas'.format(outdir)
makedir(efamfastas)
efammsas = '{0}/EFam-MSAs'.format(outdir)
makedir(efammsas)
efamhmms = '{0}/EFam-HMMs'.format(outdir)
makedir(efamhmms)
monsters = '{0}/monster-families'.format(outdir)
makedir(monsters)
singletons = '{0}/singletons'.format(outdir)
makedir(singletons)
if jackhmmer is None or not path.exists(jackhmmer):
    if jackhmmer_1 is None or not path.exists(jackhmmer_1):
        # RUN JACKHMMER ON A REPRESENTATIVE SPECIES vs. ALL
        jackhmmer_1 = run_jackhmmer(seq, seqdb, outdir, eval_threshold, inclusion_eval_threshold, msv_threshold,
                                    num_iterations, 1, cpu)
    if (jackhmmer_2 is None or not path.exists(jackhmmer_2)) and not args.nojck2:
        # FIND ALL ORFANS
        orfans = list_orfans(seqdb, jackhmmer_1, outdir)
        # RUN JACKHMMER ON ORFANS vs. ALL
        jackhmmer_2 = run_jackhmmer(orfans, seqdb, outdir, eval_threshold, inclusion_eval_threshold, msv_threshold,
                                    num_iterations, 1, cpu)
    if not args.nojck2:
        jackhmmer = merge_jackhmmers(jackhmmer_1, jackhmmer_2, outdir)
    else:
        jackhmmer = jackhmmer_1
# CLUSTER THE OUTPUT OF JACKHMMER BASED ON QUERY SEQUENCE AND REMOVE PROTEINS THAT DO NOT PASS QUALITY CONTROL
cluster_jackhmmer_output(jackhmmer, temp, seqdb, seq_identity_threshold)
## MERGE CLUSTERS THAT HAVE MORE THAN IDENTITY_THRESHOLD% IDENTITY
# merge_clusters(temp, clusters, seq_identity_threshold, clu_identity_threshold)
## MERGE OVERLAPPING SEQUENCES WITHIN A CLUSTER
# uniquify_clusters(clusters, uniseq_clusters, seq_identity_threshold)
# COLLECT MONSTER FAMILIES
monsters_cluster, monsters_fasta = collect_monsters(temp, seqdb, monsters, monster_thresh)
# BREAK MONSTER FAMILIES TO SMALLER CLUSTERS
unique_monsters_clusters = '{0}/unique-monsters-clusters'.format(monsters)
makedir(unique_monsters_clusters)
break_monster_families(temp, monsters_cluster, unique_monsters_clusters, monsters_fasta, eval_threshold**2,
                       inclusion_eval_threshold**2, msv_threshold**2, 1, 1, cpu, seq_identity_threshold)
# COLLECT THE SEQUENCES IN SHORT CLUSTERS
singletons_fasta = collect_singletons(temp, singletons, singleton_thresh, seqdb)
# CHECK IF SINGLETONS BELONG TO OUR DEFINED CLUSTERS AND ADD THEM TO THEIR CORRESPONDING CLUSTERS
merge_singletons(singletons, singletons_fasta, temp, seqdb, eval_threshold**(1/2),
                 inclusion_eval_threshold**(1/2), msv_threshold**(1/2), num_iterations, cpu, seq_identity_threshold)
# MERGE CLUSTERS THAT HAVE MORE THAN IDENTITY_THRESHOLD% IDENTITY
merge_clusters(temp, clusters, seq_identity_threshold, clu_identity_threshold)
# MERGE OVERLAPPING SEQUENCES WITHIN A CLUSTER
uniquify_clusters(clusters, uniseq_clusters, seq_identity_threshold)
# RENAME CLUSTERS BY EFAM IDS
code_clusters(uniseq_clusters, efamclusters)
# FETCH FASTA SEQUENCES FOR MEMBERS OF THE CLUSTERS
make_fasta_clusters(efamclusters, efamfastas, seqdb, 'EFam')
# DO MULTIPLE SEQUENCE ALIGNMENTS ON CLUSTERS
multiple_sequence_alignment(efamfastas, efammsas)
# MAKE PROFILE HMMS FOR CLUSTERS
build_hmms(efammsas, efamhmms)
print ('Removing temporary directories..')
removedir(temp)
removedir(clusters)
removedir(singletons)
removedir(monsters)
removedir(uniseq_clusters)
print ('Done!')
#print ('Directory \'{0}\' \t contains homologous clusters'.format(uniseq_clusters))
print ('Directory \'{0}\' \t contains homologous clusters with EFam IDs'.format(efamclusters))
print ('Directory \'{0}\' \t contains fasta sequences for homologous clusters'.format(efamfastas))
print ('Directory \'{0}\' \t contains Multiple sequence alignments for clusters'.format(efammsas))
print ('Directory \'{0}\' \t contains HMM profiles for clusters'.format(efamhmms))
