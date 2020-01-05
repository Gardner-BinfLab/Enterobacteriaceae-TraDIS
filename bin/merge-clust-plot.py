#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from sys import argv
from os import listdir, mkdir, path
from shutil import rmtree
from re import match


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


def add_insertion_index(in_clusters, in_iis, k12essentials, output):
    iis_dict = {}
    list_of_files = listdir(in_iis)
    for filename in list_of_files:
        with open(in_iis+'/'+filename) as iis_file:
            for line in iis_file:
                cells = line.split()
                iis_dict[cells[0]] = cells[1]
    list_of_files = listdir(in_clusters)
    for filename in list_of_files:
        with open('{0}/{1}'.format(in_clusters, filename)) as from_file:
            with open('{0}/{1}'.format(output, filename), 'w') as to_file:
                for line in from_file:
                    cells = line.split()
                    match_result = match('b\d\d\d\d', cells[1])
                    for item in cells:
                        to_file.write(item+'\t')
                    if cells[1] in iis_dict.keys():
                        to_file.write(iis_dict[cells[1]]+'\n')
                    elif match_result is not None:
                        if cells[1] in k12essentials:
                            to_file.write('0\n')
                        else:
                            to_file.write('7\n')
                    else:
                        to_file.write('7\n')

indir = str(argv[1])
iifile = str(argv[2])
k12file = str(argv[3])
outdir = str(argv[4])
makedir(outdir)
k12list = list()
with open (k12file, 'r') as infile:
    for line in infile:
        k12list.append(line.rstrip())
add_insertion_index(indir, iifile, k12list, outdir)