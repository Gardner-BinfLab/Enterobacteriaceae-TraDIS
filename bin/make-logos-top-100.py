#! /bin/sh
""":"
exec python3 $0 ${1+"$@"}
"""

from os import listdir, path, remove
import numpy as np

genomespath = '../data/fasta-genome/chromosome'
plotspath = '../data/plot-files/chromosome'
genomextension = 'fa'
logospath = '../results/logos/100logos.txt'
complement = {'a': 't', 'c': 'g', 'g': 'c', 't': 'a', 'u': 'a', 'm': 'k', 'k': 'm', 'r': 'y',
              'y': 'r', 'w': 'w', 's': 's', 'v': 'b', 'b': 'v', 'h': 'd', 'd': 'h', 'n': 'n'}

if path.isfile(logospath):
    input_var = 'i'
    while not (input_var == 'y' or input_var == 'Y' or input_var == 'n' or input_var == 'N'):
        input_var = input('File {0} already exists. Replace? [Y/n] '.format(logospath))
        if input_var == 'Y' or input_var == 'y':
            remove(logospath)
        else:
            raise SystemExit

list_of_files = listdir(plotspath)
for filename in list_of_files:
    plots_f = []
    plots_b = []
    sequences = ''
    name, extension = path.splitext(filename)
    with open('{0}/{1}'.format(plotspath, filename)) as plotfile:
        for line in plotfile:
            cells = line.split()
            plots_f.append(int(cells[0]))
            plots_b.append(int(cells[1]))
    with open('{0}/{1}.{2}'.format(genomespath, name, genomextension)) as genomefile:
        for line in genomefile:
            if not line.startswith('>'):
                if line.endswith('\n'):
                    line = line[:-1]
                sequences += line
    plots = np.array(plots_f + plots_b)
    top_indices = np.argpartition(plots, -100)[-100:] # returns the indices of 100 maximums in plots
    with open (logospath, 'a') as logosfile:
        # for i in range(0, len(plots_f)):
        #     if plots_f[i] > 0 or plots_b[i] > 0:
        #         if 10 <= i < len(sequences)-10:
        #             seq_to_write = sequences[i-10:i+11]
        #         elif i < 10:
        #             seq_to_write = sequences[(i-10)%len(sequences):len(sequences)]+sequences[0:i+11]
        #         else:
        #             seq_to_write = sequences[i-10:len(sequences)]+sequences[0:(i+11)%len(sequences)]
        #         if plots_f[i] > 0:
        #             logosfile.write('>{0}_{1} forward\n'.format(name, i+1))
        #             logosfile.write('{0}\n'.format(seq_to_write))
        #         if plots_b[i] > 0:
        #             bases = list(seq_to_write.lower())
        #             letters = [complement[elements] for elements in bases]
        #             seq_to_write = ''.join(letters)
        #             logosfile.write('>{0}_{1} backward\n'.format(name, i+1))
        #             logosfile.write('{0}\n'.format(seq_to_write[::-1]))

        # for item in top_indices:
        #     if plots[item] > 0:
        #         i = item%len(plots_f)
        #         if 10 <= i <= len(sequences)-19:
        #             seq_to_write = sequences[i-10:i+19]
        #         elif i < 10:
        #             seq_to_write = sequences[(i-10)%len(sequences):len(sequences)]+sequences[0:i+19]
        #         else:
        #             seq_to_write = sequences[i-10:len(sequences)]+sequences[0:(i+19)%len(sequences)]
        #         if item < len(plots_f):
        #             logosfile.write('>{0}_{1} forward\n'.format(name, i+1))
        #             logosfile.write('{0}\n'.format(seq_to_write))
        #         else:
        #             bases = list(seq_to_write.lower())
        #             letters = [complement[elements] for elements in bases]
        #             seq_to_write = ''.join(letters)
        #             logosfile.write('>{0}_{1} backward\n'.format(name, i+1))
        #             logosfile.write('{0}\n'.format(seq_to_write[::-1]))
        for item in top_indices:
            if plots[item] > 0:
                i = item % len(plots_f)
                if item < len(plots_f):
                    if 10 <= i <= len(sequences) - 19:
                        seq_to_write = sequences[i - 10:i + 19]
                    elif i < 10:
                        seq_to_write = sequences[(i - 10) % len(sequences):len(sequences)] + sequences[0:i + 19]
                    else:
                        seq_to_write = sequences[i - 10:len(sequences)] + sequences[0:(i + 19) % len(sequences)]
                    logosfile.write('>{0}_{1} forward\n'.format(name, i + 1))
                    logosfile.write('{0}\n'.format(seq_to_write))
                else:
                    if 18 <= i <= len(sequences) - 11:
                        seq_to_write = sequences[i - 18:i + 11]
                    elif i < 18:
                        seq_to_write = sequences[(i - 18) % len(sequences):len(sequences)] + sequences[0:i + 11]
                    else:
                        seq_to_write = sequences[i - 18:len(sequences)] + sequences[0:(i + 11) % len(sequences)]
                    bases = list(seq_to_write.lower())
                    letters = [complement[elements] for elements in bases]
                    seq_to_write = ''.join(letters)
                    logosfile.write('>{0}_{1} backward\n'.format(name, i + 1))
                    logosfile.write('{0}\n'.format(seq_to_write[::-1]))