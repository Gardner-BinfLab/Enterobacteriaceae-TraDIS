from os import listdir

k12 = '../results/ecogene-k12.txt'
clusters = '../results/hieranoid/clusters.txt'
genes = '../results/insertion-indices/gamma/'
output = '../results/ecogenecounterparts/'
list_of_files = listdir(genes)
genenames = []

for filename in list_of_files:
    print('#################'+filename[:-4])
    with open(k12, 'r') as k12file:
        for gene in k12file:
            gene = gene.rstrip()
            with open(clusters, 'r') as clusterfile:
                for line in clusterfile:
                    cells = line.split()
                    # counter = 0
                    if gene in cells:
                        flag = 0
                        for item in cells:
                            if item.startswith(filename[:-4]):
                                # counter += 1
                                # if counter > 1:
                                # print('heyy')
                                flag += 1
                                genenames.append(item)
                        if not flag:
                            print(gene)

    allgenes = []
    with open(genes + filename, 'r') as fromfile:
        for line in fromfile:
            cells = line.split()
            allgenes.append(cells[0])

    essentiality = []
    for item in allgenes:
        if item not in genenames:
            essentiality.append('0')
        else:
            essentiality.append('1')

    with open(output + filename, 'w') as tofile:
        for (item1, item2) in zip(allgenes, essentiality):
            tofile.write(item1 + '\t' + item2 + '\n')