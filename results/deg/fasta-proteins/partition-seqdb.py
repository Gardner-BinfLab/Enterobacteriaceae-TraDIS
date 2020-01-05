# Partitions seqdb because it's too large for submitting to eggnog-mapper and eggnog-mapper will run diamond instead of hmmer on it.
fcounter = 1
lcounter = 0
seqdb = 'seqdb.fasta'
with open(seqdb, 'r') as fromfile:
    for line in fromfile:
        if line.startswith('>'):
            lcounter += 1
            if lcounter > 4800:
                fcounter += 1
                lcounter = 1
        with open('seqdb-' + str(fcounter)+ '.fasta', 'a') as tofile:
            tofile.write(line)