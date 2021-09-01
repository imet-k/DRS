#identification TTS by intersection
tss = list()
def clustering(a):
    final = list() # clustering can produce more than TTS
    prev = a[0]
    curr = [a[0]]
    for i in a[1:]:
        if i - prev < 10: #as long as the distance between reads are less than 10 bp, take them in the same cluster
            prev = i
            curr.append(i)
        else:
            if len(curr) > 0.2*len(a): #if the current cluster count is significant
                final.append((max(set(curr), key=curr.count), len(curr))) #then report the most frequent site as TSS
            prev = i
            curr = [i] #start a new cluster
    if len(curr) > 0.2 * len(a):
        final.append((max(set(curr), key=curr.count),len(curr)))
    return(final)
#input file is obtained by bedtools intesect -a genes_file.bed -b reads_file.bed -wo
input_file = 'intersect_genemodel_n104.bed'
output_file = 'TTS_n104.tab'
with open(input_file,'r') as bed, open(output_file,'w') as new:
    prev_gene = ''
    current = list()
    for line in bed:
        read = line.strip().split()
        if read[3] == prev_gene:
            if read[5] == read[11] and (min(int(read[2]),int(read[8]))-max(int(read[1]),int(read[7])))/(int(read[2])-int(read[1])) > 0.9:
                if read[11] == '+' and int(read[8]) > int(read[2]):
                    current.append((int(read[8])))
                elif read[11] == '-' and int(read[1]) > int(read[7]):
                    current.append((int(read[7])))
        else:
            if len(current) > 0 and prev_gene !='':
                current.sort()
                for site in clustering(current):
                    new.write('%s\t%i\t%i\n'%(prev_gene,site[0],site[1]))
            prev_gene = read[3]
            current = []
    if len(current) > 0:
        current.sort()
        for site in clustering(current):
            new.write('%s\t%i\t%i\n' % (prev_gene, site[0], site[1]))