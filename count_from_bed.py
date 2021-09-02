input_files = ['LB1.bed','LB2.bed','Bile1.bed','Bile2.bed']
output_file = 'count.tab'
counts = dict()
with open(output_file,'w') as new:
    new.write('Gene\t')
    for sample in len(input_files):
        input_file = input_files[sample]
        new.write('%s\t'%(input_files))
        with open(input_files,'r') as old:
            prev_gene = ''
            count =0
            for line in old:
                read = line.strip().split()
                if read[0] == prev_gene:
                    count +=1
                else:
                    if prev_gene in counts:
                        counts[prev_gene][sample] = count
                    elif prev_gene != '':
                        counts[prev_gene] = [0] * len(input_files)
                        counts[prev_gene][sample] = count
                    prev_gene = read[0]
                    count = 1
    for gene in counts:
        new.write('%s\t'%(gene))
        for count in counts[gene]:
            new.write('%i\t'%(count))
        new.write('\n')




