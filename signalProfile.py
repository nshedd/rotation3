## inputs:
##      BW_files_list: list of /path/to/file.bigwig
##      loci_list: list of chr1:10,000,000-10,001,000
##      output_file: /path/to/output.tsv

import numpy as np
import pandas as pd
import pyBigWig
import sys

BW_files_list = sys.argv[1] 
loci_file = sys.argv[2] 
output_name = sys.argv[3]

BW_files = [] 

with open(BW_files_list, 'r') as file:
    for line in file:
        BW_files.extend(line.split()) 

loci_list = []

with open(loci_file, 'r') as file:
    for line in file:
        loci_list.extend(line.split()) 

dataframe = pd.DataFrame(data={'File': BW_files})

for locus in loci_list:

    chrom = locus.split(':')[0]
    print(chrom)
    location = locus.split(':')[1]
    loc1 = int((location.split('-')[0]).replace(',', ""))
    loc2 = int(location.split('-')[1].replace(',', ""))
    print(loc1)
    print(loc2)

    signal_values = []

    for file in BW_files:
        print(file)
        with pyBigWig.open(file) as b:
            values = [ x if not np.isnan(x) else 0 for x in b.values(chrom, loc1, loc2) ]
        mean_signal = sum(values) / float(len(values))
        signal_values.append(mean_signal)

    dataframe[locus] = signal_values
dataframe.to_csv(output_name, sep='\t', index=False)
