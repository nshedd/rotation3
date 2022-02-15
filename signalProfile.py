## run as: python3 signalProfile.py /data/rusers/sheddn/TR-ATACseq/all_files.txt /data/zusers/fankaili/ccre/ubi-rDHS_V3/GRCh38_ubi-rDHS_V3.txt /data/rusers/sheddn/ubi_rDHS/GRCh38_ubi-rDHS_V3_signals.tsv

## inputs:
##      BW_files_list: list of /path/to/file.bigwig
##      loci_list: list of chr1 10000000    10001000
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

loci_list = pd.read_csv(loci_file, header=None, sep='\t')
print(loci_list)

loci_table = loci_list[0] + ":" + loci_list[1].astype(str) + "-" + loci_list[2].astype(str)
loci_table = loci_table.tolist()
print(loci_table)

dataframe = pd.DataFrame(data={'Locus': loci_table})

for file in BW_files:
        print(file)

        signal_values = []

        for locus in (range(0,len(loci_list))):
            chrom = loci_list.iat[locus,0]
            loc1 = int(loci_list.iat[locus,1])
            loc2 = int(loci_list.iat[locus,2])

            with pyBigWig.open(file) as b:
                values = [ x if not np.isnan(x) else 0 for x in b.values(chrom, loc1, loc2) ]
            mean_signal = sum(values) / float(len(values))
            signal_values.append(mean_signal)

        dataframe[file] = signal_values

dataframe.to_csv(output_name, sep='\t', index=False)


