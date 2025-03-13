# get_values_for_chromosome.py
# This file processes each .bigWig file in the input index range.
# For each file:
#     get a df that contains 200bp values from bigWig file for a specified chromosome (script argument)
# note that chr 23 here refers to chr X and there is no chr Y in the data (ENCODE files)
# this file requires `/data/hg38length.txt` which contains the length of each chromosome for hg38 assembly
#
#
# how to run: takes two arguments: i=sys.argv[1] and j=sys.argv[2]
# It then downloads ith-jth range of files from the genome_df38.csv
#
# OUTPUT: note that computing all files (1632 total) in a single
# job might take too long, it is better to save results for every ~200 files,
# for example the output file names for chromosome 6 (chr_i=6) might look like these:
# hg38_chr6_200_0_200.h5, hg38_chr6_200_200_400.h5, ..., hg38_chr6_200_1400_1632.h5
# that is, hg38_chr6_{bin_size}_{start_sample}_{end_sample}
# Depending on the chromosome, the set of files will take approximately 2GB to 16GB per chromosome.
#
# The chunk size of 200 samples is hardcoded to later files!
# If you want to change this variable you will need to
# change the downstream scripts to read using different filenames.
# However, do not confuse this with the 200 base pair bin length.
#

import pandas as pd
import numpy as np
import pyBigWig
import sys

# Get the metadata for the samples
df = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = df.loc[df['Genome'] == 'hg38']
df38 = df38.reset_index()
del df38['index']

# Get the meta for the chromosome
chr_length = pd.read_csv("./data/hg38length.txt", sep = '\t', header = None)
chr_length.columns = ['Chromosome', 'Total length (bp)', 'Nm', 'Nmm']
chr_id = int(float(sys.argv[1])) # chromosome ID

start_p = 0
bin_length = 200

if chr_id == 23:  # This is chromosome X
    chr_id = 'X'

end_p = (int(chr_length.loc[chr_length['Chromosome'] == str(chr_id), "Total length (bp)"]) // bin_length) * bin_length

# loop over sample files
chr_vals = []
for file_num in range(int(float(sys.argv[2])), int(float(sys.argv[3]))):
    temp = []
    i = 0
    file_name = str(df38['Accession'][file_num])
    bw = pyBigWig.open(str('./hg38data/') + str(file_name) + str('.bigWig'))

    # Loop over bins in the file
    print(file_num, file_name)
    while i < end_p: 
        temp.append(np.median(bw.values('chr' + str(chr_id), i, i + bin_length)))
        i = i + bin_length
    chr_vals.append(temp)

# Reassemble all counts as a dataframe
test_data = pd.DataFrame(chr_vals)
filename = './results38/hg38_chr' + str(chr_id) + '_200_' + str(int(float(sys.argv[2]))) + '_' + str(int(float(sys.argv[3]))) + '.h5'
test_data.to_hdf(filename, key='data', mode='w')
