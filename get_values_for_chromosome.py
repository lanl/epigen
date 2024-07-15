#get a df that contains 200bp values from bigWig file for a specified chromosome (script argument)
#note that chr 23 here refers to chr X and there is no chr Y in the data (ENCODE files)
#hg38length.txt contains the length of each chromosome for hg38 assembly

import pandas as pd
import requests
import h5py
import json
import os
import numpy as np
import pyBigWig
import urllib.request 
import sys

df = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = df.loc[df['Genome'] == 'hg38']
df38 = df38.reset_index()
del df38['index']

chr_length = pd.read_csv("./data/hg38length.txt", sep = '\t', header = None)
chr_length.columns = ['Chromosome', 'Total length (bp)', 'Nm', 'Nmm']
chr_id = int(float(sys.argv[1]))
start_p = 0
bin_length = 200
if chr_id == 23: #this is chromosome X
    end_p = (int(chr_length.loc[chr_length['Chromosome'] == 'X', "Total length (bp)"])// bin_length) * bin_length
    chr_id = 'X'
else:
    end_p = (int(chr_length.loc[chr_length['Chromosome'] == str(chr_id), "Total length (bp)"])// bin_length) * bin_length


chr_vals = []
for file_num in range(int(float(sys.argv[2])), int(float(sys.argv[3]))):
    temp = []
    i = 0
    file_name = str(df38['Accession'][file_num])
    bw = pyBigWig.open(str('./hg38data/') + str(file_name) + str('.bigWig'))
    print(file_num, file_name)
    while i < end_p: 
        temp.append(np.median(bw.values('chr' + str(chr_id), i, i + bin_length)))
        i = i + bin_length
    chr_vals.append(temp)
    


test_data = pd.DataFrame(chr_vals)
filename = './results38/hg38_chr' + str(chr_id) + '_200_' + str(int(float(sys.argv[2]))) + '_' + str(int(float(sys.argv[3]))) + '.h5'
test_data.to_hdf(filename, key='data', mode='w')
#OUTPUT: note that computing all files (1632 total) in a single job might take too long, it is better to save results for every ~200 files, for example the output file names
#for chromosome 6 (chr_i=6) might look like these: 
#hg38_chr6_200_0_200.h5, hg38_chr6_200_200_400.h5, ..., hg38_chr6_200_1400_1632.h5
