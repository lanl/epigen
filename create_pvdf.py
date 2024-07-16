#input: hg38_chrID_200data files that are generated through run_get200bpvalues_for_chromosome_hg38_script
#input: genome_df38.csv contains Target variable that we used as indicies for pvdf_chr_id.pkl dataframes (any other variables or just numbers can be used, this is for tracking purposes as samples arranged in the same order as in genome_df38.csv file)

import pandas as pd
import numpy as np
import h5py
import sys


df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
for chr_id in range(int(float(sys.argv[1])), int(float(sys.argv[2]))):
    if chr_id == 23: #this is chromosome X
        chr_id = 'X'
    print(chr_id)
    filename = "./results38/hg38_chr" + str(chr_id) + "_200data"
    with h5py.File(filename + '.h5', 'r') as hdf:
        data_group = hdf.get('data')
        # Retrieve each component of the DataFrame
        axis0 = np.array(data_group.get('axis0')).astype(str)
        axis1 = np.array(data_group.get('axis1')).astype(str)
        block0_items = np.array(data_group.get('block0_items')).astype(str)
        block0_values = np.array(data_group.get('block0_values'))
    df = pd.DataFrame(data=block0_values, index=axis1, columns=axis0)
    df.index = df38['Target']
    c = pd.DataFrame(df.max(axis='rows'))
    c.columns = ['max']
    c = c.loc[c['max'] >= 1.3, 'max'] #original data are signal p-values (-log10(regular p-value)), here we'll only extract columns where regular p-value <= 0.05
    pv_df = 10**(-df.iloc[:, list(map(int, c.index))])
    pv_df.index = range(0, len(pv_df))
    pv_df.to_pickle("./results38/pvdf_" + str(chr_id) + ".pkl")

