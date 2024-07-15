#input: hg38_chrID_200data files that are generated 
import pandas as pd
import numpy as np
import h5py
import sys

start_p = 0
divisor = 200
binwidth = 200
setbp = 500 
#all
df38 = pd.read_csv("./genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
for chr_id in range(int(float(sys.argv[1])), int(float(sys.argv[2]))):
    print(chr_id)
    filename = "./results38/fin/hg38_chr" + str(chr_id) + "_200data"
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
    c = c.loc[c['max'] >= 1.3, 'max'] 
    pv_df = 10**(-df.iloc[:, list(map(int, c.index))])
    pv_df.index = range(0, len(pv_df))
    pv_df.to_pickle("./results38/fin/pvdf_" + str(chr_id) + ".pkl")

