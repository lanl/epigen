
#Section 1: the first step is to combine all h5 data file chunks for each chromosome into the one h5 file for each chromosome
import pandas as pd
import requests
import h5py
import json
import os
import numpy as np

for chr_i in range(1, 24):
	print(chr_i)
	if chr_i == 23:
	    chr_i = 'X'
	data1 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_0_200.h5", key='data')
	data2 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_200_400.h5", key='data')
	data3 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_400_600.h5", key='data')
	data4 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_600_800.h5", key='data')
	data5 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_800_1000.h5", key='data')
	data6 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1000_1200.h5", key='data')
	data7 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1200_1400.h5", key='data')
	data8 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1500_1632.h5", key='data')
	tdf = pd.concat([data1, data2], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data3], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data4], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data5], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data6], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data7], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data8], axis = 0).reset_index(drop=True)
	tdf.shape
	filename = "./hg38_chr" + str(chr_i) + "_200data.h5" #saving combined h5 files for each chr
	tdf.to_hdf(filename, key='data', mode='w')


#Section 2 needs GPU: the second step is to use cudf (to use GPU for pandas dfs, see RAPIDS library) to calculate correlation matrices
#because the data may contain ~1 million columns it will take hours to compute pairwise correlation of 1632 samples with so much columns
#GPU usage is needed
#module load python/3.10-anaconda-2023.03
#module load cuda/11.5 #make sure cuda is enabled
#conda activate rapids-23.06 #install RAPIDS library (large)

import cudf
import numpy as np
import pandas as pd
import h5py
from scipy.spatial.distance import pdist, squareform
import time
from datetime import timedelta


for chr_i in range(1, 24):
	print(chr_i)
	if chr_i == 23:
	    chr_i = 'X'
    path = "./hg38_chr"
    filename = "./results38/hg38_chr" + str(chr_i) + "_200data"
    start_time = time.monotonic()
    correlation_matrix = 0
    with h5py.File(filename + '.h5', 'r') as hdf:
        data_group = hdf.get('data')
        # Retrieve each component of the DataFrame
        axis0 = np.array(data_group.get('axis0')).astype(str)
        axis1 = np.array(data_group.get('axis1')).astype(str)
        block0_items = np.array(data_group.get('block0_items')).astype(str)
        block0_values = np.array(data_group.get('block0_values'))
    df_read = pd.DataFrame(data=block0_values, index=axis1, columns=axis0)
    # Reset index before transposing
    df_read.index = range(0,1632) #1632 sample files in the hg38-aligned dataset
    df_read = df_read.fillna(0) #rare but can happen
    df_transposed = df_read.transpose()
    df_cudf = cudf.DataFrame.from_pandas(df_transposed)
    correlation_matrix = df_cudf.corr().transpose()
    correlation_matrix = 1 - correlation_matrix
    #print(correlation_matrix)
    correlation_matrix.to_csv(filename + 'correlation.h5') #saving correlation files
    #pd.read_csv(filename + 'correlation.h5', index_col=0)
    #print(chr_i, " ", df_read.shape)
    end_time = time.monotonic()
    print(timedelta(seconds=end_time - start_time))

