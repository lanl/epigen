
#First step to get h5 with data for each chr

module load python/3.8-anaconda-2020.07
source activate
conda activate myenv
python


import pandas as pd
import requests
import h5py
import json
import os
import numpy as np

for chr_i in [8, 17]:
	print(chr_i)
	data1 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_0_200.h5", key='data')
	data2 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_200_400.h5", key='data')
	data3 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_400_600.h5", key='data')
	data4 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_600_800.h5", key='data')
	data5 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_800_1000.h5", key='data')
	data6 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1000_1200.h5", key='data')
	data7 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1200_1400.h5", key='data')
	data8 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1400_1500.h5", key='data')
	data9 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1500_1632.h5", key='data')
	tdf = pd.concat([data1, data2], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data3], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data4], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data5], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data6], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data7], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data8], axis = 0).reset_index(drop=True)
	tdf = pd.concat([tdf, data9], axis = 0).reset_index(drop=True)
	tdf.shape
	filename = "./hg38_chr" + str(chr_i) + "_200data.h5"
	tdf.to_hdf(filename, key='data', mode='w')

df_read = pd.read_hdf(filename, key='data')
df_read



#HERE to get corr GPU, note because of potables incompatibility there is weird way to read the file
#salloc --partition=gpu_debug --reservation=gpu_debug
salloc -p gpu -C gpu80 # if out of memory error; experimented with mask clusters, errors
#salloc -A w23_biothreat_g -p gpu

salloc --partition=gpu_debug --reservation=gpu_debug
module load python/3.10-anaconda-2023.03
module load cuda/11.5
conda activate rapids-23.06
python

import cudf
import numpy as np
import pandas as pd
import h5py
from scipy.spatial.distance import pdist, squareform
import time
from datetime import timedelta

#path = "/lustre/scratch5/akim/results19/fin"
#filename = path + 'hg19_200data'
#df_corfile = pd.read_csv(path + "1000bp_chr6_hg38_cor_dist.csv", sep = '\t')
#df_corfile

#17, 8, 1


for chromosome_n in [8, 17]:
    path = "/lustre/scratch5/akim/results19/fin/"
    filename = path + "hg19_chr" + str(chromosome_n) + "_200data"
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
    df_read.index = range(0,1479)
    df_read = df_read.fillna(0) #hg19, chr 1, 8, 17 
    df_transposed = df_read.transpose()
    df_cudf = cudf.DataFrame.from_pandas(df_transposed)
    correlation_matrix = df_cudf.corr().transpose()
    correlation_matrix = 1 - correlation_matrix
    #print(correlation_matrix)
    correlation_matrix.to_csv(filename + 'correlation.h5')
    pd.read_csv(filename + 'correlation.h5', index_col=0)
    print(chromosome_n, " ", df_read.shape)
    end_time = time.monotonic()
    print(timedelta(seconds=end_time - start_time))


#hg19

for chromosome_n in [22]:
    path = "/lustre/scratch5/akim/results19/fin/"
    filename = path + "hg19_chr" + str(chromosome_n) + "_200data"
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
    df_read.index = range(0,1479)
    df_read = df_read.fillna(0) #hg19, chr 1, 8, 17 
    df_read.index = pd.to_numeric(df_read.index)
    dfch = pd.read_csv("df_hg19_chrome.csv")
    dfav = pd.read_csv("df_hg19_avocado.csv")
    indices_to_select = dfch.iloc[:, 0].tolist()
    ch_df = df_read.loc[indices_to_select]
    indices_to_select = dfav.iloc[:, 0].tolist()
    av_df = df_read.loc[indices_to_select]
    df_transposed = ch_df.transpose()
    df_cudf = cudf.DataFrame.from_pandas(df_transposed)
    correlation_matrix = df_cudf.corr().transpose()
    correlation_matrix = 1 - correlation_matrix
    #print(correlation_matrix)
    correlation_matrix.to_csv(filename + 'correlation_chrome.h5')
    pd.read_csv(filename + 'correlation_chrome.h5', index_col=0)
    print(chromosome_n, " ", ch_df.shape)
    df_transposed = av_df.transpose()
    df_cudf = cudf.DataFrame.from_pandas(df_transposed)
    correlation_matrix = df_cudf.corr().transpose()
    correlation_matrix = 1 - correlation_matrix
    #print(correlation_matrix)
    correlation_matrix.to_csv(filename + 'correlation_avocado.h5')
    pd.read_csv(filename + 'correlation_avocado.h5', index_col=0)
    print(chromosome_n, " ", av_df.shape)
    end_time = time.monotonic()
    print(timedelta(seconds=end_time - start_time))(base) akim@ch-fe1:~> 
