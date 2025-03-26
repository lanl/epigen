# create_pvdf.py
#
# This reduces the signal bin counts created in the previous step.
# Since most bins do not contain any data across all samples,
# This script drops all bins which do not contain any data.
#
# The resulting files have ~50k bins instead of ~1M.
#
# Inputs to sript: range of chromosomes to analyze from i to j
# where i=sys.argv[0] and j=sys.argv[1]
#
#
# input file needed: hg38_<chrID>_200data files that are generated through
#        run_get_values_for_chromosome_script
# input file needed: genome_df38.csv contains Target variable that we used as
#        indicies for pvdf_chr_id.pkl dataframes. (Encode Accession #)
#
import pandas as pd
import numpy as np
import h5py
import sys


# Concatenate chunked files from get_values_from_chromosome.py
print("Concatenating chromosome data.")
for chr_i in range(1, 24):
    print("Chromosome:", chr_i)
    if chr_i == 23:
        chr_i = "X"
    data1 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_0_200.h5", key="data")
    data2 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_200_400.h5", key="data")
    data3 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_400_600.h5", key="data")
    data4 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_600_800.h5", key="data")
    data5 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_800_1000.h5", key="data")
    data6 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1000_1200.h5", key="data")
    data7 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1200_1400.h5", key="data")
    data8 = pd.read_hdf("./hg38_chr" + str(chr_i) + "_200_1500_1632.h5", key="data")
    tdf = pd.concat([data1, data2], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data3], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data4], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data5], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data6], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data7], axis=0).reset_index(drop=True)
    tdf = pd.concat([tdf, data8], axis=0).reset_index(drop=True)
    filename = "./hg38_chr" + str(chr_i) + "_200data.h5"  # saving combined h5 files for each chr
    tdf.to_hdf(filename, key="data", mode="w")


# Read genome metadata
df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ["Accession", "Target", "Biosample term name", "Genome"]]

# Loop over chromosomes
for chr_id in range(int(float(sys.argv[1])), int(float(sys.argv[2]))):
    if chr_id == 23:  # This is chromosome X
        chr_id = "X"

    print("Chromosome:", chr_id)
    # get file for chromosome
    filename = "./results38/hg38_chr" + str(chr_id) + "_200data"
    with h5py.File(filename + ".h5", "r") as hdf:
        data_group = hdf.get("data")
        # Retrieve each component of the DataFrame
        axis0 = np.array(data_group.get("axis0")).astype(str)
        axis1 = np.array(data_group.get("axis1")).astype(str)
        block0_items = np.array(data_group.get("block0_items")).astype(str)
        block0_values = np.array(data_group.get("block0_values"))

    # original data are signal p-values (-log10(regular p-value)).
    # here we'll only extract columns where regular p-value <= 0.05
    # this corresponds to a signal threshold of 1.3.
    # Then we convert this back to a p-value.
    signal_thresh = 1.3

    df = pd.DataFrame(data=block0_values, index=axis1, columns=axis0)
    df.index = df38["Target"]
    c = pd.DataFrame(df.max(axis="rows"))
    c.columns = ["max"]
    c = c.loc[c["max"] >= signal_thresh, "max"]
    pv_df = 10 ** (-df.iloc[:, list(map(int, c.index))])  # calculate p-value

    # assign sample numbers to each row. Order corresponds to genome_df38.csv file.
    pv_df.index = range(0, len(pv_df))
    pv_df.to_pickle("./results38/pvdf_" + str(chr_id) + ".pkl")
