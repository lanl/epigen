# produce_go_results.py
#
# Code to generate Figure 10 and 11.
#
# For each chromosome this will generate a list of important genes.
# Fig 11. Can be produced by plugging these genes into KEGG.
# https://www.genome.jp/kegg/
#
# Input required:
#     ./data/genome_df38.csv
#     Binned data for all chromosomes:
#     ./chr_files/hg38_chr{}_200datacorrelation.h5
#
# Output:
#     allchrgenes_tokobas.txt
#     subsetchrgenes_tokobas.txt
# For more information, please see figures 10 and 11 in the paper.
#
# Time to run:
# Several hours.
#
#

import pandas as pd
from collections import Counter
from scipy.cluster.hierarchy import fcluster
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform
import numpy as np
import plotly.express as px

def check_symmetric(arr):
    if arr.shape[0] != arr.shape[1]:
        raise ValueError("The given array is not square!")
    non_symmetric_indices = []
    n = arr.shape[0]
    for i in range(n):
        for j in range(i + 1, n):  # Only check the upper triangle
            if arr[i, j] != arr[j, i]:
                non_symmetric_indices.append((i, j))
    return non_symmetric_indices


def make_symmetric(mat):
    rows, cols = mat.shape
    for i in range(rows):
        for j in range(i + 1, cols):  # only consider upper triangular part
            if mat[i, j] != mat[j, i]:  # unsymmetrical
                # Take the minimum of the two unsymmetrical entries
                symmetric_value = min(mat[i, j], mat[j, i])
                mat[i, j] = symmetric_value
                mat[j, i] = symmetric_value
    return mat


def get_peaks_matching_genes(binwidth, d, ret_df, setbp):
    start_p = list(map(int, d.index))
    start_df = pd.DataFrame()
    new_df = pd.DataFrame()
    Nlimit = len(start_p)
    for k in range(Nlimit):
        if new_df.shape[0] < Nlimit:
            gdf = ret_df[(ret_df["start"] >= start_p[k] * binwidth - setbp) & (ret_df["end"] <= start_p[k] * binwidth + binwidth + setbp)]
            if gdf.shape[0] != 0:
                for t in np.unique(gdf["GeneName"]):
                    temp = gdf[gdf["GeneName"] == t]
                    data = {"GeneName": [t], "start": np.min(temp["start"]), "end": np.max(temp["end"])}
                    new_df = pd.concat([start_df, pd.DataFrame(data)], axis=0).reset_index(drop=True)
                    start_df = new_df

            if new_df.shape[0] > 0:
                new_df = new_df.sort_values("start")
                new_df.reset_index(inplace=True)
                del new_df["index"]

    return new_df

#metadata
marks = pd.read_excel("./data/target_activity_factor.xlsx")
marks.columns = ["Target", "Activity", "Factor"]
df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ["Accession", "Target", "Biosample term name", "Genome"]]
merged_df = pd.merge(df38, marks, on="Target")
merged_df

# linkages for all  chrs
list_linkages = []
for i in range(1, 24):
    chr_id_1 = i
    if chr_id_1 == 23:
        chr_id_1 = "X"
    df_corr_1 = pd.read_csv("../chr_files/hg38_chr" + str(chr_id_1) + "_200data" + "correlation.h5", index_col=0)
    cor_dist_1 = df_corr_1.to_numpy()
    np.fill_diagonal(cor_dist_1, 0)
    indices = check_symmetric(cor_dist_1)
    if len(indices) != 0:
        cor_dist_1 = make_symmetric(cor_dist_1)
    condensed_dist_1 = squareform(cor_dist_1)
    linkresult_1 = sch.linkage(condensed_dist_1, method="complete")
    linkresult_1[linkresult_1 < 0] = 0
    list_linkages.append(linkresult_1)
len(list_linkages)


# get a summary df with number of important genes for all chrs 
thrs = 0.05 #p-value threshold to define important genes
start_p = 0
divisor = 200 #binwidth
cut_dist = 0.3 #distance at which the dendrogram was cut
buffer_bp = 500 #buffer 


# get summary dataframe, takes ~10 mins per chromosome
ll = []
lu = []
lc = []
chr_length = pd.read_csv("./data/hg38length.txt", sep="\t", header=None)
chr_length.columns = ["Chromosome", "Total length (bp)", "Nm", "Nmm"]
for chr_id in range(1, 24):
    print(chr_id)
    if chr_id == 23:
        end_p = (int(chr_length.loc[chr_length["Chromosome"] == "X", "Total length (bp)"]) // divisor) * divisor
    else:
        end_p = (int(chr_length.loc[chr_length["Chromosome"] == str(chr_id), "Total length (bp)"]) // divisor) * divisor
    num_steps = (end_p - start_p) // divisor
    mapping = {i: i * divisor for i in range(int(end_p / divisor + 1))}
    pv_df = pd.read_pickle("../chr_files/pvdf_" + str(chr_id) + ".pkl")
    ret_df = pd.read_csv("../chr_files/chr" + str(chr_id) + "_ret_df.csv", sep="\t")
    linkresult = list_linkages[chr_id - 1]
    clusters = fcluster(linkresult, cut_dist, criterion="distance")
    count = 0
    start_cldf = pd.DataFrame(
        {
            "ClusterSize": [],
            "EpiMark": [],
            "BioCell": [],
            "Activity": [],
            "Factor": [],
            "PeakKnownGenes": [],
            "HowManyKnownGenes": [],
            "PeakUnknownRegions(id)": [],
            "PeakUnknownRegions(start_position_200bp_long)": [],
            "HowManyUnknownRegions(200bp_long)": [],
        }
    )
    single_start_cldf = pd.DataFrame()
    genes_list = [] # collect known important genes from all clusters
    for i in range(len(np.unique(clusters))):
        outgenes_list = []
        idx_one_cluster = np.where(clusters == i)[0]
        subset_cluster = pv_df.loc[idx_one_cluster]
        if subset_cluster.shape[0] >= 2:  # consider clusters of at least size 2
            d = pd.DataFrame(np.max(subset_cluster, axis=0))
            d.columns = ["max"]
            d = d.loc[d["max"] <= thrs, "max"]
            if len(d) > 0:
                count = count + 1
                ret_dfi = get_peaks_matching_genes(divisor, d, ret_df, buffer_bp)
                if ret_dfi.shape[0] != 0:
                    outgenes_list.append(ret_dfi["GeneName"])
                    genes_list.append(list(outgenes_list[0]))
                    dd = {
                        "ClusterSize": [len(idx_one_cluster)],
                        "EpiMark": [list(df38.iloc[idx_one_cluster, 1])],
                        "BioCell": [list(df38.iloc[idx_one_cluster, 2])],
                        "Activity": [list(merged_df.iloc[idx_one_cluster, 4])],
                        "Factor": [list(merged_df.iloc[idx_one_cluster, 5])],
                        "PeakKnownGenes": [list(outgenes_list[0])],
                        "HowManyKnownGenes": [len(outgenes_list[0])],
                        "PeakUnknownRegions(id)": [list(d.index)],
                        "PeakUnknownRegions(start_position_200bp_long)": [list(map(lambda x: mapping[int(x)], list(d.index)))],
                        "HowManyUnknownRegions(200bp_long)": [len(d.index)],
                    }
                else:
                    dd = {
                        "ClusterSize": [len(idx_one_cluster)],
                        "EpiMark": [list(df38.iloc[idx_one_cluster, 1])],
                        "BioCell": [list(df38.iloc[idx_one_cluster, 2])],
                        "Activity": [list(merged_df.iloc[idx_one_cluster, 4])],
                        "Factor": [list(merged_df.iloc[idx_one_cluster, 5])],
                        "PeakKnownGenes": [list(outgenes_list)],
                        "HowManyKnownGenes": [0],
                        "PeakUnknownRegions(id)": [0],
                        "PeakUnknownRegions(start_position_200bp_long)": [0],
                        "HowManyUnknownRegions(200bp_long)": [0],
                    }
                ret_cldf = pd.concat([start_cldf, pd.DataFrame(dd)], axis=0).reset_index(drop=True)
                start_cldf = ret_cldf
    unique_genes_in_clusters = np.unique([item for sublist in genes_list for item in sublist])
    lu.append(unique_genes_in_clusters) # list of arrays of known important genes for each chr 
    ll.append(len(unique_genes_in_clusters)) #how many known genes are in each chr
    lc.append(np.sum(ret_cldf.loc[ret_cldf["HowManyKnownGenes"] > 0, "ClusterSize"]))

#Example output of ret_cldf for chr_1, there are 14 known important genes and their names are in the "PeakKnownGenes" column
# ret_cldf.iloc[205]
# Out[39]: 
# ClusterSize                                                                                    5.0
# EpiMark                                              [H3K9me3, H3K9me3, H3K9me3, H3K9me3, H3K9me3]
# BioCell                                          [common myeloid progenitor, CD34-positive, com...
# Activity                                         [Repressive, Repressive, Repressive, Repressiv...
# Factor                                           [Histone modification, Histone modification, H...
# PeakKnownGenes                                   [MIR4255, MIR4421, MIR4421, MIR6500, MIR6500, ...
# HowManyKnownGenes                                                                             14.0
# PeakUnknownRegions(id)                           [3944, 4959, 4960, 4963, 4964, 4966, 5499, 550...
# PeakUnknownRegions(start_position_200bp_long)    [788800, 991800, 992000, 992600, 992800, 99320...
# HowManyUnknownRegions(200bp_long)                                                           8415.0
# Name: 205, dtype: object


# selected_columns = ['EpiMark', 'BioCell', 'Activity', 'Factor', 'PeakKnownGenes']
# for index, row in ret_cldf.iterrows():
#     print(index, row['ClusterSize'])
#     Counter(row['EpiMark'])
#     Counter(row['BioCell'])
#     Counter(row['Activity'])
#     Counter(row['Factor'])
#     flattened_data = row['PeakKnownGenes']
#     flattened_data = [re.match('([a-zA-Z]+)', i).group() for i in flattened_data]
#     counts = Counter(flattened_data)
#     print("Counts:", counts)
#     print("--------------------------------------------------------------------------------------------")


#input is "lu" from above code which is list of arrays, each array contains names of important genes for each chr, list size is 23 corresponding to the # of chrs

all_genes = []
combined = []

for idx, gene_array in enumerate(lu, start=1):
    for gene in gene_array:
        all_genes.append(gene)
        combined.append((gene, idx))

with open("allchr_impgenes_ontology.txt", "w") as f: #all genes
    for gene in all_genes:
        f.write(f"{gene}\n")



# important genes for each chromosome, uncomment all below, we already printed out important genes for each chr
# chr1 = [
#     "COA6-AS1",
#     "DPM3",
#     "ENO1-AS1",
#     "G0S2",
#     "GAS5-AS1",
#     "H2AC20",
#     "H2AC21",
#     "H2AW" "H2BC20P",
#     "H2BU1",
#     "H3C13",
#     "HES4",
#     "LCE3A",
#     "LOC100287049",
#     "LOC105376805" "LOC105378591",
#     "LOC105378933",
#     "MIR11399",
#     "MIR1182",
#     "MIR12133",
#     "MIR137" "MIR1537",
#     "MIR1976",
#     "MIR2682",
#     "MIR30C1",
#     "MIR30E",
#     "MIR3124",
#     "MIR3605" "MIR3620",
#     "MIR3671",
#     "MIR3917",
#     "MIR3972",
#     "MIR4257",
#     "MIR4258",
#     "MIR4259" "MIR4420",
#     "MIR4425",
#     "MIR4632",
#     "MIR4781",
#     "MIR5087",
#     "MIR5187",
#     "MIR548AC" "MIR555",
#     "MIR5581",
#     "MIR6068",
#     "MIR6077",
#     "MIR6084",
#     "MIR6727",
#     "MIR6728" "MIR6731",
#     "MIR6732",
#     "MIR6733",
#     "MIR6734",
#     "MIR6735",
#     "MIR6737",
#     "MIR6738" "MIR6741",
#     "MIR6769B",
#     "MIR760",
#     "MIR7846",
#     "MIR7852",
#     "MIR9-1",
#     "MIR92B" "PACERR",
#     "PCAT6",
#     "PYDC5",
#     "RFX5-AS1",
#     "RNF207-AS1",
#     "RNU11",
#     "RNU5D-1" "RNU5E-1",
#     "RNU5F-1",
#     "RNVU1-1",
#     "RNVU1-14",
#     "RNVU1-15",
#     "RNVU1-17",
#     "RNVU1-19" "RNVU1-20",
#     "RNVU1-2A",
#     "RNVU1-3",
#     "RNVU1-4",
#     "RNVU1-6",
#     "RNVU1-7",
#     "RNVU1-8" "SCARNA18B",
#     "SCARNA2",
#     "SCARNA4",
#     "SNORA100",
#     "SNORA103",
#     "SNORA14B" "SNORA16A",
#     "SNORA44",
#     "SNORA55",
#     "SNORA58B",
#     "SNORA61",
#     "SNORA66",
#     "SNORA73A" "SNORA73B",
#     "SNORD145",
#     "SNORD167",
#     "SNORD21",
#     "SNORD38A",
#     "SNORD38B" "SNORD44",
#     "SNORD45A",
#     "SNORD45B",
#     "SNORD45C",
#     "SNORD46",
#     "SNORD55",
#     "SNORD74" "SNORD75",
#     "SNORD76",
#     "SNORD77",
#     "SNORD78",
#     "SNORD79",
#     "SNORD80",
#     "SNORD99",
#     "ZNF593",
# ]
# chr2 = [
#     "AGBL5-AS1",
#     "BCYRN1",
#     "HAGLROS",
#     "HDAC4-AS1",
#     "IDH1-AS1",
#     "LOC105373562",
#     "LOC105373656",
#     "LOC105373759",
#     "MIR10B",
#     "MIR1244-1",
#     "MIR1244-2",
#     "MIR1244-3",
#     "MIR1244-4",
#     "MIR153-1",
#     "MIR2355",
#     "MIR26B",
#     "MIR3128",
#     "MIR3131",
#     "MIR3679",
#     "MIR375",
#     "MIR4434",
#     "MIR4437",
#     "MIR4444-1",
#     "MIR4444-2",
#     "MIR4757",
#     "MIR4777",
#     "MIR4779",
#     "MIR4783",
#     "MIR4784",
#     "MIR4785",
#     "MIR5001",
#     "MIR5192",
#     "MIR5696",
#     "MIR5703",
#     "MIR6513",
#     "MIR7704",
#     "MIR7845",
#     "MIR933",
#     "MIR9899",
#     "MRPL53",
#     "MTLN",
#     "NAT8B",
#     "RNU4ATAC",
#     "SCARNA5",
#     "SNAR-H",
#     "SNORA75",
#     "SNORA80B",
#     "SNORD20",
#     "SNORD53B",
#     "SNORD70",
#     "SNORD82",
#     "SNORD89",
#     "TMSB10",
# ]
# chr3 = [
#     "LINC01063",
#     "LOC100289361",
#     "MIR12124",
#     "MIR12127",
#     "MIR1224",
#     "MIR1226",
#     "MIR1248",
#     "MIR1284",
#     "MIR15B",
#     "MIR16-2",
#     "MIR191",
#     "MIR28",
#     "MIR3136",
#     "MIR3714",
#     "MIR425",
#     "MIR4273",
#     "MIR4442",
#     "MIR4443",
#     "MIR4787",
#     "MIR4795",
#     "MIR5193",
#     "MIR544B",
#     "MIR5588",
#     "MIR564",
#     "MIR568",
#     "MIR5787",
#     "MIR6529",
#     "MIR6823",
#     "MIR6824",
#     "MIR6825",
#     "MIR6826",
#     "MIR6872",
#     "MIR885",
#     "MIR922",
#     "MIRLET7G",
#     "NCBP2AS2",
#     "RASSF1-AS1",
#     "SCAANT1",
#     "SNORA4",
#     "SNORA6",
#     "SNORA62",
#     "SNORA63",
#     "SNORA63B",
#     "SNORA7A",
#     "SNORA7B",
#     "SNORA81",
#     "SNORD13J",
#     "SNORD2",
#     "SNORD66",
#     "TERC",
# ]
# chr4 = [
#     "FLJ20021",
#     "LOC105377590",
#     "MIR12115",
#     "MIR302A",
#     "MIR302C",
#     "MIR302D",
#     "MIR3138",
#     "MIR367",
#     "MIR3684",
#     "MIR3945",
#     "MIR4276",
#     "MIR4449",
#     "MIR4453",
#     "MIR4801",
#     "MIR5091",
#     "MIR572",
#     "MIR574",
#     "SLED1",
#     "SNHG8",
#     "SNORA24",
#     "SNORA26",
#     "SNORD144",
#     "SNORD73B",
#     "STIM2-AS1",
# ]
# chr5 = [
#     "CXXC5-AS1",
#     "HIGD2A",
#     "LINC01023",
#     "LINC01574",
#     "LINC01962",
#     "LINC02106",
#     "LOC100303749",
#     "LOC731157",
#     "MIR143",
#     "MIR145",
#     "MIR146A",
#     "MIR2277",
#     "MIR3142",
#     "MIR3650",
#     "MIR3655",
#     "MIR3661",
#     "MIR378A",
#     "MIR378H",
#     "MIR3912",
#     "MIR4281",
#     "MIR449A",
#     "MIR449B",
#     "MIR449C",
#     "MIR4634",
#     "MIR4638",
#     "MIR584",
#     "MIR6499",
#     "MIR6831",
#     "MIR8089",
#     "MIR9-2",
#     "PFN3",
#     "ROPN1L-AS1",
#     "SNORA13",
#     "SNORA74D",
#     "SNORD170",
#     "SNORD63B",
#     "SNORD72",
#     "SNORD95",
#     "SNORD96A",
#     "VTRNA1-1",
#     "VTRNA1-3",
#     "VTRNA2-1",
# ]
# chr6 = [
#     "C6orf226",
#     "CAHM",
#     "DDX39B-AS1",
#     "DINOL",
#     "H1-1",
#     "H1-2",
#     "H1-3",
#     "H1-4",
#     "H1-5",
#     "H2AC11",
#     "H2AC12",
#     "H2AC13",
#     "H2AC14",
#     "H2AC15",
#     "H2AC16",
#     "H2AC17",
#     "H2AC4",
#     "H2AC6",
#     "H2AC7",
#     "H2AC8",
#     "H2BC10",
#     "H2BC11",
#     "H2BC13",
#     "H2BC14",
#     "H2BC15",
#     "H2BC17",
#     "H2BC3",
#     "H2BC6",
#     "H2BC7",
#     "H2BC8",
#     "H2BC9",
#     "H3C1",
#     "H3C10",
#     "H3C11",
#     "H3C12",
#     "H3C2",
#     "H3C3",
#     "H3C7",
#     "H3C8",
#     "H4C1",
#     "H4C11",
#     "H4C12",
#     "H4C13",
#     "H4C2",
#     "H4C3",
#     "H4C4",
#     "H4C5",
#     "H4C6",
#     "H4C8",
#     "H4C9",
#     "HCG14",
#     "IER3-AS1",
#     "JARID2-AS1",
#     "LINC01556",
#     "LOC100270746",
#     "MIR10398",
#     "MIR1236",
#     "MIR219A1",
#     "MIR3143",
#     "MIR3145",
#     "MIR3918",
#     "MIR3925",
#     "MIR3939",
#     "MIR4282",
#     "MIR4466",
#     "MIR4639",
#     "MIR4645",
#     "MIR4646",
#     "MIR4647",
#     "MIR5690",
#     "MIR6720",
#     "MIR6832",
#     "MIR6833",
#     "MIR6834",
#     "MIR6891",
#     "MIR7111",
#     "MIR7161",
#     "MIR877",
#     "RN7SK",
#     "SF3B5",
#     "SNORA38",
#     "SNORD100",
#     "SNORD101",
#     "SNORD117",
#     "SNORD141A",
#     "SNORD141B",
#     "SNORD166",
#     "SNORD48",
#     "SNORD50A",
#     "SNORD50B",
#     "SNORD52",
#     "SNORD84",
# ]
# chr7 = [
#     "FERD3L",
#     "LOC100128325",
#     "LOC101927420",
#     "MIR10399",
#     "MIR106B",
#     "MIR12119",
#     "MIR148A",
#     "MIR196B",
#     "MIR25",
#     "MIR3907",
#     "MIR4648",
#     "MIR4651",
#     "MIR4657",
#     "MIR4658",
#     "MIR5090",
#     "MIR550A1",
#     "MIR550B1",
#     "MIR590",
#     "MIR6509",
#     "MIR6837",
#     "MIR6838",
#     "MIR6840",
#     "MIR6874",
#     "MIR6875",
#     "MIR6892",
#     "MIR93",
#     "MNX1-AS2",
#     "RNY1",
#     "RNY3",
#     "RNY4",
#     "RNY5",
#     "RPS2P32",
#     "SNORA5A",
#     "SNORA5B",
#     "SNORA5C",
#     "SNORA9",
#     "TAS2R40",
#     "UFSP1",
# ]
# chr8 = [
#     "CERNA3",
#     "MIR10400",
#     "MIR1204",
#     "MIR1205",
#     "MIR1207",
#     "MIR124-1",
#     "MIR124-2",
#     "MIR1322",
#     "MIR3148",
#     "MIR3150A",
#     "MIR3150B",
#     "MIR320A",
#     "MIR3610",
#     "MIR378D2",
#     "MIR4469",
#     "MIR4470",
#     "MIR4471",
#     "MIR4664",
#     "MIR5194",
#     "MIR5708",
#     "MIR596",
#     "MIR661",
#     "MIR6842",
#     "MIR6843",
#     "MIR6847",
#     "MIR6850",
#     "MIR6876",
#     "MIR7848",
#     "MOS",
#     "RBM12B-AS1",
#     "SMPD5",
#     "SNORA1B",
#     "SNORA99",
#     "SNORD13",
#     "SNORD54",
#     "SNORD87",
#     "YTHDF3-AS1",
# ]
# chr9 = [
#     "BNC2-AS1",
#     "CCL27",
#     "CDKN2A-DT",
#     "IFT74-AS1",
#     "LINC01230",
#     "LOC105376271",
#     "MIR101-2",
#     "MIR12126",
#     "MIR126",
#     "MIR199B",
#     "MIR2861",
#     "MIR3153",
#     "MIR3154",
#     "MIR3621",
#     "MIR3651",
#     "MIR3960",
#     "MIR4290",
#     "MIR4479",
#     "MIR4665",
#     "MIR4669",
#     "MIR4672",
#     "MIR4674",
#     "MIR548AW",
#     "MIR600",
#     "MIR601",
#     "MIR6722",
#     "MIR6851",
#     "MIR6852",
#     "MIR6853",
#     "MIR6854",
#     "MIR6855",
#     "MIR7114",
#     "MIR873",
#     "MIRLET7A1",
#     "MIRLET7F1",
#     "MRPL41",
#     "RMRP",
#     "RNU6ATAC",
#     "SNORA17A",
#     "SNORA17B",
#     "SNORA65",
#     "SNORA84",
#     "SNORD121A",
#     "SNORD24",
#     "SNORD36A",
#     "SNORD36B",
#     "SNORD36C",
# ]
# chr10 = [
#     "LINC01167",
#     "LOC107984208",
#     "LOC107984236",
#     "MIR1307",
#     "MIR146B",
#     "MIR1915",
#     "MIR202",
#     "MIR2110",
#     "MIR3155A",
#     "MIR3155B",
#     "MIR3663",
#     "MIR3941",
#     "MIR4482",
#     "MIR4675",
#     "MIR4678",
#     "MIR4683",
#     "MIR4685",
#     "MIR548AK",
#     "MIR8086",
#     "NEBL-AS1",
#     "PARD3-AS1",
#     "SNORA12",
#     "SNORD129",
#     "SNORD172",
#     "SNORD3J",
# ]
# chr11 = [
#     "CCDC85B",
#     "FBXO3-DT",
#     "LINC02716",
#     "LOC101928837",
#     "LOC644656",
#     "MASCRNA",
#     "MIR10392",
#     "MIR125B1",
#     "MIR1260B",
#     "MIR129-2",
#     "MIR130A",
#     "MIR139",
#     "MIR1908",
#     "MIR210",
#     "MIR3161",
#     "MIR3165",
#     "MIR326",
#     "MIR34B",
#     "MIR34C",
#     "MIR3654",
#     "MIR4298",
#     "MIR4487",
#     "MIR4488",
#     "MIR4489",
#     "MIR4492",
#     "MIR4687",
#     "MIR4690",
#     "MIR548AL",
#     "MIR5691",
#     "MIR6073",
#     "MIR6090",
#     "MIR611",
#     "MIR612",
#     "MIR6514",
#     "MIR6743",
#     "MIR6745",
#     "MIR6747",
#     "MIR6748",
#     "MIR6751",
#     "MIR6753",
#     "MIR7155",
#     "MIR7847",
#     "OR52W1",
#     "OR8G2P",
#     "SNORA3A",
#     "SNORA3B",
#     "SNORA40",
#     "SNORA52",
#     "SNORA54",
#     "SNORA57",
#     "SNORD13F",
#     "SNORD147",
#     "SNORD14B",
#     "SNORD14C",
#     "SNORD14D",
#     "SNORD14E",
#     "SNORD150",
#     "SNORD15A",
#     "SNORD15B",
#     "SNORD164",
#     "SNORD22",
#     "SNORD25",
#     "SNORD26",
#     "SNORD27",
#     "SNORD28",
#     "SNORD29",
#     "SNORD30",
#     "SNORD31",
#     "SNORD97",
#     "TOLLIP-AS1",
# ]
# chr12 = [
#     "ATXN2-AS",
#     "H4-16",
#     "HOXC-AS1",
#     "KRT19P2",
#     "LINC01405",
#     "LOC101928530",
#     "LOC105369728",
#     "MIR1291",
#     "MIR148B",
#     "MIR196A2",
#     "MIR200C",
#     "MIR26A2",
#     "MIR3649",
#     "MIR3652",
#     "MIR3913-1",
#     "MIR3913-2",
#     "MIR4496",
#     "MIR4498",
#     "MIR492",
#     "MIR5188",
#     "MIR548C",
#     "MIR548Z",
#     "MIR6125",
#     "MIR614",
#     "MIR615",
#     "MIR616",
#     "MIR618",
#     "MIR619",
#     "MIR620",
#     "MIR6758",
#     "MIR6760",
#     "MIR8072",
#     "MIRLET7I",
#     "PCBP2-OT1",
#     "RNU4-1",
#     "RNU4-2",
#     "RNU7-1",
#     "SBNO1-AS1",
#     "SCARNA11",
#     "SNORA120",
#     "SNORA2B",
#     "SNORA2C",
#     "SNORA49",
#     "SNORA70G",
#     "SNORD59A",
#     "SNORD59B",
# ]
# chr13 = ["FKSG29", "LOC112268114", "MIR3613", "MIR3665", "MIR623", "MIR8073", "SNORA31", "SNORA31B"]
# chr14 = [
#     "ECRP",
#     "LINC02332",
#     "LOC105370705",
#     "LOC730202",
#     "MIR12121",
#     "MIR1247",
#     "MIR151B",
#     "MIR203A",
#     "MIR203B",
#     "MIR3173",
#     "MIR342",
#     "MIR4505",
#     "MIR4507",
#     "MIR4538",
#     "MIR4539",
#     "MIR4706",
#     "MIR4707",
#     "MIR4710",
#     "MIR5580",
#     "MIR624",
#     "MIR7843",
#     "MIR9718",
#     "RN7SL1",
#     "RN7SL2",
#     "RN7SL3",
#     "RNASE3",
#     "RPPH1",
#     "SCARNA13",
#     "SEC23A-AS1",
#     "SNORA11B",
#     "SNORA28",
#     "SNORD56B",
#     "SNORD8",
#     "SNORD9",
# ]
# chr15 = [
#     "LOC100131315",
#     "MIR10393",
#     "MIR1179",
#     "MIR1282",
#     "MIR1469",
#     "MIR184",
#     "MIR3174",
#     "MIR3175",
#     "MIR422A",
#     "MIR4312",
#     "MIR4512",
#     "MIR4513",
#     "MIR4515",
#     "MIR4715",
#     "MIR626",
#     "MIR629",
#     "MIR6766",
#     "MIR7706",
#     "MIR9-3",
#     "RNU5A-1",
#     "RNU5B-1",
#     "RNU6-1",
#     "RNU6-2",
#     "RNU6-7",
#     "RNU6-8",
#     "RNU6-9",
#     "SNORD116-20",
#     "SNORD16",
#     "SNORD18A",
# ]
# chr16 = [
#     "AGRP",
#     "ATP2A1-AS1",
#     "C16orf91",
#     "CAPNS2",
#     "HBM",
#     "HBQ1",
#     "LINC00235",
#     "LINC02168",
#     "LOC101928659",
#     "LOC105371038",
#     "LOC107984813",
#     "LOC554206",
#     "MIR11401",
#     "MIR138-2",
#     "MIR140",
#     "MIR1538",
#     "MIR193B",
#     "MIR3178",
#     "MIR3180-4",
#     "MIR3180-5",
#     "MIR3181",
#     "MIR365A",
#     "MIR4518",
#     "MIR4519",
#     "MIR484",
#     "MIR5093",
#     "MIR5189",
#     "MIR548H2",
#     "MIR5587",
#     "MIR6506",
#     "MIR6771",
#     "MIR6774",
#     "MIR762",
#     "MIR7854",
#     "MIR940",
#     "MT2A",
#     "NPW",
#     "OR1F2P",
#     "SNHG19",
#     "SNHG9",
#     "SNORA10",
#     "SNORA119",
#     "SNORA30",
#     "SNORA3C",
#     "SNORA64",
#     "SNORA78",
#     "SNORD111B",
#     "SNORD60",
#     "SNORD68",
# ]
# chr17 = [
#     "CCDC182",
#     "EMC6",
#     "FAM215A",
#     "KRTAP29-1",
#     "LINC02079",
#     "LOC105371592",
#     "MIR10A",
#     "MIR1253",
#     "MIR132",
#     "MIR142",
#     "MIR152",
#     "MIR193A",
#     "MIR196A1",
#     "MIR21",
#     "MIR212",
#     "MIR22",
#     "MIR301A",
#     "MIR3064",
#     "MIR3184",
#     "MIR3185",
#     "MIR324",
#     "MIR3614",
#     "MIR3615",
#     "MIR3678",
#     "MIR423",
#     "MIR4316",
#     "MIR4520-1",
#     "MIR4520-2",
#     "MIR4521",
#     "MIR4522",
#     "MIR4523",
#     "MIR4525",
#     "MIR4723",
#     "MIR4725",
#     "MIR4727",
#     "MIR4729",
#     "MIR4733",
#     "MIR4734",
#     "MIR4736",
#     "MIR4738",
#     "MIR4740",
#     "MIR5047",
#     "MIR6080",
#     "MIR632",
#     "MIR636",
#     "MIR6516",
#     "MIR6779",
#     "MIR6780A",
#     "MIR6781",
#     "MIR6787",
#     "MIR6864",
#     "MIR6865",
#     "MIR6866",
#     "MIR6884",
#     "MTVR2",
#     "NPB",
#     "OR4D2",
#     "PITPNA-AS1",
#     "PRAC1",
#     "RPRML",
#     "SCARNA16",
#     "SCARNA20",
#     "SCARNA21",
#     "SNHG25",
#     "SNORA21",
#     "SNORA21B",
#     "SNORA48",
#     "SNORA50C",
#     "SNORA67",
#     "SNORD10",
#     "SNORD104",
#     "SNORD118",
#     "SNORD124",
#     "SNORD1A",
#     "SNORD1B",
#     "SNORD1C",
#     "SNORD3A",
#     "SNORD3D",
#     "SNORD42A",
#     "SNORD42B",
#     "SNORD49A",
#     "SNORD49B",
#     "SNORD4A",
#     "SNORD4B",
#     "SNORD65",
#     "TMEM256",
#     "TMEM88",
# ]
# chr18 = ["MIR1539", "MIR320C1", "MIR4526", "MIR4741", "MIR7153", "MIR8078", "SCARNA17", "SNORA37", "SNORD58A", "SNORD58B"]
# chr19 = [
#     "C19orf73",
#     "GNG8",
#     "LDLR-AS1",
#     "LINC02560",
#     "LOC102725254",
#     "MIR10394",
#     "MIR1181",
#     "MIR1199",
#     "MIR1227",
#     "MIR1238",
#     "MIR1470",
#     "MIR150",
#     "MIR23A",
#     "MIR24-2",
#     "MIR27A",
#     "MIR3188",
#     "MIR3189",
#     "MIR3190",
#     "MIR3191",
#     "MIR330",
#     "MIR3940",
#     "MIR4321",
#     "MIR4322",
#     "MIR4323",
#     "MIR4530",
#     "MIR4745",
#     "MIR4748",
#     "MIR4751",
#     "MIR4752",
#     "MIR4754",
#     "MIR4999",
#     "MIR5088",
#     "MIR5196",
#     "MIR5684",
#     "MIR637",
#     "MIR638",
#     "MIR639",
#     "MIR640",
#     "MIR641",
#     "MIR6515",
#     "MIR6719",
#     "MIR6789",
#     "MIR6790",
#     "MIR6791",
#     "MIR6792",
#     "MIR6797",
#     "MIR6800",
#     "MIR6801",
#     "MIR6802",
#     "MIR6803",
#     "MIR6804",
#     "MIR6805",
#     "MIR6807",
#     "MIR6886",
#     "MIR7850",
#     "MIR8085",
#     "NKG7",
#     "OR2Z1",
#     "POLR2I",
#     "PTH2",
#     "SNAR-G1",
#     "SNORA118",
#     "SNORA68",
#     "SNORD105",
#     "SNORD105B",
#     "SNORD135",
#     "SNORD152",
#     "SNORD157",
#     "SNORD175",
#     "SNORD23",
#     "SNORD32A",
#     "SNORD33",
#     "SNORD34",
#     "SNORD35A",
#     "SNORD35B",
#     "SNORD37",
#     "SNORD41",
#     "SNORD88C",
# ]
# chr20 = [
#     "FAM83C-AS1",
#     "LKAAEAR1",
#     "LOC105372493",
#     "MIR124-3",
#     "MIR1289-1",
#     "MIR1292",
#     "MIR1914",
#     "MIR3193",
#     "MIR3195",
#     "MIR4755",
#     "MIR647",
#     "MIR6813",
#     "MIR6869",
#     "OXT",
#     "PCNA-AS1",
#     "SNORA51",
#     "SNORA71B",
#     "SNORA71C",
#     "SNORA71D",
#     "SNORD110",
#     "SNORD119",
#     "SNORD12",
#     "SNORD12B",
#     "SNORD12C",
#     "SNORD17",
#     "SNORD56",
#     "SNORD86",
#     "TRERNA1",
# ]
# chr21 = [
#     "BRWD1-AS2",
#     "IL10RB-DT",
#     "KRTAP12-3",
#     "KRTAP12-4",
#     "KRTAP21-3",
#     "LOC101930094",
#     "MIR155",
#     "MIR3197",
#     "MIR6501",
#     "MIR6508",
#     "MIR6815" "URB1-AS1",
# ]
# chr22 = [
#     "LINC01311",
#     "LOC105373031",
#     "LOC91370",
#     "MIF",
#     "MIR1281",
#     "MIR1306",
#     "MIR185",
#     "MIR301B",
#     "MIR3198-1",
#     "MIR3199-1",
#     "MIR3199-2",
#     "MIR3618",
#     "MIR3928",
#     "MIR648",
#     "MIR650",
#     "MIR658",
#     "MIR6819",
#     "MIR6820",
#     "MIR6821",
#     "MIR6889",
#     "MIR7109",
#     "RNU12",
#     "SNORD139",
#     "SNORD43",
#     "SNORD83A",
#     "SNORD83B",
# ]
# chr23 = [
#     "EIF1AX-AS1",
#     "INTS6L-AS1",
#     "LINC02601",
#     "LOC105373383",
#     "MIR223",
#     "MIR3202-1",
#     "MIR3202-2",
#     "MIR374A",
#     "MIR424",
#     "MIR4536-1",
#     "MIR4536-2",
#     "MIR4767",
#     "MIR503",
#     "MIR545",
#     "MIR6089",
#     "MIR6894",
#     "MIR718",
#     "MIR766",
#     "SNORA69",
#     "SNORA70",
#     "SNORD61",
# ]  # chr X


# # Get files with the list of genes for input for Kobas and GO plots 9 and 10
# all_genes = []
# combined = []
# for index in range(1, 24):  # 1 to 23 inclusive
#     lst = locals()[f"chr{index}"]
#     for item in lst:
#         all_genes.append(item)
#         combined.append((item, index))

# with open("allchrgenes_tokobas.txt", "w") as f:
#     for item in all_genes:
#         f.write(f"{item}\n")
# element_counts = Counter(item for item, _ in combined)

# common_elements = {item for item, count in element_counts.items() if count > 1}

# element_to_lists = {}
# for element in common_elements:
#     element_to_lists[element] = [idx for item, idx in combined if item == element]

# for element, lists in element_to_lists.items():
#     print(f"Element {element} appears in lists: {', '.join(map(str, lists))}")

# filtered_list = [item for item in all_genes if not (item.startswith("SNOR") or item.startswith("MIR"))]

# with open("subsetchrgenes_tokobas.txt", "w") as f:
#     for item in filtered_list:
#         f.write(f"{item}\n")


#Figure 11 is obtained from Cytoscape
#Below is the code to get Figure 10 based on the output from gProfiler software for Gene Ontology
df = pd.read_csv("./data/gProfiler_hsapiens.csv") 
df = df.iloc[[0, 5, 7, 10, 11, 12, 13, 14, 15, 17, 36, 37, 43, 44, 47, 48, 49, 50, 56, 57], :] #only relevant output
df = df.iloc[:, [0, 1, 4, 8]]
df.columns = ['Source', 'Term Name', 'Adjusted p-value', 'Number of Genes']
replacement_map = {
    'GO:MF': 'molecular function',
    'GO:BP': 'biological process',
    'GO:CC': 'cellular component'
}

df['Source'] = df['Source'].replace(replacement_map)


data = pd.DataFrame({
    'Term': ['GO:0008150', 'GO:0003674', 'GO:0005575', 'GO:0008150', 'GO:0003674'],
    'Count': [100, 75, 50, 25, 10],
    'Type': ['Biological Process', 'Molecular Function', 'Cellular Component', 
             'Biological Process', 'Molecular Function']
})

data_sorted = data.sort_values(by='Count', ascending=True)

fig = px.bar(
    df,
    y='Term Name', 
    x='Number of Genes',
    color='Source',
    orientation='h',
    title='The Most Enriched GO Terms',
    #labels={'Count':'Number of genes', 'Term':'GO Term'},
    #color_discrete_map={'Biological Process':'green', 'Molecular Function':'red', 'Cellular Component':'blue'}
)

fig.update_layout(
    yaxis_title="GO Term",
    xaxis_title="Number of genes",
    legend_title="Source",
    yaxis=dict(autorange="reversed"), # This is to match the order in Figure 7
    plot_bgcolor='rgba(0,0,0,0)'
)

fig.show()
fig.write_image("./most_enriched.eps")
