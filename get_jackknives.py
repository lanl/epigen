#Obtain jackknives -- producing csv file for each replicate (600 for hg19 and 600 for hg38)
#input: need hg19_chr{chr_id}_200data.h5 and hg38_chr{chr_id}_200data.h5

import sys
import time
import pickle
import numpy as np
import pandas as pd
import h5py
from sklearn.metrics import pairwise_distances
from scipy.spatial.distance import squareform
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from collections import Counter


def filter_df_chr_by_accession(new_df, df_chr):
    return df_chr[df_chr.index.isin(new_df.index)]

def compute_cooccurrence(filtered_df, linkage_matrix, labels,
                         threshold=0.3, min_cluster_size=2):
    clusters = fcluster(linkage_matrix, threshold, criterion='distance')
    unique_labels = list(np.unique(labels))
    n_lbl = len(unique_labels)
    co_mat = np.zeros((n_lbl, n_lbl))
    global_counts = Counter(labels)
    for cid in np.unique(clusters):
        idxs = np.where(clusters == cid)[0]
        if len(idxs) < min_cluster_size:
            continue
        lbls = [labels[i] for i in idxs]
        for i in range(n_lbl):
            for j in range(i, n_lbl):
                c1 = lbls.count(unique_labels[i])
                c2 = lbls.count(unique_labels[j])
                if c1 and c2:
                    score = min(c1, c2) / np.sqrt(global_counts[unique_labels[i]] * global_counts[unique_labels[j]])
                else:
                    score = 0
                co_mat[i, j] += score
                co_mat[j, i] = co_mat[i, j]
    return co_mat, unique_labels

# ---------------------------
# Load and prepare metadata
# ---------------------------
marks = pd.read_excel("./data/target_activity_factor.xlsx")
marks.columns = ["Target", "Activity", "Factor"]

# load and merge genome_df38
df38 = pd.read_csv("./data/genome_df38.csv", usecols=['Accession','Target','Biosample term name','Genome'])
df38 = pd.merge(df38, marks, on='Target')

# load and merge genome_df19
df19 = pd.read_csv("./data/genome_df19.csv", usecols=['Accession','Target','Biosample term name','Genome'])
df19 = pd.merge(df19, marks, on='Target')

# restrict to shared modifiers
sample_idx = ['CTCF','EP300','H2AFZ','H2AK5ac','H2BK120ac','H2BK12ac','H2BK15ac',
              'H2BK5ac','H3K14ac','H3K18ac','H3K23ac','H3K27ac','H3K27me3',
              'H3K36me3','H3K4ac','H3K4me1','H3K4me2','H3K4me3','H3K79me1',
              'H3K79me2','H3K9ac','H3K9me3','H4K20me1','H4K8ac','H4K91ac',
              'POLR2A','RAD21','SMC3']

df38 = df38[df38['Target'].isin(sample_idx)]
df19 = df19[df19['Target'].isin(sample_idx)]
df38.index = df38['Target'] + '___' + df38['Biosample term name']
df19.index = df19['Target'] + '___' + df19['Biosample term name']

# match counts between assemblies
def match_and_extract(dfA, dfB):
    outA, outB = pd.DataFrame(), pd.DataFrame()
    for idx in dfB.index.unique():
        nB = (dfB.index == idx).sum()
        matches = dfA.loc[idx]
        nA = matches.shape[0]
        n = min(nA, nB)
        outA = pd.concat([outA, matches.head(n)])
        outB = pd.concat([outB, dfB.loc[idx].head(n)])
    return outA, outB

new_df38, new_df19 = match_and_extract(df38, df19)  # order doesn't matter for matching

# ---------------------------
# Load chromosome 6 data
# ---------------------------
def load_chr_df(genome, chr_id=6):
    fn = f"./results38/hg{genome}_chr{chr_id}_200data.h5"
    with h5py.File(fn, 'r') as hdf:
        grp = hdf['data']
        cols = np.array(grp['axis0']).astype(str)
        idxs = np.array(grp['axis1']).astype(str)
        vals = np.array(grp['block0_values'])
    df = pd.DataFrame(vals, index=idxs, columns=cols)
    # map back to Accession index
    meta = df38 if genome==38 else df19
    df.index = meta['Accession']
    return df

df38_chr6 = load_chr_df(38)
df19_chr6 = load_chr_df(19)

# ---------------------------
# Filter to jackknife inputs
# ---------------------------
filtered_df38_chr6 = filter_df_chr_by_accession(new_df38, df38_chr6)
filtered_df19_chr6 = filter_df_chr_by_accession(new_df19, df19_chr6)

# ---------------------------
# Jackknife for chunks 19 and 38
# ---------------------------
chunk_size = 100
for chunkid in (19, 38):
    # df38 jackknife
    jack38 = []
    inds38 = filtered_df38_chr6.index.tolist()
    start = chunkid * chunk_size
    for idx in inds38[start:start+chunk_size]:
        sub = filtered_df38_chr6.drop(idx)
        D = pairwise_distances(sub, metric='correlation')
        L = sch.linkage(squareform(D), method='complete')
        lbl = pd.merge(df38, sub, left_on='Accession', right_index=True)['Target'].tolist()
        mat, labs = compute_cooccurrence(sub, L, lbl)
        jack38.append({"left_out": idx, "matrix": mat, "labels": labs})
    
    rows = []
    for res in jack38:
        left_out = res["left_out"]
        mat      = res["matrix"]
        labels   = res["labels"]
        for i, li in enumerate(labels):
            for j, lj in enumerate(labels):
                rows.append({
                    "chunk":     chunkid,
                    "left_out":  left_out,
                    "label_i":   li,
                    "label_j":   lj,
                    "value":     mat[i, j]
                })
    df_out = pd.DataFrame(rows)
    df_out.to_csv(f"./data/jackknife_cooccurrence_38_chunk_{chunkid}.csv", index=False)
        
    # df19 jackknife
    jack19 = []
    inds19 = filtered_df19_chr6.index.tolist()
    for idx in inds19[start:start+chunk_size]:
        sub = filtered_df19_chr6.drop(idx)
        D = pairwise_distances(sub, metric='correlation')
        L = sch.linkage(squareform(D), method='complete')
        lbl = pd.merge(df19, sub, left_on='Accession', right_index=True)['Target'].tolist()
        mat, labs = compute_cooccurrence(sub, L, lbl)
        jack19.append({"left_out": idx, "matrix": mat, "labels": labs})
        
    rows = []
    for res in jack19:
        left_out = res["left_out"]
        mat      = res["matrix"]
        labels   = res["labels"]
        for i, li in enumerate(labels):
            for j, lj in enumerate(labels):
                rows.append({
                    "chunk":     chunkid,
                    "left_out":  left_out,
                    "label_i":   li,
                    "label_j":   lj,
                    "value":     mat[i, j]
                })
    df_out = pd.DataFrame(rows)
    df_out.to_csv(f"./data/jackknife_cooccurrence_19_chunk_{chunkid}.csv", index=False)
        
