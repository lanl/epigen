#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# produce_hg38_consistency_plot.py
#
# Code to generate Figure 8:
#   split hg38 data in two matching subsets (same (modifier, cell) pairs)
#   and analyze the difference between these sets
#
# Input required:
#     ./data/genome_df38.csv
#     ./chr_files/hg38_chr6_200data.h5
#
# To change the chromosome of interest, make sure the h5 file is available
# and change the variable:
# chr_id = 6
#
# Output:
#     hg38onlymatchdfs.eps
# For more information, please see figure 8 in the paper.
#
#
# Time to run:
# Less than an hour.
#
import numpy as np
import pandas as pd
import h5py
import scipy.cluster.hierarchy as sch
from scipy.cluster.hierarchy import fcluster
from scipy.spatial.distance import squareform
from sklearn.metrics.pairwise import pairwise_distances
import matplotlib.pyplot as plt
import matplotlib.cm as cm


def cluster_entropy(cluster_id, clusters, assigned_labels):
    assigned_labels = np.array(assigned_labels)
    n = len(np.unique(assigned_labels))
    proportions = (
        np.array([np.sum(assigned_labels[np.where(clusters == cluster_id)[0]] == label) for label in np.unique(assigned_labels)])
        / np.where(clusters == cluster_id)[0].shape[0]
    )
    entropy_terms = [-p * np.log2(p) if p > 0 else 0 for p in proportions]
    entropy = sum(entropy_terms)
    return entropy / np.log2(n)


def weighted_average_entropy(assigned_labels, clusters):
    unique_clusters = np.unique(clusters)
    total_weight = len(assigned_labels)
    total_sum = sum(
        cluster_entropy(cluster_id, clusters, assigned_labels) * np.where(clusters == cluster_id)[0].shape[0]
        for cluster_id in unique_clusters
    )
    return total_sum / total_weight


def filter_df_chr_by_accession(new_df, df_chr):
    # get the list of Accession values from new_df
    accession_values = new_df["Accession"].tolist()
    # filter df_chr based on these Accession values
    filtered_df_chr = df_chr[df_chr.index.isin(accession_values)]
    return filtered_df_chr


def plot_for_single_label(labels, ax, linkresult_1, color, label_name):
    true_wa = []
    for dist in distances:
        clusters = fcluster(linkresult_1, dist, criterion="distance")
        true_wa.append(weighted_average_entropy(labels, clusters))
    # plot the true weighted average entropy
    ax.plot(distances, true_wa, "-o", label=label_name, color=color)


# This is for chr 6 only

marks = pd.read_excel("./data/target_activity_factor.xlsx")
marks.columns = ["Target", "Activity", "Factor"]
df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ["Accession", "Target", "Biosample term name", "Genome"]]
merged_df = pd.merge(df38, marks, on="Target")
df38 = merged_df
chr_id = 6
filename = "./chr_files/hg38_chr" + str(chr_id) + "_200data"

with h5py.File(filename + ".h5", "r") as hdf:
    data_group = hdf.get("data")
    axis0 = np.array(data_group.get("axis0")).astype(str)
    axis1 = np.array(data_group.get("axis1")).astype(str)
    block0_items = np.array(data_group.get("block0_items")).astype(str)
    block0_values = np.array(data_group.get("block0_values"))
df38_chr6 = pd.DataFrame(data=block0_values, index=axis1, columns=axis0)
df38_chr6.index = df38["Accession"]


df38.index = df38["Target"] + str("___") + df38["Biosample term name"]
np.unique(df38["Target"], return_counts=True)
np.unique(df38["Biosample term name"], return_counts=True)


# df38 is the dataframe we want to split in 2 matching
# (modifier, cell type) subsets where there are the same # of modifiers and cell types
subset_A_list = []
subset_B_list = []

# Group df38 by (Target, Biosample term name (cell type))
for (target, cell), group in df38.groupby(["Target", "Biosample term name"]):
    n_rows = len(group)
    n_split = n_rows // 2  # Number of rows for each subset

    if n_split > 0:
        # shuffle to randomize the selection
        group_shuffled = group.sample(frac=1, random_state=42)
        # split the group into two equal parts
        subset_A = group_shuffled.iloc[:n_split]
        subset_B = group_shuffled.iloc[n_split : n_split * 2]
        subset_A_list.append(subset_A)
        subset_B_list.append(subset_B)
    else:
        print(f"Skipping combination {(target, cell)} due to insufficient samples: {n_rows}")


matched_df38_A = pd.concat(subset_A_list)
matched_df38_B = pd.concat(subset_B_list)

# Verify that each combination has the same count in both datasets
print("Subset A counts:")
print(matched_df38_A.groupby(["Target", "Biosample term name"]).size())
print("\nSubset B counts:")
print(matched_df38_B.groupby(["Target", "Biosample term name"]).size())

# For matched_df38_A
print("Counts in matched_df38_A:")
print("Target counts:")
print(matched_df38_A["Target"].value_counts())
print("\nCell (Biosample term name) counts:")
print(matched_df38_A["Biosample term name"].value_counts())
print("\nActivity counts:")
print(matched_df38_A["Activity"].value_counts())
print("\nFactor counts:")
print(matched_df38_A["Factor"].value_counts())

# For matched_df38_B
print("\nCounts in matched_df38_B:")
print("Target counts:")
print(matched_df38_B["Target"].value_counts())
print("\nCell (Biosample term name) counts:")
print(matched_df38_B["Biosample term name"].value_counts())
print("\nActivity counts:")
print(matched_df38_B["Activity"].value_counts())
print("\nFactor counts:")
print(matched_df38_B["Factor"].value_counts())


chr6_matched_A = df38_chr6.loc[matched_df38_A["Accession"]]
chr6_matched_B = df38_chr6.loc[matched_df38_B["Accession"]]


chr6_matched_A = df38_chr6.loc[matched_df38_A["Accession"].unique()]
chr6_matched_B = df38_chr6.loc[matched_df38_B["Accession"].unique()]


def process_subset(subset_df, df_chr6, df38_full):
    #  subset_df: the subset of df38 (e.g. matched_df38_A or matched_df38_B)
    #  df_chr6: the df38_chr6 DataFrame of chr 6 samples with 'Accession' as its index.
    #  df38_full: the full df38 DataFrame used for merging to obtain labels.

    # filter by accession values
    filtered_df = filter_df_chr_by_accession(subset_df, df_chr6)

    # compute linkage
    cor_dist = pairwise_distances(filtered_df, metric="correlation")
    condensed_dist = squareform(cor_dist)
    linkage_result = sch.linkage(condensed_dist, method="complete")

    # extract labels – merging on 'Accession' from full df38 and the index of filtered_df
    sub_target = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")["Target"].tolist()
    sub_cell = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")[
        "Biosample term name"
    ].tolist()
    labels_activity = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")["Activity"].tolist()
    labels_factor = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")["Factor"].tolist()
    labels_mark = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")["Target"].tolist()
    labels_cell = pd.merge(df38_full, filtered_df, left_on="Accession", right_on=filtered_df.index, how="inner")[
        "Biosample term name"
    ].tolist()

    return {
        "filtered_df": filtered_df,
        "cor_dist": cor_dist,
        "condensed_dist": condensed_dist,
        "linkage": linkage_result,
        "sub_target": sub_target,
        "sub_cell": sub_cell,
        "labels_activity": labels_activity,
        "labels_factor": labels_factor,
        "labels_mark": labels_mark,
        "labels_cell": labels_cell,
    }


result_A = process_subset(matched_df38_A, df38_chr6, df38)
result_B = process_subset(matched_df38_B, df38_chr6, df38)

linkage_38_A = result_A["linkage"]
linkage_38_B = result_B["linkage"]

# For subset A
df38_sub_target_A = result_A["sub_target"]
df38_labels_activity_A = result_A["labels_activity"]
df38_labels_factor_A = result_A["labels_factor"]
df38_labels_mark_A = result_A["labels_mark"]
df38_labels_cell_A = result_A["labels_cell"]

# For subset B
df38_sub_target_B = result_B["sub_target"]
df38_labels_activity_B = result_B["labels_activity"]
df38_labels_factor_B = result_B["labels_factor"]
df38_labels_mark_B = result_B["labels_mark"]
df38_labels_cell_B = result_B["labels_cell"]

# print the label lists to verify:
print("Subset A Labels:")
print("Target:", df38_sub_target_A)
print("Activity:", df38_labels_activity_A)
print("Factor:", df38_labels_factor_A)
print("Mark:", df38_labels_mark_A)
print("Cell:", df38_labels_cell_A)

print("\nSubset B Labels:")
print("Target:", df38_sub_target_B)
print("Activity:", df38_labels_activity_B)
print("Factor:", df38_labels_factor_B)
print("Mark:", df38_labels_mark_B)
print("Cell:", df38_labels_cell_B)


# plot figure 8
colors = list(cm.tab20(np.linspace(0, 1, 20))) + [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
distances = np.arange(2, 0, -0.1)

fig, ax = plt.subplots(figsize=(15, 10))
font_size = 22

# Plot for Subset A (first matched dataset)
plot_for_single_label(df38_labels_factor_A, ax, linkage_38_A, colors[0], "Factor (Subset A)")
plot_for_single_label(df38_labels_activity_A, ax, linkage_38_A, colors[1], "Activity (Subset A)")
plot_for_single_label(df38_labels_mark_A, ax, linkage_38_A, colors[2], "Modifier (Subset A)")
plot_for_single_label(df38_labels_cell_A, ax, linkage_38_A, colors[3], "Cell (Subset A)")

# Plot for Subset B (second matched dataset)
plot_for_single_label(df38_labels_factor_B, ax, linkage_38_B, colors[4], "Factor (Subset B)")
plot_for_single_label(df38_labels_activity_B, ax, linkage_38_B, colors[5], "Activity (Subset B)")
plot_for_single_label(df38_labels_mark_B, ax, linkage_38_B, colors[6], "Modifier (Subset B)")
plot_for_single_label(df38_labels_cell_B, ax, linkage_38_B, colors[7], "Cell (Subset B)")

ax.set_title("Weighted average of the normalized entropy for hg38 subsets", fontsize=30)
ax.set_xlabel("Distance", fontsize=font_size)
ax.set_ylabel("Entropy", fontsize=font_size)
ax.grid(True)
ax.tick_params(axis="both", which="major", labelsize=20)
ax.legend(fontsize=font_size)
plt.tight_layout()
plt.savefig("./hg38onlymatchdfs.eps", format="eps")
plt.show()
