# code to produce Figure 8 to get co-occurence plot, combines all chromosomes info at once
import numpy as np
import pandas as pd
import os
from collections import Counter
import matplotlib.pyplot as plt
import seaborn as sns

import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, ward, fcluster
from scipy.cluster import hierarchy
from scipy.stats import spearmanr
from sklearn.metrics.pairwise import pairwise_distances


def check_symmetric(arr):
    # Check if array is square
    if arr.shape[0] != arr.shape[1]:
        raise ValueError("The given array is not square!")

    # Check for symmetric elements
    non_symmetric_indices = []
    n = arr.shape[0]
    for i in range(n):
        for j in range(i+1, n):  # Only check the upper triangle
            if arr[i, j] != arr[j, i]:
                non_symmetric_indices.append((i, j))

    return non_symmetric_indices

def make_symmetric(mat):

    rows, cols = mat.shape
    for i in range(rows):
        for j in range(i + 1, cols):  # only consider upper triangular part
            if mat[i, j] != mat[j, i]:  # unsymmetrical
                symmetric_value = min(mat[i, j], mat[j, i])
                mat[i, j] = symmetric_value
                mat[j, i] = symmetric_value
    return mat


marks = pd.read_excel("./data/target_activity_factor.xlsx") 
marks.columns = ["Target", "Activity", "Factor"]
marks['Activity'].value_counts()
marks['Factor'].value_counts()

df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
df38['Target'].value_counts()
df38['Biosample term name'].value_counts()
merged_df = pd.merge(df38, marks, on='Target')


labels_mark = list(merged_df['Target'])
labels_cell = list(merged_df['Biosample term name'])
labels_factor = list(merged_df['Factor'])
labels_activity = list(merged_df['Activity'])


#get hierarchical clustering linkages fro all 23 chrs (23rd is X chr)
list_linkages = []
for i in range(1, 24):
    chr_id = i
    if chr_id == 23: #this is chromosome X
        chr_id = 'X'
    df_corr = pd.read_csv("./results38/hg38_chr" + str(chr_id) + "_200data" + 'correlation.h5', index_col=0)
    cor_dist = df_corr.to_numpy()
    np.fill_diagonal(cor_dist, 0)
    indices = check_symmetric(cor_dist)
    if len(indices) != 0:
        cor_dist = make_symmetric(cor_dist)
    condensed_dist = squareform(cor_dist)
    linkresult = sch.linkage(condensed_dist, method  = "complete")
    linkresult[linkresult < 0] = 0
    list_linkages.append(linkresult)



#parameters that can be changed
top_all_chrs = []
min_cluster_size = 4 #we don't consider clusters with less than 4 samples
N = 50 #visualize only top 50 co-occured items
unique_labels = list(np.unique(labels_mark)) #can be labels_cell, labels_factor, or labels_activity
co_occurrence_matrix = np.zeros((len(unique_labels), len(unique_labels)))
labels = labels_mark #can be labels_cell, labels_factor, or labels_activity
label_counts = Counter(labels)
for i in range(1, 24):
    chr_id = i 
    linkresult = list_linkages[i-1]
    clusters = fcluster(linkresult, 0.3, criterion='distance')    
    print(len(np.unique(clusters)))
    for cluster_id in np.unique(clusters):
        indices_in_cluster = np.where(clusters == cluster_id)[0]
        cluster_size = indices_in_cluster.shape[0]
        if cluster_size >= min_cluster_size:
            labels_in_cluster = [labels[index] for index in indices_in_cluster]
            unique_labels_in_cluster = list(np.unique(labels_in_cluster))
            for i, label1 in enumerate(unique_labels_in_cluster):
                for j, label2 in enumerate(unique_labels_in_cluster):
                    if i < j:
                        count_label1 = labels_in_cluster.count(label1)
                        count_label2 = labels_in_cluster.count(label2)
                        min_count = min(count_label1, count_label2)
                        co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)] += min_count/np.sqrt(label_counts[label1] * label_counts[label2])
                        co_occurrence_matrix[unique_labels.index(label2), unique_labels.index(label1)] = co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)]
                    if label1 == label2:
                        count_label1 = labels_in_cluster.count(label1)
                        count_label2 = labels_in_cluster.count(label2)
                        min_count = min(count_label1, count_label2)
                        co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)] += min_count/np.sqrt(label_counts[label1] * label_counts[label2])
                        co_occurrence_matrix[unique_labels.index(label2), unique_labels.index(label1)] = co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)]
    upper_triangle_indices = np.triu_indices_from(co_occurrence_matrix, 1)
    upper_triangle_values = co_occurrence_matrix[upper_triangle_indices]
    sorted_indices = np.argsort(upper_triangle_values)[::-1]
    N = len(np.unique(sorted_indices))
    top_indices = [upper_triangle_indices[i][sorted_indices[:N]] for i in (0, 1)]
    top_co_occurring_labels = [(unique_labels[i], unique_labels[j]) for i, j in zip(*top_indices)]
    sorted_data = sorted(top_co_occurring_labels, key=lambda x: (x[0], x[1]))
    top_all_chrs.append(sorted_data)

flattened_data = sorted_data
item_freqs = Counter(flattened_data)
item_freqs.most_common()

#annot=True will show the co-occurence counts computed above
plt.figure(figsize=(20, 20))
sns.heatmap(co_occurrence_matrix, annot=True, cmap="coolwarm", cbar=True, square=True, xticklabels=unique_labels, yticklabels=unique_labels)
plt.title("Co-occurrence Matrix Heatmap of epigenetic modifiers in clusters (GRCh38)")
plt.xlabel("Epigenetic modifiers")
plt.ylabel("Epigenetic modifiers")
plt.show()

co_occurrence_matrix = co_occurrence_matrix / 23 #normalize over 23 chromosomes
#reorder the plot for visualization purposes to get FIGURE 8
df = pd.DataFrame(co_occurrence_matrix, columns=unique_labels, index=unique_labels)

# Perform hierarchical clustering for visualization purposes to visualize most co-occured items together
linked = linkage(df, 'single')
df_ordered = df.iloc[leaves_list(linked), leaves_list(linked)]
plt.figure(figsize=(30, 30))
ax = sns.heatmap(df_ordered, cmap='coolwarm', linewidths=.5, annot=False, fmt=".1f")
ax.tick_params(labelsize=20)
plt.xlabel("Epigenetic modifiers", fontsize=28)
plt.ylabel("Epigenetic modifiers", fontsize=28)
plt.title('Co-occurrence Matrix Heatmap of Epigentic Modifiers in Clusters (GRCh38)', fontsize=30)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
plt.savefig('./results38/co_occurrence_matrix_heatmap_allchr_hg38.eps')
#plt.show()







