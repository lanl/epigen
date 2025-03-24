# produce_comparison_plots_datasets19vs38.py
#
# Code to compare hg19 and hg38 datasets, figures 12, 13, and 14.
#
# This file requires the minisom package:
#     https://github.com/JustGlowing/minisom?tab=readme-ov-file#installation
#
# Input:
#     ./data/genome_df38.csv
#     Binned data for all chromosomes:
#     ./chr_files/hg38_chr6_200data.h5
# To change the chromosome of interest (Fig 13 and 14), change the variable:
#     chr_id = 6
# and ensure that the relevant h5 file is available.
#
# Output:
#     600_entropy_paper_plot_hg38.eps
#     600_entropy_paper_plot_hg19.eps
#    ./results38/combined_entropy_paper_plot_hg19_hg38.eps
#    ./results38/marks38_600_co_occurrence_matrix_heatmap_38_N50_normalized.eps
#    ./results38/cells38_600_co_occurrence_matrix_heatmap_38_N50_normalized.eps
#    ./results38/marks19_600_co_occurrence_matrix_heatmap_19_N50_normalized.eps
#    SOM_visualization.eps
#    cells38_600.html
#    epi38_600.html
#    epi19_600.html
#    cells19_600.html
# For more information, see figures 12, 13, and 14 in the paper.
# Also included are UMAP plots (not shown in the paper)
#
# Time to run: One hour or less.
#
#
import numpy as np
import pandas as pd
from collections import Counter
from collections import OrderedDict
from ast import literal_eval
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform, pdist
from scipy.cluster.hierarchy import linkage, dendrogram, leaves_list, fcluster, ward
from sklearn.metrics.pairwise import pairwise_distances
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns
import plotly.offline as py
import plotly.graph_objs as go
from umap import UMAP


def check_symmetric(arr):
    if arr.shape[0] != arr.shape[1]:
        raise ValueError("The given array is not square!")
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
                # Take the minimum of the two unsymmetrical entries
                symmetric_value = min(mat[i, j], mat[j, i])
                mat[i, j] = symmetric_value
                mat[j, i] = symmetric_value
    return mat

def cluster_purity(cluster_id, clusters, assigned_labels):
    assigned_labels = np.array(assigned_labels)
    unique, counts = np.unique(assigned_labels[np.where(clusters == cluster_id)[0]], return_counts=True)
    return counts.max() / np.where(clusters == cluster_id)[0].shape[0]

def weighted_average_purity(assigned_labels, clusters):
    unique_clusters = np.unique(clusters)
    total_weight = len(assigned_labels)
    total_sum = sum(cluster_purity(cluster_id, clusters, assigned_labels) * np.where(clusters == cluster_id)[0].shape[0] for cluster_id in unique_clusters)
    return total_sum / total_weight

def cluster_entropy(cluster_id, clusters, assigned_labels):
    assigned_labels = np.array(assigned_labels)
    n = len(np.unique(assigned_labels))
    proportions = np.array([np.sum(assigned_labels[np.where(clusters == cluster_id)[0]] == label) for label in np.unique(assigned_labels)]) / np.where(clusters == cluster_id)[0].shape[0]
    entropy_terms = [-p * np.log2(p) if p > 0 else 0 for p in proportions]
    entropy = sum(entropy_terms)
    return entropy / np.log2(n)

def weighted_average_entropy(assigned_labels, clusters):
    unique_clusters = np.unique(clusters)
    total_weight = len(assigned_labels)
    total_sum = sum(cluster_entropy(cluster_id, clusters, assigned_labels) * np.where(clusters == cluster_id)[0].shape[0] for cluster_id in unique_clusters)
    return total_sum / total_weight


def filter_df_chr_by_accession(new_df, df_chr):
    # get the list of Accession values from new_df
    accession_values = new_df['Accession'].tolist()
    # filter df_chr based on these Accession values
    filtered_df_chr = df_chr[df_chr.index.isin(accession_values)]
    return filtered_df_chr

def match_and_extract(df19, df38): 
    new_df19 = pd.DataFrame()
    new_df38 = pd.DataFrame()
    # iterate through each unique index in df38
    for index in df38.index.unique():
        count_in_df38 = df38.loc[df38.index == index].shape[0]
        matching_samples_df19 = df19.loc[df19.index == index]

        # if df19 has fewer samples, adjust the count
        count_to_extract = min(count_in_df38, matching_samples_df19.shape[0])

        # extract the samples
        extracted_samples_df19 = matching_samples_df19.head(count_to_extract)
        extracted_samples_df38 = df38.loc[df38.index == index].head(count_to_extract)

        new_df19 = new_df19.append(extracted_samples_df19)
        new_df38 = new_df38.append(extracted_samples_df38)

    return new_df19, new_df38

def plot_for_labels(labels, ax, title, linkresult_1, chr_id, color):
    true_wa = []
    for dist in distances:
        clusters = fcluster(linkresult_1, dist, criterion='distance')
        true_wa.append(weighted_average_entropy(labels, clusters))
    # Plotting the true weighted average entropy
    ax.plot(distances, true_wa, '-o', label=f"Chr {chr_id}", color=color)
    ax.set_title(title)
    ax.set_xlabel('Distance: (1 - correlation coefficient)')
    ax.set_ylabel('Weighted average of the normalized entropy')

    ax.grid(True)
    
    return true_wa

def plot_for_single_label(labels, ax, linkresult_1, chr_id, color, label_name):
    true_wa = []
    for dist in distances:
        clusters = fcluster(linkresult_1, dist, criterion='distance')
        true_wa.append(weighted_average_entropy(labels, clusters))
    # Plotting the true weighted average entropy
    ax.plot(distances, true_wa, '-o', label=label_name, color=color)



marks = pd.read_excel("./data/target_activity_factor.xlsx") 
marks.columns = ["Target", "Activity", "Factor"]
df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
merged_df = pd.merge(df38, marks, on='Target')
df38 = merged_df
df19 = pd.read_csv("./data/genome_df19.csv", delimiter=",")
df19 = pd.DataFrame(df19)
df19 = df19.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
merged_df = pd.merge(df19, marks, on='Target')
df19 = merged_df

#get matching dfs where we have same amount of Targets (modifiers) and Cell types represented in both hg19 and hg38
sample_idx = ['CTCF', 'EP300', 'H2AFZ', 'H2AK5ac', 'H2BK120ac', 'H2BK12ac',
        'H2BK15ac', 'H2BK5ac', 'H3K14ac', 'H3K18ac', 'H3K23ac', 'H3K27ac',
        'H3K27me3', 'H3K36me3', 'H3K4ac', 'H3K4me1', 'H3K4me2', 'H3K4me3',
        'H3K79me1', 'H3K79me2', 'H3K9ac', 'H3K9me3', 'H4K20me1', 'H4K8ac',
        'H4K91ac', 'POLR2A', 'RAD21', 'SMC3'] #these modifiers present in both datasets

df38 = df38[df38['Target'].isin(sample_idx)]
df19 = df19[df19['Target'].isin(sample_idx)]
df38 = df38[df38['Biosample term name'].isin(df38['Biosample term name'])]
df19 = df19[df19['Biosample term name'].isin(df38['Biosample term name'])]
df38.index = df38['Target'] + str('___') + df38['Biosample term name'] 
df19.index = df19['Target'] + str('___') + df19['Biosample term name']

new_df19, new_df38 = match_and_extract(df19, df38)
np.unique(new_df38.index, return_counts=True)
np.unique(new_df19.index, return_counts=True)

np.unique(new_df38['Target'], return_counts=True)
np.unique(new_df19['Target'], return_counts=True)

np.unique(new_df38['Biosample term name'], return_counts=True)
np.unique(new_df19['Biosample term name'], return_counts=True)


#for chromosome 6, change on any other chromosome: 1-22 and X
chr_id = 6
filename = "./results38/hg38_chr" + str(chr_id) + "_200data"
df38_chr6 = pd.read_hdf(filename + ".h5", key='data')

filename = "./results38/hg19_chr" + str(chr_id) + "_200data"
df19_chr6 = pd.read_hdf(filename + ".h5", key='data')

df38 = pd.read_csv("./data/genome_df38.csv", delimiter=",")
df38 = pd.DataFrame(df38)
df38 = df38.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
df38_chr6.index = df38['Accession']

df19 = pd.read_csv("./data/genome_df19.csv", delimiter=",")
df19 = pd.DataFrame(df19)
df19 = df19.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]
df19_chr6.index = df19['Accession']


filtered_df38_chr6 = filter_df_chr_by_accession(new_df38, df38_chr6)
filtered_df19_chr6 = filter_df_chr_by_accession(new_df19, df19_chr6)


cor_dist_38 = pairwise_distances(filtered_df38_chr6, metric='correlation')
condensed_dist_38 = squareform(cor_dist_38)
linkage_38 = sch.linkage(condensed_dist_38, method  = "complete")

cor_dist_19 = pairwise_distances(filtered_df19_chr6, metric='correlation')
condensed_dist_19 = squareform(cor_dist_19)
linkage_19 = sch.linkage(condensed_dist_19, method  = "complete")


df38_sub_target = pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Target'].tolist()
df38_sub_cell = pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Biosample term name'].tolist()

df19_sub_target = pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Target'].tolist()
df19_sub_cell = pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Biosample term name'].tolist()


df38_labels_activity = pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Activity'].tolist()
df38_labels_factor= pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Factor'].tolist()
df38_labels_mark = pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Target'].tolist()
df38_labels_cell= pd.merge(df38, filtered_df38_chr6, left_on='Accession', right_on=filtered_df38_chr6.index, how='inner')['Biosample term name'].tolist()

df19_labels_activity = pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Activity'].tolist()
df19_labels_factor= pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Factor'].tolist()
df19_labels_mark = pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Target'].tolist()
df19_labels_cell= pd.merge(df19, filtered_df19_chr6, left_on='Accession', right_on=filtered_df19_chr6.index, how='inner')['Biosample term name'].tolist()



# get entropy plots separately for each dataset
# plot for hg38
df38_plot_results = OrderedDict()
df38_plot_results['Activity'] = df38_labels_activity
df38_plot_results['Factor'] = df38_labels_factor
df38_plot_results['Modifier'] = df38_labels_mark
df38_plot_results['Cell'] = df38_labels_cell


colors = list(cm.tab20(np.linspace(0, 1, 20)))
additional_colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  
colors += additional_colors
distances = np.arange(2, 0, -0.1)
del cor_dist_38['Unnamed: 0']
cor_dist_1 = cor_dist_38.to_numpy()
np.fill_diagonal(cor_dist_1, 0)
indices = check_symmetric(cor_dist_1)
if len(indices) != 0:
    cor_dist_1 = make_symmetric(cor_dist_1)
condensed_dist_1 = squareform(cor_dist_1)
linkresult_1 = sch.linkage(condensed_dist_1, method  = "complete")
linkresult_1[linkresult_1 < 0] = 0

colors = list(cm.tab20(np.linspace(0, 1, 20))) + [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

fig, ax = plt.subplots(figsize=(15, 10))
font_size = 22 

# plotting for each label on the same axes
plot_for_single_label(df38_labels_factor, ax, linkresult_1, 1, colors[0], 'Factor')
plot_for_single_label(df38_labels_activity, ax, linkresult_1, 1, colors[1], 'Activity')
plot_for_single_label(df38_labels_mark, ax, linkresult_1, 1, colors[2], 'Modifier')
plot_for_single_label(df38_labels_cell, ax, linkresult_1, 1, colors[3], 'Cell')

ax.set_title('GRCh38 - hg38 - Genome - Assembly', fontsize=30)
ax.set_xlabel('Distance: (1 - correlation coefficient)', fontsize=font_size)
ax.set_ylabel('Weighted average of the normalized entropy', fontsize=font_size)

ax.grid(True)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.legend(fontsize=font_size)
plt.tight_layout()
#plt.show()
plt.savefig("./results38/600_entropy_paper_plot_hg38.eps", format='eps')


# same plotting as above but for hg19
df19_plot_results = OrderedDict()
df19_plot_results['Activity'] = df19_labels_activity
df19_plot_results['Factor'] = df19_labels_factor
df19_plot_results['Modifier'] = df19_labels_mark
df19_plot_results['Cell'] = df19_labels_cell

colors = list(cm.tab20(np.linspace(0, 1, 20)))
additional_colors = [(1, 0, 0), (0, 1, 0), (0, 0, 1)]  
colors += additional_colors
distances = np.arange(2, 0, -0.1)
del cor_dist_19['Unnamed: 0']
cor_dist_1 = cor_dist_19.to_numpy()
np.fill_diagonal(cor_dist_1, 0)
indices = check_symmetric(cor_dist_1)
if len(indices) != 0:
    cor_dist_1 = make_symmetric(cor_dist_1)
condensed_dist_1 = squareform(cor_dist_1)
linkresult_1 = sch.linkage(condensed_dist_1, method  = "complete")
linkresult_1[linkresult_1 < 0] = 0

colors = list(cm.tab20(np.linspace(0, 1, 20))) + [(1, 0, 0), (0, 1, 0), (0, 0, 1)]

fig, ax = plt.subplots(figsize=(15, 10))
font_size = 22 

plot_for_single_label(df19_labels_factor, ax, linkresult_1, 1, colors[0], 'Factor')
plot_for_single_label(df19_labels_activity, ax, linkresult_1, 1, colors[1], 'Activity')
plot_for_single_label(df19_labels_mark, ax, linkresult_1, 1, colors[2], 'Modifier')
plot_for_single_label(df19_labels_cell, ax, linkresult_1, 1, colors[3], 'Cell')

ax.set_title('GRCh37 - hg19 - Genome - Assembly', fontsize=30)
ax.set_xlabel('Distance: (1 - correlation coefficient)', fontsize=font_size)
ax.set_ylabel('Weighted average of the normalized entropy', fontsize=font_size)

ax.grid(True)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.legend(fontsize=font_size)
plt.tight_layout()
#plt.show()
plt.savefig("./results38/600_entropy_paper_plot_hg19.eps", format='eps')


# get FIGURE 12: entropy plots
del cor_dist_38['Unnamed: 0']
cor_dist_38 = cor_dist_38.to_numpy()
np.fill_diagonal(cor_dist_38, 0)
indices = check_symmetric(cor_dist_38)
if len(indices) != 0:
    cor_dist_38 = make_symmetric(cor_dist_38)
condensed_dist_38 = squareform(cor_dist_38)
linkresult_38 = sch.linkage(condensed_dist_38, method  = "complete")
linkresult_38[linkresult_38 < 0] = 0


del cor_dist_19['Unnamed: 0']
cor_dist_19 = cor_dist_19.to_numpy()
np.fill_diagonal(cor_dist_19, 0)
indices = check_symmetric(cor_dist_19)
if len(indices) != 0:
    cor_dist_19 = make_symmetric(cor_dist_19)
condensed_dist_19 = squareform(cor_dist_19)
linkresult_19 = sch.linkage(condensed_dist_19, method  = "complete")
linkresult_19[linkresult_19 < 0] = 0


colors = list(cm.tab20(np.linspace(0, 1, 20))) + [(1, 0, 0), (0, 1, 0), (0, 0, 1)]
distances = np.arange(2, 0, -0.1)

fig, ax = plt.subplots(figsize=(15, 10))
font_size = 22 


plot_for_single_label(df19_labels_factor, ax, linkresult_19, colors[0], 'Factor (GRCh37)')
plot_for_single_label(df19_labels_activity, ax, linkresult_19, colors[1], 'Activity (GRCh37)')
plot_for_single_label(df19_labels_mark, ax, linkresult_19, colors[2], 'Modifier (GRCh37)')
plot_for_single_label(df19_labels_cell, ax, linkresult_19, colors[3], 'Cell (GRCh37)')
plot_for_single_label(df38_labels_factor, ax, linkresult_38, colors[4], 'Factor (GRCh38)')
plot_for_single_label(df38_labels_activity, ax, linkresult_38, colors[5], 'Activity (GRCh38)')
plot_for_single_label(df38_labels_mark, ax, linkresult_38, colors[6], 'Modifier (GRCh38)')
plot_for_single_label(df38_labels_cell, ax, linkresult_38, colors[7], 'Cell (GRCh38)')

ax.set_title('Weighted average of the normalized entropy', fontsize=30)
ax.set_xlabel('Distance', fontsize=font_size)
ax.set_ylabel('Entropy', fontsize=font_size)

ax.grid(True)
ax.tick_params(axis='both', which='major', labelsize=20)
ax.legend(fontsize=font_size)
plt.tight_layout()
plt.savefig("./results38/combined_entropy_paper_plot_hg19_hg38.eps", format='eps')
# plt.show()  



#get FIGURE 13 (2 images)
#hg38 for modifiers
top_all_chrs = []
min_cluster_size = 4
N = 50
unique_labels = list(np.unique(df38_labels_mark))
co_occurrence_matrix = np.zeros((len(unique_labels), len(unique_labels)))
labels = df38_labels_mark
label_counts = Counter(labels)
chr_id = 6
clusters = fcluster(linkresult_38, 0.3, criterion='distance')    
#print(len(np.unique(clusters)))
for cluster_id in np.unique(clusters):
    #print(np.where(clusters == cluster_id)[0])
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
                    co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)] += min_count/np.sqrt(label_counts[label1] * label_counts[label2])#/np.sqrt(count_label1*count_label2)
                    co_occurrence_matrix[unique_labels.index(label2), unique_labels.index(label1)] = co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)]
                if label1 == label2:
                    count_label1 = labels_in_cluster.count(label1)
                    count_label2 = labels_in_cluster.count(label2)
                    min_count = min(count_label1, count_label2)
                    co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)] += min_count/np.sqrt(label_counts[label1] * label_counts[label2])#/np.sqrt(count_label1*count_label2)
                    co_occurrence_matrix[unique_labels.index(label2), unique_labels.index(label1)] = co_occurrence_matrix[unique_labels.index(label1), unique_labels.index(label2)]



    upper_triangle_indices = np.triu_indices_from(co_occurrence_matrix, 1)
    upper_triangle_values = co_occurrence_matrix[upper_triangle_indices]
    sorted_indices = np.argsort(upper_triangle_values)[::-1]
    N = len(np.unique(sorted_indices))
    top_indices = [upper_triangle_indices[i][sorted_indices[:N]] for i in (0, 1)]
    top_co_occurring_labels = [(unique_labels[i], unique_labels[j]) for i, j in zip(*top_indices)]
    #print(top_co_occurring_labels)
    sorted_data = sorted(top_co_occurring_labels, key=lambda x: (x[0], x[1]))
    top_all_chrs.append(sorted_data)
    
    
flattened_data = sorted_data
item_freqs = Counter(flattened_data)
item_freqs.most_common()
plt.figure(figsize=(30, 30))
ax = sns.heatmap(co_occurrence_matrix, annot=False, cmap="coolwarm", cbar=True, square=True, xticklabels=unique_labels, yticklabels=unique_labels)
ax.tick_params(labelsize=30)
plt.xlabel("Epigenetic modifiers", fontsize=32)
plt.ylabel("Epigenetic modifiers", fontsize=32)
plt.title('Co-occurrence matrix heatmap of epigentic modifiers in clusters (GRCh38)', fontsize=35)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=30)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
#file_path = './results38/omarks38_600_co_occurrence_matrix_heatmap_38_N50.eps'
#plt.savefig(file_path, format='eps')
plt.show()

#reorder 
data = co_occurrence_matrix
df = pd.DataFrame(data, columns=unique_labels, index=unique_labels)

df38_ = pd.DataFrame(co_occurrence_matrix)
df38_.index = df.index

linked = linkage(df, 'single')
df_ordered = df.iloc[leaves_list(linked), leaves_list(linked)]
plt.figure(figsize=(30, 30))
ax = sns.heatmap(df_ordered, cmap='coolwarm', linewidths=.5, annot=False, fmt=".1f")
ax.tick_params(labelsize=20)
plt.xlabel("Epigenetic modifiers", fontsize=28)
plt.ylabel("Epigenetic modifiers", fontsize=28)
plt.title('Co-occurrence matrix heatmap of epigentic modifiers in clusters', fontsize=30)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
file_path = './results38/marks38_600_co_occurrence_matrix_heatmap_38_N50_normalized.eps'
plt.savefig(file_path, format='eps')
plt.show()


#hg38 for cell types (not shown in the paper)
top_all_chrs = []
min_cluster_size = 4
N = 50
unique_labels = list(np.unique(df38_labels_cell))
co_occurrence_matrix = np.zeros((len(unique_labels), len(unique_labels)))
labels = df38_labels_cell
label_counts = Counter(labels)
chr_id = 6
clusters = fcluster(linkresult_38, 0.3, criterion='distance')    
#print(len(np.unique(clusters)))
for cluster_id in np.unique(clusters):
    #print(np.where(clusters == cluster_id)[0])
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
    #print(top_co_occurring_labels)
    sorted_data = sorted(top_co_occurring_labels, key=lambda x: (x[0], x[1]))
    top_all_chrs.append(sorted_data)
flattened_data = sorted_data
item_freqs = Counter(flattened_data)
item_freqs.most_common()
plt.figure(figsize=(20, 20))
sns.heatmap(co_occurrence_matrix, annot=True, cmap="coolwarm", cbar=True, square=True, xticklabels=unique_labels, yticklabels=unique_labels)
plt.title("Co-occurrence Matrix Heatmap of Cell Types in Clusters")
plt.xlabel("Cell Types")
plt.ylabel("Cell Types")
plt.show()

#reorder 
data = co_occurrence_matrix
df = pd.DataFrame(data, columns=unique_labels, index=unique_labels)
linked = linkage(df, 'single')
df_ordered = df.iloc[leaves_list(linked), leaves_list(linked)]
plt.figure(figsize=(70, 70))
ax = sns.heatmap(df_ordered, cmap='coolwarm', linewidths=.5, annot=False, fmt=".1f")
ax.tick_params(labelsize=18)
plt.xlabel("Cell Types", fontsize=28)
plt.ylabel("Cell Types", fontsize=28)
plt.title('Co-occurrence Matrix Heatmap of Cell Types in Clusters', fontsize=30)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
file_path = './results38/cells38_600_co_occurrence_matrix_heatmap_38_N50_normalized.eps'
plt.savefig(file_path, format='eps')
plt.show()


#hg19 for modifiers
top_all_chrs = []
min_cluster_size = 4
N = 50
unique_labels = list(np.unique(df19_labels_mark))
co_occurrence_matrix = np.zeros((len(unique_labels), len(unique_labels)))
labels = df19_labels_mark
label_counts = Counter(labels)

chr_id = 6
clusters = fcluster(linkresult_19, 0.3, criterion='distance')    
#print(len(np.unique(clusters)))
for cluster_id in np.unique(clusters):
    #print(np.where(clusters == cluster_id)[0])
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
    #print(top_co_occurring_labels)
    sorted_data = sorted(top_co_occurring_labels, key=lambda x: (x[0], x[1]))
    top_all_chrs.append(sorted_data)
flattened_data = sorted_data
item_freqs = Counter(flattened_data)
item_freqs.most_common()
plt.figure(figsize=(30, 30))
ax = sns.heatmap(co_occurrence_matrix, annot=False, cmap="coolwarm", cbar=True, square=True, xticklabels=unique_labels, yticklabels=unique_labels)
ax.tick_params(labelsize=30)
plt.xlabel("Epigenetic modifiers", fontsize=32)
plt.ylabel("Epigenetic modifiers", fontsize=32)
plt.title('Co-occurrence matrix heatmap of epigentic modifiers in clusters (GRCh37)', fontsize=35)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=30)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
#file_path = './results38/omarks19_600_co_occurrence_matrix_heatmap_19_N50_normalized.eps'
#plt.savefig(file_path, format='eps')
#plt.show()

#reorder 
data = co_occurrence_matrix
df = pd.DataFrame(data, columns=unique_labels, index=unique_labels)

df19_ = pd.DataFrame(co_occurrence_matrix)
df19_.index = df.index

linked = linkage(df, 'single')
df_ordered = df.iloc[leaves_list(linked), leaves_list(linked)]
plt.figure(figsize=(30, 30))
ax = sns.heatmap(df_ordered, cmap='coolwarm', linewidths=.5, annot=False, fmt=".1f")
ax.tick_params(labelsize=20)
plt.xlabel("Epigenetic modifiers", fontsize=28)
plt.ylabel("Epigenetic modifiers", fontsize=28)
plt.title('Co-occurrence matrix heatmap of epigentic modifiers in clusters', fontsize=30)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=20)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
file_path = './results38/marks19_600_co_occurrence_matrix_heatmap_19_N50_normalized.eps'
plt.savefig(file_path, format='eps')
plt.show()

#hg19 for cell types (not shown in the paper)
top_all_chrs = []
min_cluster_size = 4
N = 50
unique_labels = list(np.unique(df19_labels_cell))
co_occurrence_matrix = np.zeros((len(unique_labels), len(unique_labels)))
labels = df19_labels_cell
label_counts = Counter(labels)

chr_id = 6
clusters = fcluster(linkresult_19, 0.3, criterion='distance')    
#print(len(np.unique(clusters)))
for cluster_id in np.unique(clusters):
    #print(np.where(clusters == cluster_id)[0])
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
    # print(top_co_occurring_labels)
    sorted_data = sorted(top_co_occurring_labels, key=lambda x: (x[0], x[1]))
    top_all_chrs.append(sorted_data)
flattened_data = sorted_data
item_freqs = Counter(flattened_data)
item_freqs.most_common()
plt.figure(figsize=(20, 20))
sns.heatmap(co_occurrence_matrix, annot=True, cmap="coolwarm", cbar=True, square=True, xticklabels=unique_labels, yticklabels=unique_labels)
plt.title("Co-occurrence Matrix Heatmap of Cell Types in Clusters")
plt.xlabel("Cell Types")
plt.ylabel("Cell Types")
plt.show()

# reorder
data = co_occurrence_matrix
df = pd.DataFrame(data, columns=unique_labels, index=unique_labels)
linked = linkage(df, 'single')
df_ordered = df.iloc[leaves_list(linked), leaves_list(linked)]
plt.figure(figsize=(70, 70))
ax = sns.heatmap(df_ordered, cmap='coolwarm', linewidths=.5, annot=False, fmt=".1f")
ax.tick_params(labelsize=18)
plt.xlabel("Cell Types", fontsize=28)
plt.ylabel("Cell Types", fontsize=28)
plt.title('Co-occurrence Matrix Heatmap of Cell Types in Clusters', fontsize=30)
cbar = ax.collections[0].colorbar
cbar.ax.tick_params(labelsize=18)
plt.setp(ax.get_yticklabels(), rotation=0, horizontalalignment='right')
file_path = './results38/cells19_600_co_occurrence_matrix_heatmap_19_N50_normalized.eps'
plt.savefig(file_path, format='eps')
plt.show()


from sklearn.manifold import MDS
nmds = MDS(n_components=2, dissimilarity="precomputed", random_state=42)
nmds_results = nmds.fit_transform(co_occurrence_matrix)
plt.figure(figsize=(8, 6))
plt.scatter(nmds_results[:, 0], nmds_results[:, 1], c='blue', alpha=0.6, edgecolors='k')
for i, label in enumerate(df.index):
    plt.text(nmds_results[i, 0], nmds_results[i, 1], label, fontsize=8, ha='right', color='darkblue')
plt.title("NMDS Visualization of co-occurences")
plt.xlabel("NMDS Dimension 1")
plt.ylabel("NMDS Dimension 2")
plt.grid(True)
plt.show()


# Figure 14
# co-occurence matrix and df.index
# Self-organizing maps
import numpy as np
from minisom import MiniSom
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon
from matplotlib import cm, colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable
from collections import defaultdict

som_x, som_y = 10, 10  # Adjust grid size based on your data
som = MiniSom(som_x, som_y, co_occurrence_matrix.shape[1],
              sigma=1.0, learning_rate=0.5, random_seed=42)
som.train_random(co_occurrence_matrix, 100000)

#  SOM visualization
xx, yy = som.get_euclidean_coordinates()
umatrix = som.distance_map()
weights = som.get_weights()
f = plt.figure(figsize=(12, 12))  # Increased figure size
ax = f.add_subplot(111)
ax.set_aspect('equal')

# Add hexagons 
for i in range(weights.shape[0]):
    for j in range(weights.shape[1]):
        wy = yy[(i, j)] * np.sqrt(3) / 2
        hex = RegularPolygon((xx[(i, j)], wy),
                             numVertices=6,
                             radius=.95 / np.sqrt(3),
                             facecolor=cm.Blues(umatrix[i, j]),
                             alpha=.4,
                             edgecolor='gray')
        ax.add_patch(hex)

cell_labels = defaultdict(list)
for idx, data_point in enumerate(co_occurrence_matrix):
    w = som.winner(data_point)  # Find BMU for each data point
    cell_labels[w].append(df.index[idx])

for (i, j), labels in cell_labels.items():
    wx, wy = som.convert_map_to_euclidean((i, j))
    wy = wy * np.sqrt(3) / 2

    num_labels = len(labels)
    fontsize = 9  # Adjust font size as needed
    line_height = 0.2  # Vertical spacing between labels

    # Calculate starting y-coordinate to center the labels vertically
    total_height = (num_labels - 1) * line_height
    start_y = wy + total_height / 2

    for idx, label in enumerate(labels):
        dy = -idx * line_height
        plt.text(wx, start_y + dy, label, fontsize=fontsize,
                 ha='center', va='center', color='darkblue')

padding = 2  # Adjust padding if labels are cut off
x_min, x_max = xx.min() - padding, xx.max() + padding
y_min, y_max = (yy * np.sqrt(3) / 2).min() - padding, (yy * np.sqrt(3) / 2).max() + padding
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

xrange = np.arange(weights.shape[0])
yrange = np.arange(weights.shape[1])
plt.xticks(xrange - 0.5, xrange)
plt.yticks(yrange * np.sqrt(3) / 2, yrange)

divider = make_axes_locatable(plt.gca())
ax_cb = divider.new_horizontal(size="5%", pad=0.05)
cb1 = colorbar.ColorbarBase(ax_cb, cmap=cm.Blues, orientation='vertical', alpha=.4)
cb1.ax.get_yaxis().labelpad = 16
cb1.ax.set_ylabel('Distance from neurons in the neighborhood', rotation=270, fontsize=16)
plt.gcf().add_axes(ax_cb)

plt.savefig('SOM_visualization.eps', format='eps')
plt.title("SOM Visualization with Labels")
plt.show()



# diag_elements = np.diag(df19_.values).copy()
# diag_elements[diag_elements == 0] = 1
# normalization_matrix = np.outer(diag_elements, diag_elements)
# co_occurrence_matrix_normalized = df19_.values / np.sqrt(normalization_matrix)
# np.fill_diagonal(co_occurrence_matrix_normalized, 1)
# df19_2 = co_occurrence_matrix_normalized


# diag_elements = np.diag(df38_.values).copy()
# diag_elements[diag_elements == 0] = 1
# normalization_matrix = np.outer(diag_elements, diag_elements)
# co_occurrence_matrix_normalized = df38_.values / np.sqrt(normalization_matrix)
# np.fill_diagonal(co_occurrence_matrix_normalized, 1)
# df38_2 = co_occurrence_matrix_normalized
import numpy as np
import pandas as pd
from statsmodels.stats.multitest import multipletests

def permutation_test(differences, n_permutations=10000):
    T_obs = np.mean(differences)
    permuted_stats = np.array([
        np.mean(differences * np.random.choice([-1, 1], size=len(differences)))
        for _ in range(n_permutations)
    ])
    
    if T_obs > 0:
        p_value = np.sum(permuted_stats >= T_obs) / n_permutations
    else:
        p_value = np.sum(permuted_stats <= T_obs) / n_permutations
    return p_value

labels = df19_.index.tolist()
p_values = []

# Perform permutation test for each row
for label in labels:
    row19 = df19_.loc[label].values
    row38 = df38_.loc[label].values
    differences = row19 - row38
    p = permutation_test(differences)
    p_values.append(p)

# Benjamini-Hochberg correction for multiple testing
reject, pvals_corrected, _, _ = multipletests(p_values, alpha=0.05, method='fdr_bh')

non_significant_labels = [label for label, rej in zip(labels, reject) if not rej]
significant_labels = [label for label, rej in zip(labels, reject) if rej]
print(f"Non-significant labels (p > 0.05): {non_significant_labels}")
print(f"Significant labels (p ≤ 0.05): {significant_labels}")




#UMAP plots (not shown in the paper), play with the n_neighbors=10 parameter (10, 20, ..., 100)
#hg38 for modifiers
embedding_38 = UMAP(n_neighbors=10, n_components=3, metric="precomputed").fit_transform(cor_dist_38)
dfu = pd.DataFrame(embedding_38, columns=('x', 'y', 'z'))
dfu.index = df38_sub_target
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])

colors = [
    'rgb(255, 0, 0)', 
    'rgb(0, 255, 0)', 
    'rgb(0, 0, 255)', 
    'rgb(255, 255, 0)', 
    'rgb(0, 255, 255)', 
    'rgb(255, 0, 255)', 
    'rgb(192, 192, 192)', 
    'rgb(128, 128, 128)', 
    'rgb(128, 0, 0)', 
    'rgb(128, 128, 0)', 
    'rgb(0, 128, 0)', 
    'rgb(128, 0, 128)', 
    'rgb(0, 128, 128)', 
    'rgb(0, 0, 128)', 
    'rgb(255, 165, 0)', 
    'rgb(255, 20, 147)', 
    'rgb(75, 0, 130)', 
    'rgb(240, 128, 128)', 
    'rgb(255, 255, 224)', 
    'rgb(124, 252, 0)', 
    'rgb(173, 216, 230)', 
    'rgb(255, 105, 180)', 
    'rgb(255, 218, 185)', 
    'rgb(219, 112, 147)', 
    'rgb(245, 222, 179)', 
    'rgb(32, 178, 170)', 
    'rgb(255, 228, 225)', 
    'rgb(218, 165, 32)', 
    'rgb(95, 158, 160)', 
    'rgb(175, 238, 238)'
]
#colors = np.unique(colors)
#colors = random.choices(colors, k=len(np.unique(dfu.index)))

data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
    x = dfu[dfu["class"] == name]["x"],
    y = dfu[dfu["class"] == name]["y"],
    z = dfu[dfu["class"] == name]["z"],
    name = labels[no],
    mode = 'markers',
    marker = dict(
        size = 8,
        color = '#%02x%02x%02x' %literal_eval(colors[no][3:]),#'#%02x%02x%02x' % literal_eval(colors[no][3:]),
        line = dict(
            color = literal_eval(colors[no][3:]), #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            width = 0.5
            ),
        opacity = 1
        )
    )
    data_graph.append(graph)
    
layout = go.Layout(
    scene = dict(
        camera = dict(
            eye = dict(
            x = 0.5,
            y = 0.5,
            z = 0.5
            )
        )
    ),
    margin = dict(
        l = 0,
        r = 0,
        b = 0,
        t = 0
    )
)
fig = go.Figure(data = data_graph, layout = layout)
py.iplot(fig, filename = '3d-scatter');
fig.write_html("epi38_600.html")

#hg38 for cell types
embedding_38 = UMAP(n_neighbors=10, n_components=3, metric="precomputed").fit_transform(cor_dist_38)
dfu = pd.DataFrame(embedding_38, columns=('x', 'y', 'z'))
dfu.index = df38_sub_cell
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])

colors = [
    'rgb(255, 0, 0)', 
    'rgb(0, 255, 0)', 
    'rgb(0, 0, 255)', 
    'rgb(255, 255, 0)', 
    'rgb(0, 255, 255)', 
    'rgb(255, 0, 255)', 
    'rgb(192, 192, 192)', 
    'rgb(128, 128, 128)', 
    'rgb(128, 0, 0)', 
    'rgb(128, 128, 0)', 
    'rgb(0, 128, 0)', 
    'rgb(128, 0, 128)', 
    'rgb(0, 128, 128)', 
    'rgb(0, 0, 128)', 
    'rgb(255, 165, 0)', 
    'rgb(255, 20, 147)', 
    'rgb(75, 0, 130)', 
    'rgb(240, 128, 128)', 
    'rgb(255, 255, 224)', 
    'rgb(124, 252, 0)', 
    'rgb(173, 216, 230)', 
    'rgb(255, 105, 180)', 
    'rgb(255, 218, 185)', 
    'rgb(219, 112, 147)', 
    'rgb(245, 222, 179)', 
    'rgb(32, 178, 170)', 
    'rgb(255, 228, 225)', 
    'rgb(218, 165, 32)', 
    'rgb(95, 158, 160)', 
    'rgb(175, 238, 238)'
]
# colors = np.unique(colors)
# colors = random.choices(colors, k=len(np.unique(dfu.index)))

data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
    x = dfu[dfu["class"] == name]["x"],
    y = dfu[dfu["class"] == name]["y"],
    z = dfu[dfu["class"] == name]["z"],
    name = labels[no],
    mode = 'markers',
    marker = dict(
        size = 8,
        color = '#%02x%02x%02x' %literal_eval(colors[no][3:]),#'#%02x%02x%02x' % literal_eval(colors[no][3:]),
        line = dict(
            color = literal_eval(colors[no][3:]), #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            width = 0.5
            ),
        opacity = 1
        )
    )
    data_graph.append(graph)
    
layout = go.Layout(
    scene = dict(
        camera = dict(
            eye = dict(
            x = 0.5,
            y = 0.5,
            z = 0.5
            )
        )
    ),
    margin = dict(
        l = 0,
        r = 0,
        b = 0,
        t = 0
    )
)
fig = go.Figure(data = data_graph, layout = layout)
py.iplot(fig, filename = '3d-scatter');
fig.write_html("cells38_600.html")




# hg19 for modifiers
embedding_19 = UMAP(n_neighbors=10, n_components=3, metric="precomputed").fit_transform(cor_dist_19)
dfu = pd.DataFrame(embedding_19, columns=('x', 'y', 'z'))
dfu.index = df19_sub_target
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])

colors = [
    'rgb(255, 0, 0)', 
    'rgb(0, 255, 0)', 
    'rgb(0, 0, 255)', 
    'rgb(255, 255, 0)', 
    'rgb(0, 255, 255)', 
    'rgb(255, 0, 255)', 
    'rgb(192, 192, 192)', 
    'rgb(128, 128, 128)', 
    'rgb(128, 0, 0)', 
    'rgb(128, 128, 0)', 
    'rgb(0, 128, 0)', 
    'rgb(128, 0, 128)', 
    'rgb(0, 128, 128)', 
    'rgb(0, 0, 128)', 
    'rgb(255, 165, 0)', 
    'rgb(255, 20, 147)', 
    'rgb(75, 0, 130)', 
    'rgb(240, 128, 128)', 
    'rgb(255, 255, 224)', 
    'rgb(124, 252, 0)', 
    'rgb(173, 216, 230)', 
    'rgb(255, 105, 180)', 
    'rgb(255, 218, 185)', 
    'rgb(219, 112, 147)', 
    'rgb(245, 222, 179)', 
    'rgb(32, 178, 170)', 
    'rgb(255, 228, 225)', 
    'rgb(218, 165, 32)', 
    'rgb(95, 158, 160)', 
    'rgb(175, 238, 238)'
]
#colors = np.unique(colors)
#colors = random.choices(colors, k=len(np.unique(dfu.index)))

data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
    x = dfu[dfu["class"] == name]["x"],
    y = dfu[dfu["class"] == name]["y"],
    z = dfu[dfu["class"] == name]["z"],
    name = labels[no],
    mode = 'markers',
    marker = dict(
        size = 8,
        color = '#%02x%02x%02x' %literal_eval(colors[no][3:]),#'#%02x%02x%02x' % literal_eval(colors[no][3:]),
        line = dict(
            color = literal_eval(colors[no][3:]), #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            width = 0.5
            ),
        opacity = 1
        )
    )
    data_graph.append(graph)
    
layout = go.Layout(
    scene = dict(
        camera = dict(
            eye = dict(
            x = 0.5,
            y = 0.5,
            z = 0.5
            )
        )
    ),
    margin = dict(
        l = 0,
        r = 0,
        b = 0,
        t = 0
    )
)
fig = go.Figure(data = data_graph, layout = layout)
py.iplot(fig, filename = '3d-scatter');
fig.write_html("epi19_600.html")


#hg38 for cell types
embedding_19 = UMAP(n_neighbors=10, n_components=3, metric="precomputed").fit_transform(cor_dist_19)
dfu = pd.DataFrame(embedding_19, columns=('x', 'y', 'z'))
dfu.index = df19_sub_cell
dfu["class"] = dfu.index
labels = np.unique(dfu["class"])

colors = [
    'rgb(255, 0, 0)', 
    'rgb(0, 255, 0)', 
    'rgb(0, 0, 255)', 
    'rgb(255, 255, 0)', 
    'rgb(0, 255, 255)', 
    'rgb(255, 0, 255)', 
    'rgb(192, 192, 192)', 
    'rgb(128, 128, 128)', 
    'rgb(128, 0, 0)', 
    'rgb(128, 128, 0)', 
    'rgb(0, 128, 0)', 
    'rgb(128, 0, 128)', 
    'rgb(0, 128, 128)', 
    'rgb(0, 0, 128)', 
    'rgb(255, 165, 0)', 
    'rgb(255, 20, 147)', 
    'rgb(75, 0, 130)', 
    'rgb(240, 128, 128)', 
    'rgb(255, 255, 224)', 
    'rgb(124, 252, 0)', 
    'rgb(173, 216, 230)', 
    'rgb(255, 105, 180)', 
    'rgb(255, 218, 185)', 
    'rgb(219, 112, 147)', 
    'rgb(245, 222, 179)', 
    'rgb(32, 178, 170)', 
    'rgb(255, 228, 225)', 
    'rgb(218, 165, 32)', 
    'rgb(95, 158, 160)', 
    'rgb(175, 238, 238)'
]
#colors = np.unique(colors)
#colors = random.choices(colors, k=len(np.unique(dfu.index)))

data_graph = []
for no, name in enumerate(np.unique(dfu["class"])):
    graph = go.Scatter3d(
    x = dfu[dfu["class"] == name]["x"],
    y = dfu[dfu["class"] == name]["y"],
    z = dfu[dfu["class"] == name]["z"],
    name = labels[no],
    mode = 'markers',
    marker = dict(
        size = 8,
        color = '#%02x%02x%02x' %literal_eval(colors[no][3:]),#'#%02x%02x%02x' % literal_eval(colors[no][3:]),
        line = dict(
            color = literal_eval(colors[no][3:]), #            color = '#%02x%02x%02x' % literal_eval(colors[no][3:]),
            width = 0.5
            ),
        opacity = 1
        )
    )
    data_graph.append(graph)
    
layout = go.Layout(
    scene = dict(
        camera = dict(
            eye = dict(
            x = 0.5,
            y = 0.5,
            z = 0.5
            )
        )
    ),
    margin = dict(
        l = 0,
        r = 0,
        b = 0,
        t = 0
    )
)
fig = go.Figure(data = data_graph, layout = layout)
py.iplot(fig, filename = '3d-scatter');
fig.write_html("cells19_600.html")


