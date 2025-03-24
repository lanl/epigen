# get_newick_string.py
#
#
# The newick string is a format for describing trees.
# This format is used by phylogenetic software to describe trees.
# https://en.wikipedia.org/wiki/Newick_format
#
# The string is used for
#
# Input file needed: correlation file  for chromosome of interest.
# The file name will be something like:
#     ./results38/hg38_chr6_200datacorrelation.h5
#
#
# To change the chromosome of interest, change the
# variable chr_id = 6
#
# This file takes approximately 1-2hrs to run; it is not instantaneous.
#
import numpy as np
import pandas as pd
import scipy.cluster.hierarchy as sch
from scipy.spatial.distance import squareform


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


def get_peaks_matching_genes(binwidth, d, ret_df):
    start_p = list(map(int, d.index))
    start_df = pd.DataFrame()
    new_df = pd.DataFrame()
    Nlimit = len(start_p)
    for k in range(Nlimit):
        if new_df.shape[0] < Nlimit:
            gdf = ret_df[(ret_df['start'] >= start_p[k] * binwidth - setbp) & (ret_df['end'] <= start_p[k] * binwidth + binwidth + setbp)]
            if gdf.shape[0] != 0:
                for t in np.unique(gdf['GeneName']):
                    temp = gdf[gdf['GeneName'] == t]
                    data = {'GeneName': [t], 'start': np.min(temp['start']), 'end': np.max(temp['end'])}
                    new_df = pd.concat([start_df, pd.DataFrame(data)], axis=0).reset_index(drop=True)
                    start_df = new_df

            if (new_df.shape[0] > 0):
                new_df = new_df.sort_values('start')
                new_df.reset_index(inplace=True)
                del new_df['index']

    return new_df


def get_samples_from_cluster(cluster_idx, linkresult, n_samples):
    # recursively fetches all original sample indices for a given cluster index
    if cluster_idx < n_samples:
        return [cluster_idx]
    else:
        row_idx = cluster_idx - n_samples
        left = int(linkresult[row_idx, 0])
        right = int(linkresult[row_idx, 1])
        return get_samples_from_cluster(left, linkresult, n_samples) + get_samples_from_cluster(right, linkresult, n_samples)


def linkage_to_newick(node, labels, gene_counts, linkresult):
    if node < len(labels):
        return f"{labels[int(node)]}:{gene_counts.get(node, 1.0)}"
    else:
        left, right, _, _ = linkresult[int(node) - len(labels)]
        left_str = linkage_to_newick(left, labels, gene_counts, linkresult)
        right_str = linkage_to_newick(right, labels, gene_counts, linkresult)
        left_count = gene_counts.get(left, 0)
        right_count = gene_counts.get(right, 0)
        current_count = gene_counts.get(node, 0)
        if current_count > (left_count + right_count):
            print(f"Discrepancy found at node {node}: current={current_count}, left={left_count}, right={right_count}")
        return f"({left_str},{right_str}):{current_count}"


def get_newick_from_dendro(pv_df, linkresult, thrs, divisor, ret_df):
    gene_counts = {}
    for idx in range(len(pv_df)):
        subset_cluster = pv_df.iloc[[idx]]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(divisor, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[idx] = ret_dfi.shape[0]
    current_id = len(pv_df)
    for row in linkresult:
        left, right = int(row[0]), int(row[1])
        left_samples = get_samples_from_cluster(left, linkresult, len(pv_df))
        right_samples = get_samples_from_cluster(right, linkresult, len(pv_df))
        combined_samples = left_samples + right_samples
        subset_cluster = pv_df.iloc[combined_samples]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))  # Use min instead of max for all samples condition
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(divisor, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[current_id] = ret_dfi.shape[0]
        current_id = current_id + 1
    root_node = len(pv_df) + len(linkresult) - 1
    newick_str = linkage_to_newick(root_node, pv_df.index, gene_counts, linkresult)
    return gene_counts, newick_str


def get_gene_count_from_dendro(pv_df, linkresult, thrs, divisor, ret_df):
    gene_counts = {}
    current_id = 0
    for row in linkresult:
        left, right = int(row[0]), int(row[1])
        left_samples = get_samples_from_cluster(left, linkresult, len(pv_df))
        right_samples = get_samples_from_cluster(right, linkresult, len(pv_df))
        combined_samples = left_samples + right_samples
        subset_cluster = pv_df.iloc[combined_samples]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))  # Use min instead of max for all samples condition
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(divisor, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[current_id] = ret_dfi.shape[0]
        current_id = current_id + 1
    return gene_counts


def get_linkage_with_gene_counts_and_samples(pv_df, linkresult, thrs, divisor, ret_df):
    gene_counts = {}
    for idx in range(len(pv_df)):
        subset_cluster = pv_df.iloc[[idx]]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(divisor, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[idx] = ret_dfi.shape[0]

    samples_list = []

    for row_idx, row in enumerate(linkresult):
        left, right = int(row[0]), int(row[1])
        left_samples = get_samples_from_cluster(left, linkresult, len(pv_df))
        right_samples = get_samples_from_cluster(right, linkresult, len(pv_df))
        combined_samples = left_samples + right_samples

        # Append the samples to the samples list
        samples_list.append(','.join(map(str, combined_samples)))

        subset_cluster = pv_df.iloc[combined_samples]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(divisor, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[row_idx + len(pv_df)] = ret_dfi.shape[0]
                linkresult[row_idx, 2] = ret_dfi.shape[0]

    samples_array = np.array(samples_list)[:, np.newaxis]
    linkresult = np.hstack((linkresult, samples_array))

    return linkresult


def get_filtered_clusters(linkage_matrix, num_min_imp_genes, pv_df_length):
    # Convert the sample string into sets for each cluster
    sample_sets = [set(map(int, row[4].split(','))) for row in linkage_matrix]

    valid_clusters = {}

    # Reverse iterate through the linkage matrix
    for idx in range(len(linkage_matrix) - 1, -1, -1):
        row = linkage_matrix[idx]
        gene_count = int(float(row[2]))

        # If gene count is >= 6 or the cluster is of size 2 (smallest possible size)
        if gene_count >= 6 or len(sample_sets[idx]) == 2:
            current_samples = sample_sets[idx]

            # Check if any of the samples in current cluster already have a valid cluster
            overlapping_samples = {sample for sample in current_samples if sample in valid_clusters}

            if overlapping_samples:
                # Check if the current cluster is larger than the stored valid cluster for those overlapping samples
                for sample in overlapping_samples:
                    if len(current_samples) > len(valid_clusters[sample]):
                        valid_clusters[sample] = current_samples
            else:
                # If no overlapping samples store the current cluster for all its samples
                for sample in current_samples:
                    valid_clusters[sample] = current_samples

    # Convert valid clusters from sets to original format (comma-separated string)
    unique_clusters = {",".join(map(str, sorted(cluster))): cluster for cluster in valid_clusters.values()}.keys()

    # Get rows from the linkage matrix that match the unique clusters
    filtered_clusters = [linkage_matrix[idx] for idx, samples in enumerate(sample_sets) if
                         ",".join(map(str, sorted(samples))) in unique_clusters]

    all_original_leaves = set(range(int(float(linkage_matrix[-1][0])) + 1))
    classified_leaves = {int(leaf) for cluster in unique_clusters for leaf in cluster.split(',')}

    # Add missing leaves as individual clusters
    for missing_leaf in all_original_leaves - classified_leaves:
        filtered_clusters.append([missing_leaf, missing_leaf, 0, 1, str(missing_leaf)])

    cleaned_clusters = []
    for row in filtered_clusters:
        samples = list(map(int, str(row[4]).split(',')))
        if (float(row[2]) != 0 or max(samples) < pv_df_length) and max(samples) < pv_df_length:
            cleaned_clusters.append(row)

    return sorted(cleaned_clusters, key=lambda x: (len(x[4].split(',')), float(x[0]), float(x[1])))


if __name__ == "__main__":

    # INPUT CHROMOSOME NUMBER
    chr_id = 6  # for chromosome 6; could be 1-22 or X #
    # Significance threshold
    thrs = 0.05
    # Threshold for identifying relevance to genes (buffer size)
    setbp = 500
    bin_width = 200  # length of bins in base pairs
    df = pd.read_csv("./data/genome_df38.csv", delimiter=",")
    df = pd.DataFrame(df)
    df38 = df.loc[:, ['Accession', 'Target', 'Biosample term name', 'Genome']]

    # We defined region of the genome as important in a certain cluster if its p-value was less than 0.05 in all samples
    # within that cluster. For every 200 bp region marked as important, we expanded the region by adding a 500 bp buffer to each
    # region. These expanded regions of genome were then matched with known gene positions from the Genome Browser, allowing
    # us to associate certain regions with specific genes.

    # read file with features where at least one sample has p-value <= 0.05
    pv_df = pd.read_pickle("./results38/pvdf_" + str(chr_id) + ".pkl")
    pv_df.index = df38['Target']
    # load reference file which contains known genes and their location in hg38 assembly
    ret_df = pd.read_csv("./results38/chr" + str(chr_id) + '_ret_df.csv', sep="\t")
    # load correlation file and obtain linkage from which we can get a dendrogram
    df_corr = pd.read_csv("./results38/hg38_chr" + str(chr_id) + "_200data" + 'correlation.h5', index_col=0)
    cor_dist = df_corr.to_numpy()
    np.fill_diagonal(cor_dist, 0)
    indices = check_symmetric(cor_dist)
    if len(indices) != 0:
        cor_dist = make_symmetric(cor_dist)
    condensed_dist = squareform(cor_dist)
    linkresult = sch.linkage(condensed_dist, method="complete")
    linkresult[linkresult < 0] = 0

    # with tip branch lengths
    gene_counts = {}
    # tip nodes
    for idx in range(len(pv_df)):
        subset_cluster = pv_df.iloc[[idx]]
        d = pd.DataFrame(np.max(subset_cluster, axis=0))
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(bin_width, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[idx] = ret_dfi.shape[0]

    # all other nodes
    current_id = len(pv_df)
    for row in linkresult:
        left, right = int(row[0]), int(row[1])
        left_samples = get_samples_from_cluster(left, linkresult, len(pv_df))
        right_samples = get_samples_from_cluster(right, linkresult, len(pv_df))
        combined_samples = left_samples + right_samples
        subset_cluster = pv_df.iloc[combined_samples]
        print(subset_cluster)
        d = pd.DataFrame(np.max(subset_cluster, axis=0))  # Use min instead of max for all samples condition
        d.columns = ['max']
        d = d.loc[d['max'] <= thrs, 'max']
        if len(d) > 0:
            ret_dfi = get_peaks_matching_genes(bin_width, d, ret_df)
            if ret_dfi.shape[0] != 0:
                gene_counts[current_id] = ret_dfi.shape[0]
                print("#genes", ret_dfi.shape[0])
                print(ret_dfi['GeneName'])
        current_id = current_id + 1
    # Convert linkage matrix to Newick format using the function
    root_node = len(pv_df) + len(linkresult) - 1
    newick_str = linkage_to_newick(root_node, pv_df.index, gene_counts, linkresult)
    print(newick_str)  # save chr6 newick string, we will use it to plot circular dendrogram later
    with open('./results38/chr6_newick.txt', 'w') as file:
        file.write(newick_str)






