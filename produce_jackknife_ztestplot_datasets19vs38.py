#Collects jackknife replicate files (600 total); here just read from already combined file: "./data/jackknife_cooccurrence.npz"
#Figure 14: heatmap of p-values (jackknife z-test). Test which co-occurrence patterns of pairs of modifiers differ between hg19 and hg38 datasets.

import os
import glob
import re
import numpy as np
import pandas as pd
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import seaborn as sns
import matplotlib.pyplot as plt

# 1. Gather files
files_38 = sorted(f for f in glob.glob("../epigen_data/jackknife/jackknife_cooccurrence_38_chunk_*_iter_*.csv")
                  if re.search(r"chunk_[0-5]_iter_\d+\.csv$", os.path.basename(f)))
files_19 = sorted(f for f in glob.glob("../epigen_data/jackknife/jackknife_cooccurrence_19_chunk_*_iter_*.csv")
                  if re.search(r"chunk_[0-5]_iter_\d+\.csv$", os.path.basename(f)))

# 2. Load jackknife matrices
def load_jackknife(files, labels):
    mats = []
    for f in files:
        df = pd.read_csv(f, index_col=0)
        df = df.loc[labels, labels]
        mats.append(df.values)
    return np.stack(mats)

# 3. Load first file to get label order
sample_df = pd.read_csv(files_38[0], index_col=0)
unique_labels = sample_df.index.tolist()

# 4. Load all matrices
jackknife_cooccurrence_38 = load_jackknife(files_38, unique_labels)
jackknife_cooccurrence_19 = load_jackknife(files_19, unique_labels)
# Save jackknife matrices
np.savez_compressed(
    "./data/jackknife_cooccurrence.npz",
    cooc38=jackknife_cooccurrence_38,
    cooc19=jackknife_cooccurrence_19,
    labels=unique_labels  # Optional: Save label order for reuse
)
print("Saved to ./data/jackknife_cooccurrence.npz")



# Load from .npz
data = np.load("./data/jackknife_cooccurrence.npz", allow_pickle=True)
jackknife_cooccurrence_38 = data["cooc38"]
jackknife_cooccurrence_19 = data["cooc19"]
unique_labels = data["labels"].tolist()  # Convert from array to list if needed
print("Loaded shapes:",
      "38:", jackknife_cooccurrence_38.shape,
      "19:", jackknife_cooccurrence_19.shape)

# 5. Compute statistics (z-test)
results = []
L = len(unique_labels)
for i in range(L):
    for j in range(i, L):
        rep38 = jackknife_cooccurrence_38[:, i, j]
        rep19 = jackknife_cooccurrence_19[:, i, j]
        mean38, mean19 = rep38.mean(), rep19.mean()
        std38, std19 = rep38.std(ddof=1), rep19.std(ddof=1)
        delta = mean38 - mean19
        se = np.sqrt(std38**2 + std19**2)
        if se == 0:
            if mean38 == mean19 == 1.0:
                p_value = 1.0
            elif mean38 == mean19 == 0.0:
                p_value = np.nan
            else:
                p_value = 1.0
            z_score = 0.0
        else:
            z_score = delta / se
            p_value = 2 * (1 - norm.cdf(abs(z_score)))
        results.append({
            'label_i': unique_labels[i],
            'label_j': unique_labels[j],
            'p_value': p_value
        })

df_results = pd.DataFrame(results)
# 6. Adjust p-values
mask_valid = df_results['p_value'].notna()
rej, p_adj, _, _ = multipletests(df_results.loc[mask_valid, 'p_value'], method='fdr_bh')
df_results['p_adjusted'] = np.nan
df_results.loc[mask_valid, 'p_adjusted'] = p_adj

# 7. Build p-value matrix
labels = sorted(set(df_results['label_i']).union(df_results['label_j']))
p_mat = pd.DataFrame(np.nan, index=labels, columns=labels)
for _, row in df_results.iterrows():
    i, j = row['label_i'], row['label_j']
    p = row['p_adjusted']
    p_mat.at[i, j] = p
    p_mat.at[j, i] = p

# 8. Plot heatmap, mask NaNs (correspond to 0 co-occurences)
mask = p_mat.isna()
cmap = sns.color_palette("coolwarm", as_cmap=True)
cmap.set_bad("lightgray")

plt.figure(figsize=(30, 30))
ax = sns.heatmap(
    p_mat, mask=mask, cmap=cmap, vmin=0, vmax=1, square=True,
    linewidths=2, cbar_kws={'label': 'FDR-adjusted p-value'}
)
ax.tick_params(labelsize=30)
plt.title("Jackknife mean-difference: adjusted p-values", fontsize=30)
plt.xlabel("Epigenetic modifiers", fontsize=28)
plt.ylabel("Epigenetic modifiers", fontsize=28)
plt.setp(ax.get_yticklabels(), rotation=0, ha='right')
ax.collections[0].colorbar.ax.tick_params(labelsize=28)
plt.savefig("jackknife_mean_difference_heatmap_adjusted_nan.png", dpi=150, bbox_inches='tight')
print("Heatmap produced")
plt.show()
