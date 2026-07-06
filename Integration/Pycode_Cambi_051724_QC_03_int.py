import os
import glob
import warnings
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import sys
import importlib

#sc.settings.verbosity = 3
warnings.simplefilter(action='ignore', category=FutureWarning)

script_dirs = ['E:/multiomatic/code/']
for d in script_dirs:
    path = os.path.abspath(d)
    if path not in sys.path:
        sys.path.append(path)

from Pycode_iter_clust_scvi import iter_clust_scvi
from Pycode_scanpy_helper_funtions import downsample_by_clusters

adata = sc.read_h5ad('F:/Colab/Cambi_lab/Integration/QC_02/Cambi_051724_QC_02_int_filt.h5ad')
adata.layers["counts"] = adata.X.copy()

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

adata = iter_clust_scvi(
    adata,
    categorical_covariates=['sample'],
    continuous_covariates=['n_genes_by_counts', 'pct_counts_mt'],
    n_latent=10,
    num_iterations=2
)

sc.tl.umap(adata)

scVI_embed = pd.DataFrame(adata.obsm['X_scVI'], index=adata.obs.index)
umap_embed = pd.DataFrame(adata.obsm['X_umap'], index=adata.obs.index)

scVI_embed.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_scVI_embed.csv')
umap_embed.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_umap_embed.csv')

adata.write_h5ad("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int.h5ad", compression="gzip")
adata.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_obs.csv')

print(len(adata.obs["level_1"].unique()))
print(len(adata.obs["level_2"].unique()))

#####
#####

import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
import IPython
IPython.display.set_matplotlib_formats = set_matplotlib_formats

## Vis
#adata = sc.read_h5ad('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int.h5ad')
os.chdir('F:/Colab/Cambi_lab/Integration/QC_03')
sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=5, figsize=(5, 5))

def plot_umap(adata, obs_key, filename_suffix):
    n_cls = len(adata.obs[obs_key].unique())
    colormap = plt.get_cmap('tab20', n_cls)
    colors = [colormap(i) for i in range(n_cls)]
    sc.pl.umap(adata, color=obs_key, legend_loc='on data', palette=colors, save=f'_{filename_suffix}.pdf')

plot_umap(adata, 'level_1', 'Cambi_051724_QC_03_int_level_1')
plot_umap(adata, 'level_2', 'Cambi_051724_QC_03_int_level_2')
plot_umap(adata, 'class_name', 'Cambi_051724_QC_03_int_class_name')
plot_umap(adata, 'subclass_name', 'Cambi_051724_QC_03_int_subclass_name')

## Umap by conditions
sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=12, figsize=(5, 5))
conditions = adata.obs['genotype'].unique()
num_cols = 2
num_rows = (len(conditions) + num_cols - 1) // num_cols
fig, axs = plt.subplots(num_rows, num_cols, figsize=(5 * num_cols, 5 * num_rows))

n_cls = len(adata.obs['level_1'].unique())
colormap = plt.get_cmap('tab20', n_cls)
colors = [colormap(i) for i in range(n_cls)]

for i, cond in enumerate(conditions):
    ax = axs[i // num_cols, i % num_cols] if num_rows > 1 else axs[i % num_cols]
    adata_subset = adata[adata.obs['genotype'] == cond].copy()
    sc.pl.umap(adata_subset, color='level_1', palette=colors, ax=ax, show=False, legend_loc='on data')
    ax.set_title(cond)

for i in range(len(conditions), num_cols * num_rows):
    fig.delaxes(axs[i // num_cols, i % num_cols] if num_rows > 1 else axs[i % num_cols])

plt.tight_layout()
plt.savefig('figures/umap_Cambi_051724_QC_03_int_level_1_by_genotype.pdf', dpi=300)
plt.show()

sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=5, figsize=(5, 5))
sc.pl.umap(adata, color="total_counts", save='_Cambi_051724_QC_03_int_n_genes_by_counts.pdf')
sc.pl.umap(adata, color="pct_counts_mt", save='_Cambi_051724_QC_03_int_pct_counts_mt.pdf')
sc.pl.umap(adata, color="scDblFinder.score", save='_Cambi_051724_QC_03_int_scDblFinder.score.pdf')

##########
##########

## Compile level_2 cluster stats
cl_stats = adata.obs.groupby('level_2', observed=False)[
    ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'scDblFinder.score']
].mean()
cl_stats['cell_count'] = adata.obs.groupby('level_2', observed=False).size()

pct_counts_mt_99th = cl_stats['pct_counts_mt'].quantile(0.99)
cl_stats['mt_flag'] = cl_stats['pct_counts_mt'] > pct_counts_mt_99th
cl_stats['dbl_flag'] = cl_stats['scDblFinder.score'] > 0.1
cl_stats['cl_flag'] = np.select(
    [cl_stats['scDblFinder.score'] > 0.1, cl_stats['total_counts'] < 1000],
    ['doublet', 'low_quality'],
    default='singlet'
)

def top_annotations(series):
    counts = series.value_counts(normalize=True) * 100
    counts = counts[counts > 0]
    return ', '.join([f"{label}: {perc:.2f}%" for label, perc in counts.head(3).items()])

for col in ['class_name', 'subclass_name', 'supertype_name']:
    cl_stats[col] = (
        adata.obs.groupby('level_2', observed=False)[col]
        .apply(top_annotations))

cl_stats.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_level_2_cluster_stats.csv')

#####
#####

## marker genes
cluster_key = 'level_2'
adata_sub = downsample_by_clusters(adata, cluster_key=cluster_key, ncells_per_cluster=200)

warnings.filterwarnings("ignore", category=pd.errors.PerformanceWarning)
sc.tl.rank_genes_groups(adata_sub, groupby=cluster_key, method='wilcoxon', pts=False)
de_genes = sc.get.rank_genes_groups_df(adata_sub, group=None).copy()
de_genes = de_genes.query("pvals_adj < 0.05 and logfoldchanges > 0.25").copy()
de_genes.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_level_2_marker_genes.csv')

top_df = (
    de_genes.groupby('group')[['names', 'scores', 'logfoldchanges']]
    .apply(lambda g: pd.Series({
        'top_score_genes': ', '.join(g.nlargest(10, 'scores')['names']),
        'top_logfc_genes': ', '.join(g.nlargest(10, 'logfoldchanges')['names'])
    })))

cl_stats['top_score_genes'] = cl_stats.index.map(lambda c: top_df['top_score_genes'].get(c, ''))
cl_stats['top_logfc_genes'] = cl_stats.index.map(lambda c: top_df['top_logfc_genes'].get(c, ''))
cl_stats.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_level_2_cluster_stats_combined.csv')

#####
#####

## cl_flag
adata.obs['cl_flag'] = adata.obs['level_2'].map(cl_stats['cl_flag'])
adata.obs["cl_flag"].value_counts()

sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=4, figsize=(5, 5))
sc.pl.umap(adata, color='cl_flag', save= '_Cambi_051724_QC_03_int_cl_flag.pdf')

sc.tl.leiden(adata, resolution=2, key_added="level_1_2")
sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=5, figsize=(5, 5))
plot_umap(adata, 'level_1_2', 'Cambi_051724_QC_03_int_level_1_2')

## remove flag
adata.obs['cl_flag'] = adata.obs.apply(lambda row: 'singlet' if row['level_1'] in ['18', '22', '28', '30'] and row['cl_flag'] == 'doublet' else row['cl_flag'], axis=1)
adata.obs["cl_flag"].value_counts()

sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=4, figsize=(5, 5))
sc.pl.umap(adata, color='cl_flag', save= '_Cambi_051724_QC_03_int_cl_flag.pdf')

## flag additional low quality clusters
adata.obs['cl_flag'] = adata.obs.apply(lambda row: 'low_quality' if row['level_2'] in ['7_13', '11_13', '12_12', '14_15', '20_8', '21_9', '24_7',] else row['cl_flag'], axis=1)
adata.obs["cl_flag"].value_counts()

sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=4, figsize=(5, 5))
sc.pl.umap(adata, color='cl_flag', save= '_Cambi_051724_QC_03_int_cl_flag.pdf')

## save data
adata.write_h5ad("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int.h5ad", compression= "gzip")
adata.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_obs.csv', index=True)

#####
#####

import os
import anndata as ad
import scanpy as sc
import pandas as pd
import numpy as np

adata = sc.read_h5ad('F:/Colab/Cambi_lab/Data/h5ad/Cambi_051724_dbl_filt.h5ad')
cl_QC_info = pd.read_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_obs.csv', index_col=0)
cl_QC_info = cl_QC_info[['cl_flag']].copy()
adata.obs = adata.obs.join(cl_QC_info)
adata.obs["cl_flag"].value_counts()

adata = adata[adata.obs["cl_flag"] == "singlet"].copy()
adata.write_h5ad("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_filt.h5ad", compression= "gzip")
adata.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_filt_obs.csv', index=True)

#####
#####
