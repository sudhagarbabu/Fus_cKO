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

adata = sc.read_h5ad('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int.h5ad')
adata = adata[adata.obs["cl_flag"] == "singlet"].copy()
adata

## MCT anno
MCT_anno = pd.read_csv("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_level_2_cluster_stats_combined_with_MCT_anno.csv")
MCT_anno = MCT_anno[['level_2', 'MCT']].copy()

adata.obs['level_2'] = adata.obs['level_2'].astype(str)
MCT_anno['level_2'] = MCT_anno['level_2'].astype(str)
adata.obs['MCT'] = adata.obs['level_2'].map(dict(zip(MCT_anno['level_2'], MCT_anno['MCT'])))

## check
adata.obs['MCT'].isna().sum()

#####
#####

## Vis
import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
import IPython
IPython.display.set_matplotlib_formats = set_matplotlib_formats

os.chdir('F:/Colab/Cambi_lab/Integration/QC_03')
sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=5, figsize=(5, 5))

def plot_umap(adata, obs_key, filename_suffix):
    n_cls = len(adata.obs[obs_key].unique())
    colormap = plt.get_cmap('tab20', n_cls)
    colors = [colormap(i) for i in range(n_cls)]
    sc.pl.umap(adata, color=obs_key, legend_loc='on data', palette=colors, save=f'_{filename_suffix}.pdf')

plot_umap(adata, 'MCT', 'Cambi_051724_QC_03_int_MCT')

#####
#####

## refine anno based on supertype_name
adata.obs['level_1'] = adata.obs['level_1'].astype(str)
adata.obs['supertype_name'] = adata.obs['supertype_name'].astype(str)
adata.obs['MCT'] = adata.obs['MCT'].astype(str)

level_1_mappings = {
    ('18', '22'): {'OPC': 'OPC', 'COP': 'COP', 'NFOL': 'NFOL', 'MFOL': 'MFOL', 'MOL': 'MOL'},
    ('13',): {'Peri': 'Peri', 'SMC': 'SMC'}
}

for levels, mapping in level_1_mappings.items():
    mask = adata.obs['level_1'].isin(levels)
    for key, val in mapping.items():
        adata.obs.loc[mask & adata.obs['supertype_name'].str.contains(key, case=False, na=False), 'MCT'] = val

plot_umap(adata, 'MCT', 'Cambi_051724_QC_03_int_MCT')

## save data
adata.write_h5ad("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_MCT_anno.h5ad", compression= "gzip")
adata.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_MCT_anno_obs.csv', index=True)

#####
#####

## merge_clusters
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

script_dirs = ['E:/multiomatic/code/', 'E:/transcriptomic_clustering']
for d in script_dirs:
    path = os.path.abspath(d)
    if path not in sys.path:
        sys.path.append(path)

from Pycode_iter_clust_scvi import iter_clust_scvi
from Pycode_scanpy_helper_funtions import downsample_by_clusters
import transcriptomic_clustering as tc
from transcriptomic_clustering.merging import merge_clusters, get_cluster_assignments
from transcriptomic_clustering import hclust, get_cluster_means

adata = sc.read_h5ad('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_MCT_anno.h5ad')
adata = adata[adata.obs["cl_flag"] == "singlet"].copy()
adata

import matplotlib.pyplot as plt
from matplotlib_inline.backend_inline import set_matplotlib_formats
import IPython
IPython.display.set_matplotlib_formats = set_matplotlib_formats

os.chdir('F:/Colab/Cambi_lab/Integration/QC_03')
sc.set_figure_params(dpi=300, dpi_save=300, frameon=False, fontsize=5, figsize=(5, 5))

def plot_umap(adata, obs_key, filename_suffix):
    n_cls = len(adata.obs[obs_key].unique())
    colormap = plt.get_cmap('tab20', n_cls)
    colors = [colormap(i) for i in range(n_cls)]
    sc.pl.umap(adata, color=obs_key, legend_loc='on data', palette=colors, save=f'_{filename_suffix}.pdf')

#####
#####

## Microglia
adata_mic = adata[adata.obs["MCT"] == "Microglia"].copy()

cluster_key = "level_2"
adata_reduced = ad.AnnData(X=adata_mic.obsm["X_scVI"], obs=adata_mic.obs.copy())
adata_mic.obs[cluster_key] = adata_mic.obs[cluster_key].astype("category")
cluster_assignments = get_cluster_assignments(adata_mic, cluster_key)
cluster_by_obs = adata_mic.obs[cluster_key].cat.codes.to_numpy()

thresh = {
    'q1_thresh': 0.50,
    'q2_thresh': None,
    'cluster_size_thresh': 10,
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.01,
    'lfc_thresh': 1.0,
    'score_thresh': 150,
    'low_thresh': 1,
    'min_genes': 1
}

merged_clusters, markers = merge_clusters(
    adata_norm=adata_mic,
    adata_reduced=adata_reduced,
    cluster_assignments=cluster_assignments,
    cluster_by_obs=cluster_by_obs,
    thresholds=thresh,
    de_method="ebayes",
    n_markers=0
)

cell_to_merged_clusters = {
    adata_mic.obs_names[i]: merged_cl
    for merged_cl, idxs in merged_clusters.items()
    for i in idxs
}

adata_mic.obs[f"{cluster_key}_merged"] = adata_mic.obs_names.map(cell_to_merged_clusters)
adata_mic.obs[f"{cluster_key}_merged"] = adata_mic.obs[f"{cluster_key}_merged"].astype("category")
plot_umap(adata_mic, 'level_2_merged', 'Cambi_051724_QC_03_int_level_2_merged_Microglia')

#####

cluster_key = "level_2_merged"
adata_reduced = ad.AnnData(X=adata_mic.obsm["X_scVI"], obs=adata_mic.obs.copy())
adata_mic.obs[cluster_key] = adata_mic.obs[cluster_key].astype("category")
cluster_assignments = get_cluster_assignments(adata_mic, cluster_key)
cluster_by_obs = adata_mic.obs[cluster_key].cat.codes.to_numpy()

thresh = {
    'q1_thresh': 0.50,
    'q2_thresh': None,
    'cluster_size_thresh': 100,
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.01,
    'lfc_thresh': 1.0,
    'score_thresh': 200,
    'low_thresh': 1,
    'min_genes': 0
}

merged_clusters, markers = merge_clusters(
    adata_norm=adata_mic,
    adata_reduced=adata_reduced,
    cluster_assignments=cluster_assignments,
    cluster_by_obs=cluster_by_obs,
    thresholds=thresh,
    de_method="ebayes",
    n_markers=0
)

cell_to_merged_clusters = {
    adata_mic.obs_names[i]: merged_cl
    for merged_cl, idxs in merged_clusters.items()
    for i in idxs
}

adata_mic.obs[f"{cluster_key}_final"] = adata_mic.obs_names.map(cell_to_merged_clusters)
adata_mic.obs[f"{cluster_key}_final"] = adata_mic.obs[f"{cluster_key}_final"].astype("category")
#adata_mic.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_Microglia_obs.csv', index=True)
plot_umap(adata_mic, 'level_2_merged_final', 'Cambi_051724_QC_03_int_level_2_merged_final_Microglia')

adata_mic_clusters = sorted(adata_mic.obs['level_2_merged_final'].unique(), key=lambda x: tuple(map(int, x.split("_"))))
adata_mic_map = {cl: f"Microglia_{i+1}" for i, cl in enumerate(adata_mic_clusters)}
adata_mic.obs["HCT"] = adata_mic.obs['level_2_merged_final'].map(adata_mic_map).astype("category")
plot_umap(adata_mic, 'HCT', 'Cambi_051724_QC_03_int_HCT_Microglia')

#####
#####

## MOL
adata_mol = adata[adata.obs["MCT"] == "MOL"].copy()

cluster_key = "level_2"
adata_reduced = ad.AnnData(X=adata_mol.obsm["X_scVI"], obs=adata_mol.obs.copy())
adata_mol.obs[cluster_key] = adata_mol.obs[cluster_key].astype("category")
cluster_assignments = get_cluster_assignments(adata_mol, cluster_key)
cluster_by_obs = adata_mol.obs[cluster_key].cat.codes.to_numpy()

thresh = {
    'q1_thresh': 0.50,
    'q2_thresh': None,
    'cluster_size_thresh': 10,
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.01,
    'lfc_thresh': 1.0,
    'score_thresh': 200,
    'low_thresh': 1,
    'min_genes': 3
}

merged_clusters, markers = merge_clusters(
    adata_norm=adata_mol,
    adata_reduced=adata_reduced,
    cluster_assignments=cluster_assignments,
    cluster_by_obs=cluster_by_obs,
    thresholds=thresh,
    de_method="ebayes",
    n_markers=0
)

cell_to_merged_clusters = {
    adata_mol.obs_names[i]: merged_cl
    for merged_cl, idxs in merged_clusters.items()
    for i in idxs
}

adata_mol.obs[f"{cluster_key}_merged"] = adata_mol.obs_names.map(cell_to_merged_clusters)
adata_mol.obs[f"{cluster_key}_merged"] = adata_mol.obs[f"{cluster_key}_merged"].astype("category")
plot_umap(adata_mol, 'level_2_merged', 'Cambi_051724_QC_03_int_level_2_merged_MOL')

#####

cluster_key = "level_2_merged"
adata_reduced = ad.AnnData(X=adata_mol.obsm["X_scVI"], obs=adata_mol.obs.copy())
adata_mol.obs[cluster_key] = adata_mol.obs[cluster_key].astype("category")
cluster_assignments = get_cluster_assignments(adata_mol, cluster_key)
cluster_by_obs = adata_mol.obs[cluster_key].cat.codes.to_numpy()

thresh = {
    'q1_thresh': 0.50,
    'q2_thresh': None,
    'cluster_size_thresh': 100,
    'qdiff_thresh': 0.7,
    'padj_thresh': 0.01,
    'lfc_thresh': 1.0,
    'score_thresh': 200,
    'low_thresh': 1,
    'min_genes': 0
}

merged_clusters, markers = merge_clusters(
    adata_norm=adata_mol,
    adata_reduced=adata_reduced,
    cluster_assignments=cluster_assignments,
    cluster_by_obs=cluster_by_obs,
    thresholds=thresh,
    de_method="ebayes",
    n_markers=0
)

cell_to_merged_clusters = {
    adata_mol.obs_names[i]: merged_cl
    for merged_cl, idxs in merged_clusters.items()
    for i in idxs
}

adata_mol.obs[f"{cluster_key}_final"] = adata_mol.obs_names.map(cell_to_merged_clusters)
adata_mol.obs[f"{cluster_key}_final"] = adata_mol.obs[f"{cluster_key}_final"].astype("category")
#adata_mol.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_MOL_obs.csv', index=True)
plot_umap(adata_mol, 'level_2_merged_final', 'Cambi_051724_QC_03_int_level_2_merged_final_MOL')

adata_mol_clusters = sorted(adata_mol.obs['level_2_merged_final'].unique(), key=lambda x: tuple(map(int, x.split("_"))))
adata_mol_map = {cl: f"MOL_{i+1}" for i, cl in enumerate(adata_mol_clusters)}
adata_mol.obs["HCT"] = adata_mol.obs['level_2_merged_final'].map(adata_mol_map).astype("category")
plot_umap(adata_mol, 'HCT', 'Cambi_051724_QC_03_int_HCT_MOL')

#####
#####

## HCT anno
adata.obs["HCT"] = pd.NA
adata.obs.loc[adata_mic.obs_names, "HCT"] = adata_mic.obs["HCT"]
adata.obs.loc[adata_mol.obs_names, "HCT"] = adata_mol.obs["HCT"]
adata.obs["HCT"] = adata.obs["HCT"].fillna(adata.obs["MCT"]).astype("category")

plot_umap(adata, 'HCT', 'Cambi_051724_QC_03_int_HCT')

## save data
adata.write_h5ad("F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_HCT_anno.h5ad", compression= "gzip")
adata.obs.to_csv('F:/Colab/Cambi_lab/Integration/QC_03/Cambi_051724_QC_03_int_HCT_anno_obs.csv', index=True)
