import scprinter as scp
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import time
import pandas as pd
import numpy as np
import os
import pickle
import torch
import matplotlib as mpl
mpl.rcParams['pdf.fonttype'] = 42
from scanpy.plotting.palettes import zeileis_28
from tqdm.contrib.concurrent import *
from tqdm.auto import *
import anndata
import scanpy as sc
import statistics as stat
import json
import csv
import re
import copy
from sklearn.preprocessing import OneHotEncoder


# Specify directories we will use. Make sure you provide the full absolute path and not the relative path
main_dir = '/mnt/disks/data/hipsc-ec/broad-nnfc-hiPSC-to-EC/results/'
work_dir = f'{main_dir}/2026-01-23_scprinter/'
if not os.path.exists(work_dir):
    os.system("mkdir -p " + work_dir)
frag_dir = f'{work_dir}/fragment_file_merged'
if not os.path.exists(frag_dir):
    os.system("mkdir -p " + frag_dir)

import snapatac2 as snap
fragment_file = f'{frag_dir}/merged_fragments.namesort.bed.gz'
genome = scp.genome.hg38

# Initialize the scPrinter object
# When you finish using the object, run printer.close() otherwise you won't be able to load it properly next time.
import time
start = time.time()
printer_path = os.path.join(work_dir, 'hipsc_ec_scATAC_scprinter.h5ad')
if os.path.exists(printer_path):
    printer = scp.load_printer(f'{work_dir}/hipsc_ec_scATAC_scprinter.h5ad', genome)
else:
    # Load data from the fragments file. Needs at least theses four columns: chromosome, start, end, cell barcode
    printer = scp.pp.import_fragments(
                            path_to_frags=f'{fragment_file}',
                            barcodes=None,
                            savename=printer_path,
                            genome=genome,
                            sorted_by_barcode=True,
                            low_memory=False,
                            )

scp.pp.call_peaks(printer=printer,
                  frag_file=fragment_file,
                  cell_grouping=[None], # here we call peaks on the cells that are included in the final analyses
                  group_names=['all'],
                  preset='seq2PRINT', n_jobs=1)
# Fetched the cleaned peaks, save, it will be used in the next step
cleaned_peaks = pd.DataFrame(printer.uns["peak_calling"]['all_cleaned'][:])
cleaned_peaks.to_csv(f'{work_dir}/seq2print_cleaned_narrowPeak.bed',
                     sep='\t', header=False, index=False)
    
# Call peaks using chromvar preset, this set of peak are recommended to be use as cell x peak for scATAC-seq data, or analysis
scp.pp.call_peaks(printer=printer,
                  frag_file=fragment_file,
                  cell_grouping=[None], # here we call peaks on the cells that are included in the final analyses
                  group_names=['chromvar_all'],
                  preset='chromvar',
                  overwrite=False)

# Fetched the cleaned peaks, save, it will be used in the next step
cleaned_peaks = pd.DataFrame(printer.uns["peak_calling"]['chromvar_all_cleaned'][:])
cleaned_peaks.to_csv(f'{work_dir}/regions.bed',
                     sep='\t', header=False, index=False)

# First construct a peak-by-cell matrix of ATAC counts
peak_path = f'{work_dir}/regions.bed'
adata = scp.pp.make_peak_matrix(printer,
                       regions=peak_path,
                       region_width=300,
                       cell_grouping=None,
                       group_names=None,
                       sparse=True)

adata.write(f'{work_dir}/cell_peak.h5ad')


# Only keep peaks with > 0 coverage
adata = anndata.read_h5ad(f'{work_dir}/cell_peak.h5ad')
coverage = adata.X.sum(axis=0)
adata = adata[:, coverage > 0]


# We can calculate chromVAR motif scores using either GPU (device = "cuda", much faster) or CPU (device = "cpu", slower)
device = "cpu"

if device == "cuda":
    import warnings
    warnings.filterwarnings("ignore")
    import scanpy as sc
    import anndata
    import cupy as cp
    import cupyx as cpx
    import time
    import rmm
    from rmm.allocators.cupy import rmm_cupy_allocator
    rmm.reinitialize(
        managed_memory=True, # Allows oversubscription
        pool_allocator=True, # default is False
        devices=0, # GPU device IDs to register. By default registers only GPU 0.
    )
    cp.cuda.set_allocator(rmm_cupy_allocator)

# Sample background peaks for each peak
scp.chromvar.sample_bg_peaks(adata,
                             genome=genome,
                             method='chromvar',
                             niterations=250)

# Scan motifs
motif = scp.motifs.FigR_Human_Motifs(genome,
                                     bg=list(adata.uns['bg_freq']),
                                     n_jobs=100,
                                     pvalue=5e-5, mode='motifmatchr')
motif.prep_scanner(None, pvalue=5e-5)
motif.chromvar_scan(adata)

# Compute motif scores for single cells
chromvar = scp.chromvar.compute_deviations(adata, chunk_size=50000, device=device)


# Save for later use
chromvar.write(f'{work_dir}/chromvar_cisbp.h5ad')

chromvar = anndata.read_h5ad(f'{work_dir}/chromvar_cisbp.h5ad')
print(chromvar.X.shape)

# First perform PCA to get PC loadings
sc.tl.pca(chromvar)

# Using PCA embedding for UMAP projection
from cuml import UMAP
vec = UMAP(metric='cosine').fit_transform(chromvar.obsm["X_pca"])
chromvar.obsm['X_umap'] = vec

# Visualize clusters and cell-type-specific TF motif scores
sc.pl.umap(chromvar, color=['CEBPA', 'TCF7', "PAX5"],cmap='RdBu_r',vmin=-3, vmax=3)

chromvar


sc.pp.neighbors(chromvar)
sc.tl.leiden(chromvar, flavor="igraph", n_iterations=10, resolution=0.2, random_state=2)
sc.pl.umap(chromvar, color="leiden", size=2)


chromvar


mapping = pd.read_csv(f'{work_dir}/barcode_to_celltype.csv')
obs = chromvar.obs.to_pandas()
obs["celltype"] = obs.index.map(
    dict(zip(mapping["barcode"], mapping["cell_type"]))
)
chromvar.obs = obs

cell_depth = np.array(np.sum(adata.X, axis=1)).squeeze()
threshold = 5e6 # 5M totdal depth per pseudobulk. It is recommended to use the same depth for pseudobulks to prevent global biases.
pbulk_centers = []
barcode_groups = []
cell_barcodes = np.array(chromvar.obs.index)

# Fit a KNN on our cells
from sklearn.neighbors import NearestNeighbors
nbrs = NearestNeighbors(n_neighbors=chromvar.shape[0], algorithm='ball_tree').fit(chromvar.obsm["X_pca"])

# Go through each cell type
for celltype in np.unique(chromvar.obs["celltype"]):
    celltype_center_count = 0
    cell_inds = np.where(chromvar.obs["celltype"] == celltype)[0]

    # Keep sampling new pseudobulk centers untill we have 10 centers per cell type
    while celltype_center_count < 10:

        # Generate a new pseudobulk center by random sampling
        new_center = np.random.choice(cell_inds, 1)[0]

        # Sort cells by their distance to the center cell
        distances, indices = nbrs.kneighbors([chromvar.obsm["X_pca"][new_center, :]])

        # Keep adding the next nearest neighbor from the center cell to the pseudo bulk until we reach a depth threshold
        nbr_inds = indices[0, :]
        cumulative_depth = np.cumsum(cell_depth[nbr_inds])
        n_members = np.min(np.where(cumulative_depth > threshold)[0])
        pbulk_members = nbr_inds[:n_members]

        # Calculate purity of pseudobulk: percentage of the dominant cell type
        celltype_labels = chromvar.obs["celltype"][pbulk_members]
        purity = sum(celltype_labels == stat.mode(celltype_labels))/len(celltype_labels)

        # Discard the current pseudobulk if it's impure (mixture of different cell types)
        if purity < 0.99:
            continue

        # Append the new pseudobulk to the list
        pbulk_centers.append(new_center)
        celltype_center_count += 1
        center_ind = len(pbulk_centers)
        chromvar.obs[f"pbulk_{center_ind}_member"] = np.zeros(chromvar.shape[0])
        chromvar.obs[f"pbulk_{center_ind}_member"][pbulk_members] = 1
        new_bc_group = pd.DataFrame({"barcode":cell_barcodes[pbulk_members], "group":f"{celltype}_pbulk_{celltype_center_count}"})
        barcode_groups.append(new_bc_group)
barcode_groups = pd.concat(barcode_groups, axis=0)
np.array(pbulk_centers)

barcode_groups

# Plot pseudo-bulk centers on UMAP (black dots are pseudobulk centers)
chromvar.obs["pbulk_center"] = np.zeros(chromvar.shape[0])
chromvar.obs["pbulk_center"][pbulk_centers] = 1
import matplotlib.colors as mcolors
cmap = mcolors.ListedColormap(['#D3D3D3', 'black'])
sc.pl.umap(chromvar, color="pbulk_center", size=5, cmap=cmap)

sc.pl.umap(chromvar, color=["pbulk_5_member", "pbulk_15_member", "pbulk_25_member"], size=5, cmap=cmap)

printer.close()





