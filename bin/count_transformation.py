#!/usr/bin/env python3

#0.4.1 QC for mitochondrial contamination, gene counts, cell number and Co
#kernel scVelo
import scanpy as sc
import os
import fire
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import seaborn as sns
import math
from scipy.stats import median_abs_deviation


def count_trans(dimensionality_out: str):
    adata= sc.read(dimensionality_out)
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    sc.pp.highly_variable_genes(adata, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, batch_key='sanger_id')
    adata.var.highly_variable.value_counts()

    adata_normalised = adata.copy() ##rather than saving and re-reeading in the same file! made a copy

    sc.pp.scale(adata_normalised, max_value=10)

    sc.tl.pca(adata_normalised, svd_solver='arpack')
    sc.pl.pca(adata_normalised, color='CST3', save="_plot24_normalised.png", show=False)#plot24
    #plot25
    sc.pl.umap(adata_normalised, color=['pcw', 'sanger_id', 'gender', 'sorting', 'chemistry_sorting'], save="_metadata.png", show=False)
    #plot26
    sc.pl.violin(adata_normalised, ['n_genes_by_counts'],rotation=90, groupby="sanger_id", save="_n_genes_by_counts_sangerid.png", show=False)
    #plot27
    sc.pl.violin(adata_normalised, ['total_counts'],rotation=90, groupby="sanger_id", save="_total_counts_sangerid.png", show=False)
    #plot28
    sc.pl.violin(adata_normalised, [ 'pct_counts_mt'],rotation=90,groupby="sanger_id", save="_pct_counts_mt_sangerid.png", show=False)

    sc.tl.umap(adata_normalised, min_dist=0.3, spread=1.0,init_pos='spectral') # Using the default distance scaling from umap-learn
    adata_normalised.obsm['X_umap_disp_03'] = adata_normalised.obsm['X_umap']

    adata.write("adata_transformation.h5ad")

if __name__ == "__main__":
    fire.Fire(count_trans)
