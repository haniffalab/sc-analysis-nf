#!/usr/bin/env python3

######## Filtering low quality reads ########
#0.4.1 QC for mitochondrial contamination, gene counts, cell number and Co - kernel scVelo
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

def filter_reads(filter_in: str, samples: str, filter_out: str):
    adata= sc.read(filter_in)
    adata.raw = adata
    sample_list = [samples]

    adata.var["mt"] = adata.var_names.str.startswith("MT-") # mitochondrial genes
    adata.var["ribo"] = adata.var_names.str.startswith(("RPS", "RPL")) # ribosomal genes
    adata.var["hb"] = adata.var_names.str.contains(("^HB[^(P)]")) # hemoglobin genes.

    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True
    )
    p1 = sns.displot(adata.obs["total_counts"], bins=100, kde=False)
    p1.savefig("qcmetrics_totalcounts.png") #plot1
    sc.pl.violin(adata, "pct_counts_mt", save="_qcmetrics_pctcountsmt.png", show=False ) #plot2
    sc.pl.scatter(adata, "total_counts", "n_genes_by_counts", color="pct_counts_mt", save="_qcmetrics.png", show=False) #plot3
    
    # Automatic threshholding based on median absolute deviation (MAD)
    def is_outlier(adata, metric: str, nmads: int):
        M = adata.obs[metric]
        outlier = (M < np.median(M) - nmads * median_abs_deviation(M)) | (
            np.median(M) + nmads * median_abs_deviation(M) < M
        )
        return outlier

    samples_QC = {}
    QC = {}
    for sample in sample_list:
        adata_byid = adata[adata.obs.sanger_id == sample]
        adata_byid.obs["outlier"] = (
            is_outlier(adata_byid, "log1p_total_counts", 8)
            | is_outlier(adata_byid, "log1p_n_genes_by_counts", 8)
            )
        adata_byid.obs.outlier.value_counts()
        samples_QC[sample] = adata_byid
        QC[sample] = adata_byid
    ###### RuntimeWarning: Mean of empty slice.
    adata = ad.concat(samples_QC, index_unique= None)
    adata.obs.outlier.value_counts()

    samples_QC = {}
    Q = {}
    for sample in sample_list:
        adata_byid = adata[adata.obs.sanger_id == sample]
        adata_byid.obs["mt_outlier"] = is_outlier(adata_byid, "pct_counts_mt", 5) | (
            adata_byid.obs["pct_counts_mt"] > 30
            )
        adata_byid.obs.mt_outlier.value_counts()
        samples_QC[sample] = adata_byid
        QC[sample] = adata_byid
    #ImplicitModificationWarning: Trying to modify attribute `.obs` of view, initializing view as actual.
    adata = ad.concat(samples_QC, index_unique= None)
    adata.obs.mt_outlier.value_counts()

    ###### would be good to print this in the html page ######
    print(f"Total number of cells: {adata.n_obs}")
    adataQC = adata[(~adata.obs.outlier) & (~adata.obs.mt_outlier)].copy()
    print(f"Number of cells after filtering of low quality cells: {adataQC.n_obs}")

    #### plots4/5/6
    sc.pl.violin(adata, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True,
            groupby="sanger_id", save="_all_counts.png", show=False)
    sc.pl.violin(adataQC, ['n_genes_by_counts', 'total_counts', 'pct_counts_mt'], multi_panel=True,
            groupby="sanger_id", save="_qc_all_counts.png", show=False)

    conditions = [
        (adata.obs['predicted_doublets'] == True),
        (adata.obs['outlier'] == True),
            (adata.obs['mt_outlier'] == True),
        (adata.obs['n_genes_by_counts'] < 500),
    (adata.obs['n_genes_by_counts'] >= 500) & (adata.obs['predicted_doublets'] != True) & (adata.obs['mt_outlier'] != True) & (adata.obs['outlier'] != True) 
    ]

    values = ['Doublet', 'Count_outlier','MT_outlier','Low_nFeature', 'Pass']
    adata.obs['QC'] = np.select(conditions, values)
    adata.obs['QC'] = adata.obs['QC'].astype('category')
    adata.obs[['predicted_doublets','outlier', 'mt_outlier','n_genes_by_counts','QC']].head(10)
    
    adata.write(filter_out)

if __name__ == "__main__":
    fire.Fire(filter_reads)