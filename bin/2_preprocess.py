#!/usr/bin/env python3

import fire
import scanpy as sc
import pandas as pd
import numpy as np


def preprocess(anndata_path: str, results_file: str):
    adata = sc.read_h5ad(anndata_path)
    print(adata)
    # Show those genes that yield the highest fraction of counts in each single cell, across all cells.
    ####### eventually want to add wildcard of JOBID or something like that #######
    sc.pl.highest_expr_genes(adata, n_top=20, save="_2.png", show=False)

    # Basic filtering:
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=3)

    ###Let's assemble some information about mitochondrial genes, which are important for quality control.###
    # High proportions are indicative of poor-quality cells (Islam et al. 2014; Ilicic et al. 2016), possibly because of loss of cytoplasmic RNA from perforated cells.
    # The reasoning is that mitochondria are larger than individual transcript molecules and less likely to escape through tears in the cell membrane.
    adata.var["mt"] = adata.var_names.str.startswith(
        "MT-"
    )  # annotate the group of mitochondrial genes as 'mt'
    sc.pp.calculate_qc_metrics(
        adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True
    )

    # have to make 3 separate vioin plots as 'multi_panel=True' causing attribute error..
    # first violin plot
    sc.pl.violin(
        adata,
        "n_genes_by_counts",
        jitter=0.4,
        save="_n_genes_by_counts_2.png",
        show=False,
    )
    # second violin plot
    sc.pl.violin(
        adata, "total_counts", jitter=0.4, save="_total_counts_2.png", show=False
    )
    # third violin plot
    sc.pl.violin(
        adata, "pct_counts_mt", jitter=0.4, save="_pct_counts_mt_2.png", show=False
    )

    # Remove cells that have too many mitochondrial genes expressed or too many total counts:
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="pct_counts_mt",
        save="_pct_counts_mt_2.png",
        show=False,
    )
    sc.pl.scatter(
        adata,
        x="total_counts",
        y="n_genes_by_counts",
        save="_n_genes_by_counts_2.png",
        show=False,
    )

    # Actually do the filtering by slicing the AnnData object.
    adata = adata[adata.obs.n_genes_by_counts < 2500, :]
    adata = adata[adata.obs.pct_counts_mt < 5, :]

    # Total-count normalize (library-size correct) the data matrix X to 10,000 reads per cell, so that counts become comparable among cells.
    sc.pp.normalize_total(adata, target_sum=1e4)

    # Logarithmize the data:
    sc.pp.log1p(adata)

    # Identify highly-variable genes.
    sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5)

    sc.pl.highly_variable_genes(adata, save="scatter_2.png", show=False)

    # Set the '.raw' attribute of the AnnData object to the normalized and logarithmized raw gene expression for later use in differential testing and visualizations of gene expression.
    # This simply freezes the state of the AnnData object.
    adata.raw = adata

    ####################### commenting on whether some things are relevant ##########################################
    # If you don’t proceed below with correcting the data with sc.pp.regress_out and scaling it via sc.pp.scale,
    # you can also get away without using .raw at all.

    # The result of the previous highly-variable-genes detection is stored as an annotation in .var.highly_variable
    # and auto-detected by PCA and hence, sc.pp.neighbors and subsequent manifold/graph tools.
    # In that case, the step actually do the filtering below is unnecessary, too.
    #################################################################################################################

    # Actually do the filtering:
    adata = adata[:, adata.var.highly_variable]

    # Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
    # Scale the data to unit variance.
    sc.pp.regress_out(adata, ["total_counts", "pct_counts_mt"])

    # Scale each gene to unit variance. Clip values exceeding standard deviation 10.
    sc.pp.scale(adata, max_value=10)

    ## make this adata file into output:
    # To allow for input of the processed adata object to be utilised into further scripts.
    adata.write_h5ad(results_file)


if __name__ == "__main__":
    fire.Fire(preprocess)
