#!/usr/bin/env python3

# 0.3.1 Scrublet
# Run Scrublet to identify doublets
# Kernel scVeloIsaac

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
import scrublet as scr
from scipy.stats import median_abs_deviation

def calculate_expected_doublet_rate(adata_object):
    # expected values from https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x
    expected_rates = {
        1000: 0.008,
        2000: 0.016,
        3000: 0.023,
        4000: 0.031,
        5000: 0.039,
        6000: 0.046,
        7000: 0.054,
        8000: 0.061,
        9000: 0.069,
        10_000: 0.076,
    }
    # number of cells (rounded)
    real_cells = math.ceil(adata_object.shape[0] / 1000) * 1000
    if real_cells > 10_000:
        real_cells = 10_000
        # Handle the case for 0 real_cells
    expected_rate = expected_rates[real_cells]
    return expected_rate

def scrublet(samples: str):
    ##reading in data
    adata = sc.read("metadata_adata.h5ad")

    calculate_expected_doublet_rate(adata) # would be good to printed in html output page
    sample_list = [samples] # how is best to over come this? same in the add_metadata.py script

    for sample in sample_list: 
        adata_byid = adata[adata.obs.sanger_id == sample]
        real_cells = math.ceil(adata_byid.shape[0] / 1000) * 1000
    #### this is and issue (seemingly only worls at the start of a cell- is there alternative?)
    samples_scrubbed = {}
    scrubs = {}
    for sample in sample_list: #
        adata_byid = adata[adata.obs.sanger_id == sample]
        scrub = scr.Scrublet(
            adata_byid.X, expected_doublet_rate=calculate_expected_doublet_rate(adata_byid)
        )
        adata_byid.obs["doublet_scores"], adata_byid.obs["predicted_doublets"] = (
            scrub.scrub_doublets(
                min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30
            )
        )
        samples_scrubbed[sample] = adata_byid
        scrubs[sample] = scrub
    #  Run per lane, each sample as identified by sanger_id was sequenced on separate lane here
    ## Plots 
    for sample in scrubs:
        scrubs[sample].plot_histogram()
    plt.savefig("scrublet_histogram.png")

    for sample in scrubs:
        scrubs[sample].set_embedding(
            "UMAP", scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3)
        )
    ########################### plot above isnt saved? 
    for sample in scrubs:
        print(scrubs[sample])
        scrubs[sample].plot_embedding("UMAP", order_points=True)
    plt.savefig("scrublet_umap.png")

    adata = ad.concat(samples_scrubbed, index_unique=None)
    #writing out data 
    adata.write("scrublet_adata.h5ad")

if __name__ == "__main__":
    fire.Fire(scrublet)

###########################################
# 3 total outputs (2 plots & 1 h5ad file) #
###########################################
