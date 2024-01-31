#!/usr/bin/env python3

#0.3.1 Scrublet
#Run Scrublet to identify doublets
#Kernel scVeloIsaac

import scanpy as sc
import os
#import anndata as ad
import anndata as ad
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

import seaborn as sns
import math
#from plotnine import *
import scrublet as scr
from scipy.stats import median_abs_deviation
#import sctk as sk

# import tables
# import scipy.sparse as sp
# from typing import Dict, Optional

# Mark doublets using Scrublet
# ~10000 cells per replicates so expected doublet rate is 7.6% 
# Needs to be done on a sample by sample basis

def calculate_expected_doublet_rate(adata_object):
    """
    for a given adata object, using the number of cells in the object to return 
    the expected doublet rate and the number of cells in the object
    """
    #expected values from https://uofuhealth.utah.edu/huntsman/shared-resources/gba/htg/single-cell/genomics-10x
    expected_rates = {1000: 0.008, 
                      2000: 0.016,
                    3000: 0.023,
                     4000: 0.031,
                     5000: 0.039,
                     6000: 0.046,
                     7000: 0.054,
                     8000: 0.061,
                     9000: 0.069,
                     10_000: 0.076}
    #number of cells (rounded)
    real_cells = math.ceil(adata_object.shape[0] / 1000) * 1000
    if real_cells > 10_000:
        real_cells = 10_000
            # Handle the case for 0 real_cells
    expected_rate = expected_rates[real_cells]
    return expected_rate

adata = sc.read('/sc-analysis-nf/output/FCAImmP7241240/Fskin_obj_0_2_1.h5ad')

adata.obs['sanger_id'].unique()

####would be good to hav this as a figure in html page 
calculate_expected_doublet_rate(adata)

sample_list = ['FCAImmP7241240']
for param in ['sanger_id']:
    print(param)
    print(adata.obs[param].value_counts(), '\n')

for sample in sample_list:
    adata_byid = adata[adata.obs.sanger_id == sample]
    real_cells = math.ceil(adata_byid.shape[0] / 1000) * 1000
    print(f"Sample {sample}, real_cells: {real_cells}")


#### this is and issue (seemingly only worls at the start of a cell- is there alternative?)
#%%capture
samples_scrubbed = {}
scrubs = {}
for sample in sample_list:
    adata_byid = adata[adata.obs.sanger_id == sample]
    scrub = scr.Scrublet(adata_byid.X, expected_doublet_rate = calculate_expected_doublet_rate(adata_byid))
    adata_byid.obs['doublet_scores'], adata_byid.obs['predicted_doublets'] = scrub.scrub_doublets(min_counts=2, min_cells=3, min_gene_variability_pctl=85, n_prin_comps=30)
    samples_scrubbed[sample] = adata_byid
    scrubs[sample] = scrub
    
#  Run per lane, each sample as identified by sanger_id was sequenced on separate lane here

#%matplotlib inline
for sample in scrubs:
    scrubs[sample].plot_histogram()
plt.savefig("scrublet_histogram.png")
#%matplotlib inline
for sample in scrubs:
    scrubs[sample].set_embedding('UMAP', scr.get_umap(scrub.manifold_obs_, 10, min_dist=0.3))

#%matplotlib inline
for sample in scrubs:
    print(scrubs[sample])
    scrubs[sample].plot_embedding('UMAP', order_points=True);
plt.savefig("scrublet_umap.png")

adata = ad.concat(samples_scrubbed, index_unique= None)

adata.obs.sanger_id.unique()
adata.write('Fskin_obj_0_3_1.h5ad')

###########################################
# 3 total outputs (2 plots & 1 h5ad file) #
###########################################