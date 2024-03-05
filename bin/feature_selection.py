#!/usr/bin/env python3

#########  8. Feature selection  ###########

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

### add fire definition 
adata=sc.read("adata_filtered.h5ad")
sc.pl.highest_expr_genes(adata, n_top=20, save="_FCAImmP7241240.png", show=False ) #Change name 

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key='sanger_id')
sc.pl.highly_variable_genes(adata, save="_scatter_FCAImmP7241240.png", show=False) # plot 7 

adata.write("adata_feature_selection.h5ad")