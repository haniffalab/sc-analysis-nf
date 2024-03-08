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

def dimensionality(dimensionality_in: str, dimensionality_out: str):

    adata= sc.read(dimensionality_in)

    unique_donors = list(adata.obs['donor'].unique())
    print(unique_donors)

    for gender in ['male', 'female']:
        gender_data = adata[adata.obs['gender'] == gender]
        unique_donors = gender_data.obs['donor'].nunique()
        unique_pcw = gender_data.obs['pcw'].unique()

        print(f"Gender: {gender}")
        print(f"Number of Unique Donors: {unique_donors}\n")
        print(f"Unique pcw values: {unique_pcw}\n")
############### Principal componant analysis (PCA)
    sc.pp.scale(adata, max_value=10)
    sc.tl.pca(adata, svd_solver='arpack')
    sc.pl.pca_scatter(adata, color='KRT4', save="_FCAImmP7241240.png", show=False) #plot 8
    sc.pl.pca_variance_ratio(adata, log=True, save="_FCAImmP7241240.png", show=False) #plot 9 
##############nearest neighbour graph
    sc.pp.neighbors(adata, n_neighbors=15, n_pcs=50)
    sc.tl.umap(adata)
    ##plot 10
    sc.pl.umap(adata, color=['LYZ','MS4A1','CD3D'], save="_neighbour_FCAImmP7241240.png", show=False) #color need to be custom?
    ## leiden is a clustering algorithm- Simone mentions it in clustering and annotation section. 
    sc.tl.leiden(adata) #- 8. clustering 
   
    # This information can be used to identify & characterize cells at different stages of the cell cycle, 
    # providing insights into cell proliferation, differentiation, & other biological processes.
    s_genes = ['MCM5','PCNA','TYMS','FEN1','MCM7','MCM4','RRM1','UNG','GINS2','MCM6','CDCA7','DTL','PRIM1',
            'UHRF1','CENPU','HELLS','RFC2','POLR1B','NASP','RAD51AP1','GMNN','WDR76','SLBP','CCNE2','UBR7',
            'POLD3','MSH2','ATAD2','RAD51','RRM2','CDC45','CDC6','EXO1','TIPIN','DSCC1','BLM','CASP8AP2',
            'USP1','CLSPN','POLA1','CHAF1B','MRPL36','E2F8']
    g2m_genes = ['HMGB2','CDK1','NUSAP1','UBE2C','BIRC5','TPX2','TOP2A','NDC80','CKS2','NUF2','CKS1B',
                'MKI67','TMPO','CENPF','TACC3','PIMREG','SMC4','CCNB2','CKAP2L','CKAP2','AURKB','BUB1',
                'KIF11','ANP32E','TUBB4B','GTSE1','KIF20B','HJURP','CDCA3','JPT1','CDC20','TTK','CDC25C',
                'KIF2C','RANGAP1','NCAPD2','DLGAP5','CDCA2','CDCA8','ECT2','KIF23','HMMR','AURKA','PSRC1',
                'ANLN','LBR','CKAP5','CENPE','CTCF','NEK2','G2E3','GAS2L3','CBX5','CENPA']
    cell_cycle_genes = s_genes + g2m_genes

    adata= adata.raw.to_adata()
    cell_cycle_genes = [x for x in cell_cycle_genes if x in adata.var_names]
    len(cell_cycle_genes)

    sc.tl.score_genes_cell_cycle(adata, s_genes=s_genes, g2m_genes=g2m_genes)
    #plot11
    sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_doubletscores.png", show=False)

    adata = adata[adata.obs['QC'] == 'Pass']
    adata.layers["counts"] = adata.X.copy()
    #plot12 
    sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_qc_doubletscores.png", show=False)
    #plot13
    sc.pl.umap(adata, color=['leiden'], save="_leiden.png", show=False) #UserWarning: No data for colormapping provided via 'c'. Parameters 'cmap' will be ignored
    #plot14
    sc.pl.umap(adata, color=['pct_counts_ribo','pct_counts_mt'], save="_pct_counts.png", show=False)
    #plot15
    with plt.rc_context({'figure.figsize': (5.5, 4)}): 
        sc.pl.violin(adata, ['pct_counts_ribo', 'pct_counts_mt'], groupby='leiden', stripplot = False, save="_leiden_pct_counts.png", show=False) #userwarning
    #plot16 
    sc.pl.umap(adata, color=['S_score','G2M_score','leiden'], save="_cellcycle.png", show=False)
    #plot17
    with plt.rc_context({'figure.figsize': (5.5, 4)}):
        sc.pl.violin(adata, ['S_score','G2M_score'], groupby = 'leiden', stripplot = False, save="_cellcycle.png", show=False)
    #plot18
    with plt.rc_context({'figure.figsize': (5.5, 4)}):
        sc.pl.violin(adata, ['outlier','mt_outlier'], groupby = 'leiden', stripplot = False, save="_outlier.png", show=False)
    #plot19
    sc.pl.umap(adata, color=['doublet_scores','pct_counts_mt','n_genes_by_counts'], save="_plot19.png", show=False)
    #plot20
    sc.pl.umap(adata, color=['pct_counts_ribo','pct_counts_mt','leiden'], save="_plot20_leiden.png", show=False)
    #plot21
    with plt.rc_context({'figure.figsize': (5.5, 4)}):
        sc.pl.violin(adata, ['pct_counts_ribo', 'pct_counts_mt'], groupby='leiden', stripplot = False, save="_pctcounts_leiden.png", show=False )
    #plot22
    with plt.rc_context({'figure.figsize': (5.5, 4)}):
        sc.pl.violin(adata, ['outlier','mt_outlier'], groupby = 'leiden', stripplot = False, save="_outlier_leiden.png", show=False)
    #plot23
    sc.pl.umap(adata, color=['S_score','G2M_score','leiden'], save="_cellcycle_leiden.png", show=False)
    #adata.raw = adata

    adata.write(dimensionality_out)


if __name__ == "__main__":
    fire.Fire(dimensionality)